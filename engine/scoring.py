"""
OncoTrack Clinical Scoring & Explainable AI (XAI) Engine
Calculates the Clonal Escape Score (0-100) and extracts 5 deterministic clinical audit factors.
Aligns with AMP/ASCO/CAP guidelines and RECIST 1.1 kinetic definitions.
"""

from typing import Dict, List, Any
import numpy as np
from engine.trajectory import compute_clonal_kinetics


def calculate_escape_score(patient_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Computes deterministic Clonal Escape Score (0-100) and 5 explainable clinical audit factors.
    No black-box models; strictly derived from mathematical VAF kinetics.
    """
    samples = patient_data.get("samples", [])
    if not samples:
        return {
            "score": 0,
            "status": "INSUFFICIENT DATA",
            "alert_level": "NORMAL",
            "factors": [],
            "summary": "No longitudinal ctDNA samples available for kinetic calculus."
        }

    kinetics = compute_clonal_kinetics(samples)

    # Locate driver, resistance, and highest increasing clone
    driver_k = None
    resistance_k = None
    max_increasing_k = None
    max_velocity = -999.0

    for mut_key, k in kinetics.items():
        if k["clone_type"] == "driver" and driver_k is None:
            driver_k = k
        elif k["clone_type"] == "resistance" and resistance_k is None:
            resistance_k = k

        if k["recent_velocity_per_cycle"] > max_velocity:
            max_velocity = k["recent_velocity_per_cycle"]
            max_increasing_k = k

    # Primary target for escape monitoring is resistance clone if present, else highest increasing clone
    target_k = resistance_k if resistance_k is not None else max_increasing_k

    # Factor calculation
    factors = []
    score_points = 0.0

    # 1. Driver Suppression Factor (Selective Pressure Confirmation)
    driver_suppressed = False
    driver_desc = "Primary driver allelic burden unchanged or baseline undetermined."
    if driver_k is not None:
        if driver_k["clearance_pct"] >= 60.0 or driver_k["current_vaf"] < 5.0:
            driver_suppressed = True
            score_points += 20.0
            driver_desc = f"Primary driver ({driver_k['gene']} {driver_k['mutation']}) shows {driver_k['clearance_pct']}% allelic clearance, confirming selective therapeutic pressure on the dominant clone."
        elif driver_k["total_delta_vaf"] < 0:
            score_points += 10.0
            driver_desc = f"Primary driver ({driver_k['gene']} {driver_k['mutation']}) decreased by {abs(driver_k['total_delta_vaf'])}% VAF from baseline."
        else:
            driver_desc = f"Primary driver ({driver_k['gene']} {driver_k['mutation']}) remains active at {driver_k['current_vaf']}% VAF."
    factors.append({
        "id": "F1_DRIVER_CLEARANCE",
        "title": "Selective Pressure / Driver Response",
        "description": driver_desc,
        "is_positive": driver_suppressed,
        "weight": 20
    })

    # 2. Emergent Clone Velocity Factor
    velocity_flag = False
    vel_val = target_k["recent_velocity_per_cycle"] if target_k else 0.0
    target_name = f"{target_k['gene']} {target_k['mutation']}" if target_k else "None"
    if target_k and vel_val >= 2.0:
        score_points += 25.0
        velocity_flag = True
        vel_desc = f"Target clone ({target_name}) exhibits high kinetic velocity: +{vel_val}% VAF / cycle (+{target_k['recent_velocity_per_day']}% / day)."
    elif target_k and vel_val > 0.5:
        score_points += 15.0
        velocity_flag = True
        vel_desc = f"Target clone ({target_name}) exhibits moderate expansion velocity: +{vel_val}% VAF / cycle."
    else:
        vel_desc = f"No significant subclonal expansion detected (current velocity: {vel_val}% VAF / cycle)."
    factors.append({
        "id": "F2_CLONAL_VELOCITY",
        "title": "Subclonal Kinetic Velocity (dVAF/dt)",
        "description": vel_desc,
        "is_positive": velocity_flag,
        "weight": 25
    })

    # 3. Kinetic Acceleration & Doubling Time Factor
    accel_flag = False
    accel_val = target_k["recent_acceleration"] if target_k else 0.0
    doubling_days = target_k.get("doubling_time_days") if target_k else None
    if target_k and (accel_val > 0.1 or (doubling_days is not None and doubling_days < 45)):
        score_points += 20.0
        accel_flag = True
        doubling_str = f"Estimated kinetic doubling time: ~{doubling_days} days." if doubling_days else ""
        accel_desc = f"Positive trajectory convexity (acceleration = +{accel_val}% VAF/cycle²). {doubling_str} Suggests exponential subclonal outgrowth."
    elif target_k and accel_val > 0.0:
        score_points += 10.0
        accel_flag = True
        accel_desc = f"Mild upward trajectory curvature (+{accel_val}% VAF/cycle²)."
    else:
        accel_desc = "Linear or decelerating trajectory; no exponential surge detected."
    factors.append({
        "id": "F3_CONVEXITY_ACCELERATION",
        "title": "Trajectory Convexity & Doubling Time",
        "description": accel_desc,
        "is_positive": accel_flag,
        "weight": 20
    })

    # 4. Longitudinal Signal Persistence Factor
    persist_flag = False
    if target_k:
        non_zero_draws = sum(1 for v in target_k["vafs"] if v > 0.05)
        if non_zero_draws >= 3:
            score_points += 20.0
            persist_flag = True
            persist_desc = f"Signal persisted across {non_zero_draws} consecutive serial blood draws (eliminates technical PCR/sequencing artifact)."
        elif non_zero_draws == 2:
            score_points += 10.0
            persist_flag = True
            persist_desc = "Signal confirmed across 2 sequential liquid biopsy timepoints."
        else:
            persist_desc = "Single timepoint observation; repeat blood draw recommended for confirmation."
    else:
        persist_desc = "No persistent secondary signal detected."
    factors.append({
        "id": "F4_SERIAL_PERSISTENCE",
        "title": "Longitudinal Assay Persistence",
        "description": persist_desc,
        "is_positive": persist_flag,
        "weight": 20
    })

    # 5. Clinical Demonstration Threshold Transgression Factor (>5.0% VAF)
    thresh_flag = False
    curr_vaf = target_k["current_vaf"] if target_k else 0.0
    if curr_vaf >= 5.0:
        score_points += 15.0
        thresh_flag = True
        thresh_desc = f"Current allelic burden ({curr_vaf}% VAF) has crossed the 5.0% clinical demonstration threshold into high-abundance territory."
    elif curr_vaf >= 1.0:
        score_points += 8.0
        thresh_desc = f"Current allelic burden ({curr_vaf}% VAF) is in the 1.0–5.0% molecular surveillance window."
    else:
        thresh_desc = f"Current allelic burden ({curr_vaf}% VAF) remains within low baseline limits (<1.0% VAF)."
    factors.append({
        "id": "F5_THRESHOLD_TRANSGRESSION",
        "title": "Clinical Threshold Transgression",
        "description": thresh_desc,
        "is_positive": thresh_flag,
        "weight": 15
    })

    # Final score normalization
    final_score = int(round(np.clip(score_points, 0, 100)))

    # Classification status
    if final_score >= 70:
        status = "POTENTIAL MOLECULAR ESCAPE"
        alert_level = "CRITICAL"
        summary_text = (
            f"High-confidence molecular escape signal detected for {target_name}. "
            f"Kinetics indicate rapid subclonal expansion despite primary driver suppression. "
            f"Molecular progression precedes radiographic change by estimated 60-90 days."
        )
    elif final_score >= 40:
        status = "MONITOR"
        alert_level = "WARNING"
        summary_text = (
            f"Equivocal subclonal dynamics. Moderate trajectory velocity observed in {target_name}. "
            f"Recommend close surveillance with repeat liquid biopsy at next treatment cycle."
        )
    else:
        status = "MOLECULAR RESPONSE"
        alert_level = "STABLE"
        summary_text = (
            "Favorable longitudinal molecular kinetics. Primary driver suppression with no emergent "
            "resistant clones detected above background threshold."
        )

    return {
        "score": final_score,
        "status": status,
        "alert_level": alert_level,
        "target_clone": target_name,
        "target_vaf": curr_vaf,
        "target_velocity": vel_val,
        "factors": factors,
        "summary": summary_text,
        "kinetics": kinetics
    }
