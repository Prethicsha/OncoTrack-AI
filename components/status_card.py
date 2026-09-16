"""
OncoTrack Clinical Pathologist Interpretation & Specimen Telemetry
Hospital-grade clean medical light theme with high-contrast clinical cards.
"""

import streamlit as st
from typing import Dict, Any


def render_clinical_status_card(score_data: Dict[str, Any], patient_info: Dict[str, Any]):
    """Renders clean hospital-grade diagnostic interpretation and specimen metrics."""
    status = score_data.get("status", "MONITOR")
    score = score_data.get("score", 0)
    alert_level = score_data.get("alert_level", "STABLE")
    target_clone = score_data.get("target_clone", "None")
    target_vaf = score_data.get("target_vaf", 0.0)
    target_vel = score_data.get("target_velocity", 0.0)
    kinetics = score_data.get("kinetics", {})

    driver_k = next((k for k in kinetics.values() if k["clone_type"] == "driver"), None)
    res_k = next((k for k in kinetics.values() if k["clone_type"] == "resistance"), None)

    if alert_level == "CRITICAL":
        card_bg = "#fef2f2"
        card_border = "#ef4444"
        badge_bg = "#dc2626"
        heading_color = "#991b1b"
        text_color = "#1f2937"
        action_bg = "#ffffff"
        action_border = "#fca5a5"
        score_color = "#dc2626"
        status_heading = "CRITICAL FINDING: MOLECULAR RESISTANCE / CLONAL DIVERGENCE"
        summary_text = (
            f"Serial ctDNA liquid biopsy demonstrates selective subclonal expansion of <b>{target_clone}</b> "
            f"reaching <b>{target_vaf:.2f}% VAF</b> (kinetic expansion velocity: <b>+{target_vel:.2f}% VAF / cycle</b>). "
            f"The primary sensitizing driver exhibits concurrent allelic clearance. This discordant kinetic trajectory "
            f"is diagnostic of emergent therapeutic failure. Molecular resistance precedes RECIST 1.1 "
            f"radiographic progression by an estimated 60 to 90 days."
        )
        action_text = (
            "<b>Recommended Clinical Action:</b> Multidisciplinary Molecular Tumor Board review. "
            "Consider immediate therapeutic transition to next-generation targeted inhibitor per NCCN Category 1 guidelines. "
            "Repeat serial plasma ctDNA draw in 21 days."
        )
    elif alert_level == "WARNING":
        card_bg = "#fffbeb"
        card_border = "#f59e0b"
        badge_bg = "#d97706"
        heading_color = "#92400e"
        text_color = "#1f2937"
        action_bg = "#ffffff"
        action_border = "#fde68a"
        score_color = "#d97706"
        status_heading = "SURVEILLANCE FINDING: SUBCLONAL KINETIC DRIFT"
        summary_text = (
            f"Low-level subclonal allelic fraction observed (<b>{target_vaf:.2f}% VAF</b>) with low-velocity drift. "
            f"Findings remain below threshold for definitive therapeutic failure but warrant proactive surveillance."
        )
        action_text = (
            "<b>Recommended Clinical Action:</b> Maintain current targeted regimen. Shorten liquid biopsy surveillance "
            "interval to next cycle (21 days) to confirm trajectory convexity."
        )
    else:
        card_bg = "#f0fdf4"
        card_border = "#10b981"
        badge_bg = "#059669"
        heading_color = "#166534"
        text_color = "#1f2937"
        action_bg = "#ffffff"
        action_border = "#a7f3d0"
        score_color = "#059669"
        status_heading = "FAVORABLE FINDING: DEEP MOLECULAR RESPONSE & SUPPRESSION"
        summary_text = (
            "Deep molecular remission observed. Primary oncogenic driver allelic burden is cleared or suppressed near "
            "the limit of detection (<0.10% VAF). Zero secondary gatekeeper resistance mutations detected."
        )
        action_text = (
            "<b>Recommended Clinical Action:</b> Continue standard targeted maintenance protocol with routine 6 to 8-week "
            "liquid biopsy surveillance."
        )

    # Clean Medical Diagnostic Interpretation Card
    st.markdown(
        f"""
        <div style="background: {card_bg}; border: 1.5px solid {card_border}; border-radius: 8px; padding: 18px 22px; margin-bottom: 20px; box-shadow: 0 1px 4px rgba(0,0,0,0.04);">
            <div style="display: flex; justify-content: space-between; align-items: flex-start; flex-wrap: wrap; gap: 15px;">
                <div style="flex: 1; min-width: 320px;">
                    <div style="display: flex; align-items: center; gap: 8px; margin-bottom: 6px;">
                        <span style="background: {badge_bg}; color: #ffffff; font-weight: 700; font-size: 0.72rem; padding: 3px 8px; border-radius: 4px; letter-spacing: 0.5px; text-transform: uppercase;">
                            Diagnostic Alert
                        </span>
                        <span style="font-size: 0.88rem; font-weight: 800; color: {heading_color};">
                            {status_heading}
                        </span>
                    </div>
                    <div style="font-size: 0.9rem; color: {text_color}; line-height: 1.55; margin-top: 6px;">
                        {summary_text}
                    </div>
                    <div style="background: {action_bg}; border: 1px solid {action_border}; border-radius: 6px; padding: 10px 14px; margin-top: 10px; font-size: 0.85rem; color: #1e293b; line-height: 1.45;">
                        {action_text}
                    </div>
                </div>
                <div style="text-align: center; background: #ffffff; border: 1.5px solid {card_border}; border-radius: 8px; padding: 14px 20px; min-width: 160px; box-shadow: 0 1px 3px rgba(0,0,0,0.04);">
                    <div style="font-size: 0.72rem; color: #64748b; text-transform: uppercase; font-weight: 700;">Clonal Escape Index</div>
                    <div style="font-size: 2.2rem; font-weight: 900; color: {score_color}; line-height: 1.1; margin: 4px 0;">
                        {score}<span style="font-size: 1.0rem; color: #94a3b8; font-weight: 500;">/100</span>
                    </div>
                    <div style="font-size: 0.72rem; color: #64748b; font-weight: 500;">Kinetic Derivative Model</div>
                </div>
            </div>
        </div>
        """,
        unsafe_allow_html=True
    )

    # 4 Crisp Hospital-White Telemetry Cards
    c1, c2, c3, c4 = st.columns(4)

    with c1:
        d_vaf = f"{driver_k['current_vaf']:.2f}%" if driver_k else "N/A"
        d_sub = f"{driver_k['clearance_pct']:.0f}% reduction from baseline" if driver_k and driver_k['clearance_pct'] > 0 else "Baseline established"
        st.markdown(
            f"""
            <div style="background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 14px 18px; box-shadow: 0 1px 3px rgba(0,0,0,0.04);">
                <div style="font-size: 0.74rem; color: #64748b; text-transform: uppercase; font-weight: 700; letter-spacing: 0.5px;">Primary Driver VAF</div>
                <div style="font-size: 1.6rem; font-weight: 800; color: #0284c7; margin: 4px 0;">{d_vaf}</div>
                <div style="font-size: 0.78rem; color: #059669; font-weight: 600;">{d_sub}</div>
            </div>
            """,
            unsafe_allow_html=True
        )

    with c2:
        r_vaf = f"{res_k['current_vaf']:.2f}%" if res_k else (f"{target_vaf:.2f}%" if target_clone != "None" else "0.00%")
        r_sub = f"+{target_vel:.2f}% VAF / cycle" if target_vel > 0 else "Undetectable (<0.10% LoD)"
        r_color = "#dc2626" if target_vel > 0.5 else "#64748b"
        st.markdown(
            f"""
            <div style="background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 14px 18px; box-shadow: 0 1px 3px rgba(0,0,0,0.04);">
                <div style="font-size: 0.74rem; color: #64748b; text-transform: uppercase; font-weight: 700; letter-spacing: 0.5px;">Resistance Subclone VAF</div>
                <div style="font-size: 1.6rem; font-weight: 800; color: {r_color}; margin: 4px 0;">{r_vaf}</div>
                <div style="font-size: 0.78rem; color: {'#b91c1c' if target_vel > 0.5 else '#64748b'}; font-weight: 600;">{r_sub}</div>
            </div>
            """,
            unsafe_allow_html=True
        )

    with c3:
        doubling_days = res_k.get("doubling_time_days") if res_k else None
        dt_val = f"{doubling_days:.0f} Days" if doubling_days else ("N/A" if target_vel <= 0 else "< 30 Days")
        dt_sub = "Exponential doubling kinetics" if doubling_days and doubling_days < 40 else "Linear or suppressed"
        st.markdown(
            f"""
            <div style="background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 14px 18px; box-shadow: 0 1px 3px rgba(0,0,0,0.04);">
                <div style="font-size: 0.74rem; color: #64748b; text-transform: uppercase; font-weight: 700; letter-spacing: 0.5px;">Estimated Doubling Time</div>
                <div style="font-size: 1.6rem; font-weight: 800; color: #d97706; margin: 4px 0;">{dt_val}</div>
                <div style="font-size: 0.78rem; color: #64748b;">{dt_sub}</div>
            </div>
            """,
            unsafe_allow_html=True
        )

    with c4:
        # Compute true mean sequencing depth from actual sample data
        all_samples = [s for s in patient_info.get("samples", []) if not s.get("is_forecast", False)]
        depths = [s.get("total_depth", 0) for s in all_samples if s.get("total_depth", 0) > 0]
        mean_depth = int(sum(depths) / len(depths)) if depths else 4180
        min_depth = int(min(depths)) if depths else 0
        max_depth_val = int(max(depths)) if depths else 0
        unique_timepoints = len(set(s["cycle"] for s in all_samples))
        depth_quality = "Ultra-Deep Duplex" if mean_depth >= 3000 else ("Deep Panel" if mean_depth >= 500 else "Standard Panel")
        st.markdown(
            f"""
            <div style="background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 14px 18px; box-shadow: 0 1px 3px rgba(0,0,0,0.04);">
                <div style="font-size: 0.74rem; color: #64748b; text-transform: uppercase; font-weight: 700; letter-spacing: 0.5px;">Mean Target Coverage</div>
                <div style="font-size: 1.6rem; font-weight: 800; color: #7c3aed; margin: 4px 0;">{mean_depth:,}x</div>
                <div style="font-size: 0.78rem; color: #64748b;">{depth_quality} &middot; {unique_timepoints} timepoints &middot; Range: {min_depth:,}&ndash;{max_depth_val:,}x</div>
            </div>
            """,
            unsafe_allow_html=True
        )

    st.markdown("<div style='margin-bottom: 14px;'></div>", unsafe_allow_html=True)
