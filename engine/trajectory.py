"""
OncoTrack Trajectory Kinetics Engine
Performs numerical differentiation (velocity & acceleration calculus),
clonal doubling time estimation, and smooth spline interpolation for 3D trajectories.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Any
from scipy.interpolate import PchipInterpolator


def compute_clonal_kinetics(samples: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """
    Computes velocity (dVAF/dt), acceleration (d2VAF/dt2), and doubling times
    for each distinct mutation / clonal lineage.
    """
    df = pd.DataFrame(samples)
    if df.empty:
        return {}

    kinetics = {}
    mutations = df[["gene", "mutation", "clone_type"]].drop_duplicates().to_dict("records")

    for mut_info in mutations:
        mut_key = f"{mut_info['gene']} {mut_info['mutation']}"
        sub_df = df[(df["gene"] == mut_info["gene"]) & (df["mutation"] == mut_info["mutation"])].sort_values("cycle")

        cycles = sub_df["cycle"].values
        days = sub_df["days"].values if "days" in sub_df.columns else (cycles - 1) * 28
        vafs = sub_df["vaf"].values

        if len(cycles) < 2:
            continue

        # Numerical first derivative (velocity in % VAF per cycle and % per day)
        delta_vaf = np.diff(vafs)
        delta_cycles = np.diff(cycles)
        velocity_per_cycle = delta_vaf / np.maximum(delta_cycles, 1)

        delta_days = np.diff(days)
        velocity_per_day = delta_vaf / np.maximum(delta_days, 1)

        # Numerical second derivative (acceleration / convexity)
        if len(velocity_per_cycle) >= 2:
            acceleration = np.diff(velocity_per_cycle) / np.maximum(delta_cycles[1:], 1)
            recent_accel = float(acceleration[-1])
        else:
            recent_accel = 0.0

        recent_velocity = float(velocity_per_cycle[-1])
        recent_vaf = float(vafs[-1])
        baseline_vaf = float(vafs[0])
        total_delta_vaf = recent_vaf - baseline_vaf

        # Clonal doubling time calculation: T_double = dt * ln(2) / ln(VAF_t / VAF_t-1)
        doubling_time_days = None
        if len(vafs) >= 2 and vafs[-1] > vafs[-2] and vafs[-2] > 0.05:
            dt_days = float(days[-1] - days[-2])
            ratio = vafs[-1] / vafs[-2]
            if ratio > 1.0:
                doubling_time_days = round(dt_days * np.log(2.0) / np.log(ratio), 1)

        # Clearance rate for drivers
        clearance_pct = 0.0
        if baseline_vaf > 0:
            clearance_pct = max(0.0, min(100.0, ((baseline_vaf - recent_vaf) / baseline_vaf) * 100.0))

        kinetics[mut_key] = {
            "gene": mut_info["gene"],
            "mutation": mut_info["mutation"],
            "clone_type": mut_info["clone_type"],
            "cycles": cycles.tolist(),
            "vafs": vafs.tolist(),
            "baseline_vaf": baseline_vaf,
            "current_vaf": recent_vaf,
            "total_delta_vaf": round(total_delta_vaf, 2),
            "recent_velocity_per_cycle": round(recent_velocity, 3),
            "recent_velocity_per_day": round(float(velocity_per_day[-1]), 4),
            "recent_acceleration": round(recent_accel, 4),
            "doubling_time_days": doubling_time_days,
            "clearance_pct": round(clearance_pct, 1),
            "is_increasing": recent_velocity > 0.2,
            "is_accelerating": recent_accel > 0.05
        }

    return kinetics


def interpolate_3d_spline(cycles: List[float], vafs: List[float], num_points: int = 50) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generates smooth, biologically monotonic PCHIP interpolated curves for 3D rendering.
    Prevents artificial oscillations/overshoots seen with naive cubic splines.
    """
    if len(cycles) < 2:
        return np.array(cycles), np.array(vafs)

    cycles_arr = np.array(cycles, dtype=float)
    vafs_arr = np.array(vafs, dtype=float)

    # Sort strictly
    idx = np.argsort(cycles_arr)
    cycles_arr = cycles_arr[idx]
    vafs_arr = vafs_arr[idx]

    # Remove duplicates if any
    unique_cycles, u_indices = np.unique(cycles_arr, return_index=True)
    unique_vafs = vafs_arr[u_indices]

    if len(unique_cycles) < 2:
        return unique_cycles, unique_vafs

    pchip = PchipInterpolator(unique_cycles, unique_vafs)
    smooth_x = np.linspace(unique_cycles.min(), unique_cycles.max(), num_points)
    smooth_z = np.clip(pchip(smooth_x), a_min=0.0, a_max=100.0)

    return smooth_x, smooth_z
