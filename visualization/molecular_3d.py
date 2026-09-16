"""
OncoTrack 3D Molecular Evolution Engine
Renders high-resolution 3D WebGL trajectories with active camera presets,
uirevision handling, and high-contrast clinical light styling.
"""

import plotly.graph_objects as go
import numpy as np
import pandas as pd
from typing import Dict, List, Any, Optional
from engine.trajectory import interpolate_3d_spline

# Clinical Color Semantics
COLOR_MAP = {
    "driver": {
        "observed": "#0284c7",     # Deep Cobalt Blue
        "forecast": "#38bdf8",     # Sky Blue
        "marker": "#0369a1"
    },
    "resistance": {
        "observed": "#dc2626",     # Vivid Medical Crimson
        "forecast": "#f87171",     # Soft Red
        "marker": "#b91c1c"
    },
    "subclone": {
        "observed": "#475569",     # Slate Gray
        "forecast": "#94a3b8",     # Light Slate
        "marker": "#334155"
    }
}

# Mathematical Camera Vectors (with explicit up and center vectors to prevent gimbal lock)
CAMERA_PRESETS = {
    "Overview (Perspective View)": {
        "eye": dict(x=1.65, y=1.65, z=1.15),
        "center": dict(x=0.0, y=0.0, z=-0.1),
        "up": dict(x=0.0, y=0.0, z=1.0)
    },
    "Temporal Flow (Top View)": {
        "eye": dict(x=0.001, y=0.001, z=2.8),
        "center": dict(x=0.0, y=0.0, z=0.0),
        "up": dict(x=0.0, y=1.0, z=0.0)
    },
    "VAF Kinetic Velocity (Side View)": {
        "eye": dict(x=0.001, y=-2.5, z=0.35),
        "center": dict(x=0.0, y=0.0, z=-0.05),
        "up": dict(x=0.0, y=0.0, z=1.0)
    },
    "Resistance Focus (Target Clone)": {
        "eye": dict(x=1.1, y=0.75, z=0.55),
        "center": dict(x=0.2, y=0.15, z=0.1),
        "up": dict(x=0.0, y=0.0, z=1.0)
    },
    "Current Assessment Timepoint": {
        "eye": dict(x=0.35, y=-1.9, z=0.8),
        "center": dict(x=0.1, y=0.0, z=0.0),
        "up": dict(x=0.0, y=0.0, z=1.0)
    }
}


def build_3d_molecular_landscape(
    samples: List[Dict[str, Any]],
    patient_id: str = "MSK-LC-003",
    camera_view: str = "Overview (Perspective View)",
    show_threshold_plane: bool = True,
    show_current_cycle_plane: bool = True,
    show_uncertainty_envelope: bool = True,
    max_cycle_filter: Optional[int] = None,
    filter_mutation: str = "All Clones"
) -> go.Figure:
    """
    Constructs the 3D longitudinal molecular evolution space.
    Axes: X = Treatment Cycle, Y = Mutation / Clonal Lane, Z = VAF (%)
    """
    df = pd.DataFrame(samples)
    if df.empty:
        fig = go.Figure()
        fig.update_layout(title="No molecular samples available.")
        return fig

    if max_cycle_filter is not None:
        df = df[df["cycle"] <= max_cycle_filter]

    mut_records = df[["gene", "mutation", "clone_type"]].drop_duplicates().to_dict("records")
    
    def sort_key(m):
        if m["clone_type"] == "driver":
            return (0, m["gene"])
        elif m["clone_type"] == "resistance":
            return (1, m["gene"])
        return (2, m["gene"])

    mut_records = sorted(mut_records, key=sort_key)
    mutation_y_map = {f"{m['gene']} {m['mutation']}": idx for idx, m in enumerate(mut_records)}
    y_tick_labels = [f"<b>{m['gene']}</b><br>{m['mutation']}" for m in mut_records]

    fig = go.Figure()

    observed_df = df[~df.get("is_forecast", False)]
    forecast_df = df[df.get("is_forecast", False)]

    max_observed_cycle = int(observed_df["cycle"].max()) if not observed_df.empty else 1
    max_total_cycle = int(df["cycle"].max()) if not df.empty else 1
    max_vaf = max(45.0, float(df["vaf"].max()) * 1.15)

    # 1. Plot Trajectories for each mutation
    for mut_info in mut_records:
        mut_key = f"{mut_info['gene']} {mut_info['mutation']}"
        clone_type = mut_info["clone_type"]
        y_pos = mutation_y_map[mut_key]

        if filter_mutation != "All Clones" and filter_mutation != mut_key:
            continue

        c_palette = COLOR_MAP.get(clone_type, COLOR_MAP["subclone"])

        # Observed data points
        obs_m = observed_df[(observed_df["gene"] == mut_info["gene"]) & (observed_df["mutation"] == mut_info["mutation"])].sort_values("cycle")
        if not obs_m.empty:
            obs_cycles = obs_m["cycle"].tolist()
            obs_vafs = obs_m["vaf"].tolist()
            obs_dates = obs_m.get("date", obs_m.get("cycle_label", "")).tolist()
            obs_depths = obs_m.get("total_depth", [4000]*len(obs_cycles)).tolist()

            # Smooth curve
            if len(obs_cycles) >= 2:
                sx, sz = interpolate_3d_spline(obs_cycles, obs_vafs, num_points=40)
                sy = [y_pos] * len(sx)

                fig.add_trace(go.Scatter3d(
                    x=sx, y=sy, z=sz,
                    mode="lines",
                    line=dict(color=c_palette["observed"], width=5.5),
                    showlegend=False,
                    hoverinfo="skip"
                ))

            # Marker points with clinical hover text
            hover_texts = []
            for c, v, dt, dp in zip(obs_cycles, obs_vafs, obs_dates, obs_depths):
                hover_texts.append(
                    f"<b>Patient:</b> {patient_id}<br>"
                    f"<b>Timepoint:</b> Cycle {c} ({dt})<br>"
                    f"<b>Locus:</b> {mut_info['gene']} {mut_info['mutation']}<br>"
                    f"<b>Lineage:</b> {clone_type.upper()}<br>"
                    f"<b>Allelic Fraction (VAF):</b> {v:.2f}%<br>"
                    f"<b>Duplex Coverage:</b> {dp:,}x<br>"
                    f"<b>Observation Type:</b> Confirmed Laboratory Assay"
                )

            lineage_label = "Primary Sensitizing Driver" if clone_type == "driver" else ("Acquired Resistance Locus" if clone_type == "resistance" else "Passenger Subclone")
            fig.add_trace(go.Scatter3d(
                x=obs_cycles,
                y=[y_pos] * len(obs_cycles),
                z=obs_vafs,
                mode="markers",
                name=f"{mut_key} ({lineage_label})",
                marker=dict(
                    size=7.5,
                    color=c_palette["marker"],
                    symbol="circle",
                    line=dict(color="#ffffff", width=1.5)
                ),
                text=hover_texts,
                hoverinfo="text"
            ))

        # Forecast data points
        fc_m = forecast_df[(forecast_df["gene"] == mut_info["gene"]) & (forecast_df["mutation"] == mut_info["mutation"])].sort_values("cycle")
        if not fc_m.empty and not obs_m.empty:
            fc_cycles = [obs_cycles[-1]] + fc_m["cycle"].tolist()
            fc_vafs = [obs_vafs[-1]] + fc_m["vaf"].tolist()
            fc_lowers = [obs_vafs[-1]] + fc_m.get("vaf_lower_ci", fc_m["vaf"]).tolist()
            fc_uppers = [obs_vafs[-1]] + fc_m.get("vaf_upper_ci", fc_m["vaf"]).tolist()

            sx_fc, sz_fc = interpolate_3d_spline(fc_cycles, fc_vafs, num_points=25)
            sy_fc = [y_pos] * len(sx_fc)

            fig.add_trace(go.Scatter3d(
                x=sx_fc, y=sy_fc, z=sz_fc,
                mode="lines",
                line=dict(color=c_palette["forecast"], width=3.5, dash="dash"),
                name=f"{mut_key} (Kinetic Projection)",
                hoverinfo="skip"
            ))

            fc_hover = [
                f"<b>Model Kinetic Forecast</b><br>"
                f"<b>Locus:</b> {mut_info['gene']} {mut_info['mutation']}<br>"
                f"<b>Projected Cycle:</b> {c}<br>"
                f"<b>Expected VAF:</b> {v:.2f}%<br>"
                f"<b>95% Poisson Interval:</b> [{l:.1f}% - {u:.1f}%]"
                for c, v, l, u in zip(fc_m["cycle"], fc_m["vaf"], fc_m["vaf_lower_ci"], fc_m["vaf_upper_ci"])
            ]

            fig.add_trace(go.Scatter3d(
                x=fc_m["cycle"].tolist(),
                y=[y_pos] * len(fc_m),
                z=fc_m["vaf"].tolist(),
                mode="markers",
                marker=dict(
                    size=6.5,
                    color=c_palette["forecast"],
                    symbol="diamond-open",
                    line=dict(color=c_palette["forecast"], width=2)
                ),
                text=fc_hover,
                hoverinfo="text",
                showlegend=False
            ))

            if show_uncertainty_envelope and clone_type == "resistance":
                ribbon_x = fc_cycles + fc_cycles[::-1]
                ribbon_y = [y_pos] * len(ribbon_x)
                ribbon_z = fc_lowers + fc_uppers[::-1]
                fig.add_trace(go.Scatter3d(
                    x=ribbon_x, y=ribbon_y, z=ribbon_z,
                    mode="lines",
                    surfaceaxis=1,
                    surfacecolor="rgba(220, 38, 38, 0.14)",
                    line=dict(color="rgba(220, 38, 38, 0.3)", width=1),
                    name="Poisson Uncertainty Envelope (95% CI)",
                    hoverinfo="skip"
                ))

    # 2. Orthogonal Plane: Current Assessment Boundary
    if show_current_cycle_plane and max_observed_cycle > 1:
        plane_y = np.linspace(-0.3, len(mut_records) - 0.7, 4)
        plane_z = np.linspace(0, max_vaf, 4)
        PY, PZ = np.meshgrid(plane_y, plane_z)
        PX = np.full_like(PY, max_observed_cycle)

        fig.add_trace(go.Surface(
            x=PX, y=PY, z=PZ,
            showscale=False,
            colorscale=[[0, "rgba(2, 132, 199, 0.16)"], [1, "rgba(2, 132, 199, 0.16)"]],
            name="Current Evaluation Boundary",
            hoverinfo="skip"
        ))

    # 3. Orthogonal Plane: Demonstration Threshold (5.0% VAF)
    if show_threshold_plane:
        thresh_x = np.linspace(0.8, max_total_cycle + 0.2, 4)
        thresh_y = np.linspace(-0.3, len(mut_records) - 0.7, 4)
        TX, TY = np.meshgrid(thresh_x, thresh_y)
        TZ = np.full_like(TX, 5.0)

        fig.add_trace(go.Surface(
            x=TX, y=TY, z=TZ,
            showscale=False,
            colorscale=[[0, "rgba(245, 158, 11, 0.14)"], [1, "rgba(245, 158, 11, 0.14)"]],
            name="Actionable Demonstration Threshold (5.0% VAF)",
            hoverinfo="skip"
        ))

    cam_cfg = CAMERA_PRESETS.get(camera_view, CAMERA_PRESETS["Overview (Perspective View)"])

    # High-Resolution Clinical Light Workstation Theme with uirevision linked to camera_view
    fig.update_layout(
        uirevision=camera_view,  # Forces Plotly to reapply camera when camera_view changes
        template="plotly_white",
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        margin=dict(l=0, r=0, t=30, b=0),
        height=620,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="center",
            x=0.5,
            font=dict(size=11, color="#334155"),
            bgcolor="rgba(248, 250, 252, 0.95)",
            bordercolor="#cbd5e1",
            borderwidth=1
        ),
        scene=dict(
            camera=dict(
                eye=cam_cfg["eye"],
                center=cam_cfg.get("center", dict(x=0, y=0, z=0)),
                up=cam_cfg.get("up", dict(x=0, y=0, z=1))
            ),
            xaxis=dict(
                title=dict(text="<b>Treatment Timeline (Cycles)</b>", font=dict(size=11, color="#475569")),
                tickmode="linear",
                tick0=1,
                dtick=1,
                range=[0.8, max_total_cycle + 0.3],
                gridcolor="#e2e8f0",
                zerolinecolor="#cbd5e1",
                backgroundcolor="#f8fafc",
                showbackground=True
            ),
            yaxis=dict(
                title=dict(text="<b>Clonal Subpopulation / Gene Locus</b>", font=dict(size=11, color="#475569")),
                tickvals=list(range(len(mut_records))),
                ticktext=y_tick_labels,
                range=[-0.4, len(mut_records) - 0.6],
                gridcolor="#e2e8f0",
                zerolinecolor="#cbd5e1",
                backgroundcolor="#f8fafc",
                showbackground=True
            ),
            zaxis=dict(
                title=dict(text="<b>Variant Allele Fraction (VAF %)</b>", font=dict(size=11, color="#475569")),
                range=[0, max_vaf],
                gridcolor="#e2e8f0",
                zerolinecolor="#cbd5e1",
                backgroundcolor="#f8fafc",
                showbackground=True
            ),
            aspectmode="manual",
            aspectratio=dict(x=1.8, y=1.2, z=1.0)
        )
    )

    return fig
