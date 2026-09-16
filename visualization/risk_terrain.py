"""
OncoTrack Analytical Risk Terrain & Clonal Architecture Engine
Renders the 3D Kinetic Analytical Space and 2D Muller Clonal Architecture in clean clinical light theme.
"""

import plotly.graph_objects as go
import numpy as np
import pandas as pd
from typing import Dict, List, Any


def build_3d_risk_terrain(patient_score_data: Dict[str, Any], patient_id: str) -> go.Figure:
    """
    Renders the 3D Analytical Space in clean clinical light theme:
    X = Velocity (% VAF / cycle), Y = Persistence (Cycles detected), Z = Clonal Escape Score
    """
    vel_grid = np.linspace(-5.0, 10.0, 25)
    pers_grid = np.linspace(1.0, 6.0, 25)
    V, P = np.meshgrid(vel_grid, pers_grid)

    Z = np.clip((V * 4.5) + (P * 8.0) + 15.0, 0, 100)

    fig = go.Figure()

    fig.add_trace(go.Surface(
        x=V, y=P, z=Z,
        colorscale=[
            [0.0, "rgba(16, 185, 129, 0.5)"],   # Green (Response / Low Risk)
            [0.45, "rgba(245, 158, 11, 0.5)"],  # Amber (Surveillance)
            [1.0, "rgba(220, 38, 38, 0.65)"]    # Red (Critical Escape Alert)
        ],
        showscale=False,
        name="Analytical Risk Manifold",
        hoverinfo="skip"
    ))

    pt_vel = patient_score_data.get("target_velocity", 0.0)
    pt_pers = 5.0 if patient_score_data.get("score", 0) > 60 else 3.0
    pt_score = patient_score_data.get("score", 0)
    pt_status = patient_score_data.get("status", "MONITOR")

    marker_color = "#dc2626" if pt_score >= 70 else ("#d97706" if pt_score >= 40 else "#059669")

    fig.add_trace(go.Scatter3d(
        x=[pt_vel],
        y=[pt_pers],
        z=[pt_score],
        mode="markers+text",
        name=f"Patient {patient_id}",
        marker=dict(
            size=12,
            color=marker_color,
            symbol="diamond",
            line=dict(color="#ffffff", width=2)
        ),
        text=[f"<b>{patient_id}</b><br>Score: {pt_score}/100<br>{pt_status}"],
        textposition="top center",
        hovertext=[
            f"<b>Patient:</b> {patient_id}<br>"
            f"<b>Kinetic Velocity:</b> +{pt_vel:.2f}% / cycle<br>"
            f"<b>Persistence:</b> {pt_pers:.0f} cycles<br>"
            f"<b>Clonal Escape Index:</b> {pt_score}/100<br>"
            f"<b>Analytical Regime:</b> {pt_status}"
        ],
        hoverinfo="text"
    ))

    fig.update_layout(
        template="plotly_white",
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        margin=dict(l=0, r=0, t=30, b=0),
        height=520,
        scene=dict(
            xaxis=dict(
                title="<b>Resistance Velocity (ΔVAF/cycle)</b>",
                gridcolor="#e2e8f0",
                backgroundcolor="#f8fafc",
                showbackground=True
            ),
            yaxis=dict(
                title="<b>Persistence (Cycles)</b>",
                gridcolor="#e2e8f0",
                backgroundcolor="#f8fafc",
                showbackground=True
            ),
            zaxis=dict(
                title="<b>Clonal Escape Index (0-100)</b>",
                range=[0, 100],
                gridcolor="#e2e8f0",
                backgroundcolor="#f8fafc",
                showbackground=True
            ),
            camera=dict(
                eye=dict(x=1.5, y=-1.5, z=1.1)
            )
        )
    )

    return fig


def build_clonal_muller_plot(samples: List[Dict[str, Any]]) -> go.Figure:
    """
    Renders standard oncology clonal evolution area chart in clean clinical light theme.
    """
    df = pd.DataFrame(samples)
    if df.empty:
        return go.Figure()

    df = df[~df.get("is_forecast", False)]
    mutations = df[["gene", "mutation", "clone_type"]].drop_duplicates().to_dict("records")
    cycles = sorted(df["cycle"].unique())

    fig = go.Figure()

    color_dict = {
        "driver": "#0284c7",
        "resistance": "#dc2626",
        "subclone": "#64748b"
    }

    for mut in mutations:
        mut_key = f"{mut['gene']} {mut['mutation']}"
        sub_df = df[(df["gene"] == mut["gene"]) & (df["mutation"] == mut["mutation"])].sort_values("cycle")

        vafs_map = dict(zip(sub_df["cycle"], sub_df["vaf"]))
        y_vals = [vafs_map.get(c, 0.0) for c in cycles]

        fig.add_trace(go.Scatter(
            x=[f"Cycle {c}" for c in cycles],
            y=y_vals,
            mode="lines+markers",
            name=f"{mut_key} ({mut['clone_type'].capitalize()})",
            line=dict(color=color_dict.get(mut["clone_type"], "#64748b"), width=3),
            marker=dict(size=8),
            stackgroup="one"
        ))

    fig.update_layout(
        template="plotly_white",
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        margin=dict(l=40, r=20, t=30, b=40),
        height=320,
        title=dict(
            text="<b>Longitudinal Clonal Composition & Allelic Sweep (Stacked VAF %)</b>",
            font=dict(size=13, color="#0f172a")
        ),
        xaxis=dict(gridcolor="#e2e8f0", title="<b>Treatment Evaluation Cycles</b>"),
        yaxis=dict(gridcolor="#e2e8f0", title="<b>Cumulative VAF (%)</b>"),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1.0,
            font=dict(size=10, color="#334155")
        )
    )

    return fig
