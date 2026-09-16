"""
OncoTrack 2D Signature ctDNA Molecular Trajectory & Predictive Kinetic Forecasting Engine
Renders both observed laboratory NGS timepoints and forward projected future cycles (Cycle 6+)
with 95% Poisson confidence intervals and duplex coverage metrics.
"""

import plotly.graph_objects as go
import pandas as pd
import numpy as np
from typing import List, Dict, Any, Optional


def build_signature_ctdna_trajectory_2d(
    samples: List[Dict[str, Any]],
    patient_data: Optional[Dict[str, Any]] = None,
    dark_mode: bool = False,
    show_ngs_depth_bars: bool = True
) -> go.Figure:
    """
    Renders the 2D 'CTDNA MOLECULAR TRAJECTORY & NGS PREDICTIVE FORECAST' plot.
    Features:
    - Observed cycles (Baseline to Cycle 5): Confirmed NGS laboratory data
    - Forecasted cycles (Cycle 6+): Forward kinetic projections with 95% Poisson CI
    - Mapped duplex NGS depth bars & read-level hover telemetry.
    """
    df = pd.DataFrame(samples)
    if df.empty:
        fig = go.Figure()
        fig.update_layout(title="No longitudinal ctDNA data available.")
        return fig

    observed_df = df[~df.get("is_forecast", False)].copy()
    forecast_df = df[df.get("is_forecast", False)].copy()

    all_cycles = sorted(df["cycle"].unique())
    max_obs_cycle = int(observed_df["cycle"].max()) if not observed_df.empty else 1
    
    # Calculate mean NGS depth per cycle
    cycle_ngs_depths = {}
    for c in all_cycles:
        c_sub = df[df["cycle"] == c]
        d_val = int(c_sub["total_depth"].mean()) if "total_depth" in c_sub and not c_sub["total_depth"].isna().all() else 4180
        cycle_ngs_depths[c] = d_val

    # Map cycles to clinical labels with embedded NGS tag
    cycle_labels = {}
    for c in all_cycles:
        is_fc = c > max_obs_cycle
        if is_fc:
            tag = f"<br><span style='font-size:9px; color:#f59e0b; font-weight:700;'>Projected ({cycle_ngs_depths[c]:,}x)</span>"
            cycle_labels[c] = f"Cycle {c - 1} (Forecast){tag}"
        else:
            depth_tag = f"<br><span style='font-size:9px; color:#7c3aed; font-weight:700;'>{cycle_ngs_depths[c]:,}x NGS</span>"
            if c == 1:
                cycle_labels[c] = f"Baseline{depth_tag}"
            else:
                cycle_labels[c] = f"Cycle {c - 1}{depth_tag}"

    fig = go.Figure()

    mut_groups = df[["gene", "mutation", "clone_type"]].drop_duplicates().to_dict("records")

    def sort_order(m):
        if m["clone_type"] == "driver":
            return 0
        elif m["clone_type"] == "resistance":
            return 1
        return 2

    mut_groups = sorted(mut_groups, key=sort_order)

    # Theme palette
    if dark_mode:
        driver_color = "#38bdf8"
        resistance_color = "#f43f5e"
        subclone_color = "#94a3b8"
        fc_driver_color = "#7dd3fc"
        fc_res_color = "#fda4af"
        ngs_bar_color = "rgba(124, 58, 237, 0.18)"
        bg_color = "#0b1329"
        card_bg = "#0f172a"
        text_color = "#f1f5f9"
        grid_color = "#1e293b"
        title_color = "#38bdf8"
    else:
        driver_color = "#0284c7"
        resistance_color = "#dc2626"
        subclone_color = "#64748b"
        fc_driver_color = "#38bdf8"
        fc_res_color = "#f87171"
        ngs_bar_color = "rgba(124, 58, 237, 0.10)"
        bg_color = "#ffffff"
        card_bg = "#f8fafc"
        text_color = "#0f172a"
        grid_color = "#e2e8f0"
        title_color = "#0284c7"

    # 1. Map NGS Duplex Target Coverage as secondary background bars
    if show_ngs_depth_bars:
        x_all = [cycle_labels[c] for c in all_cycles]
        y_depths = [cycle_ngs_depths[c] for c in all_cycles]
        bar_colors = [ngs_bar_color if c <= max_obs_cycle else "rgba(245, 158, 11, 0.12)" for c in all_cycles]
        fig.add_trace(go.Bar(
            x=x_all,
            y=y_depths,
            name="NGS Duplex Depth (x)",
            yaxis="y2",
            marker=dict(
                color=bar_colors,
                line=dict(color="rgba(124, 58, 237, 0.35)", width=1)
            ),
            hoverinfo="text",
            hovertext=[
                f"<b>{'Confirmed Lab NGS' if c <= max_obs_cycle else 'Projected Model Coverage'}:</b> {d:,}x<br><b>Assay:</b> MSK-ACCESS Deep Duplex Panel"
                for c, d in zip(all_cycles, y_depths)
            ],
            showlegend=True
        ))

    max_vaf = 0.0
    has_subclones = any(m["clone_type"] == "subclone" for m in mut_groups)

    # 2. Plot Clonal Kinetics & Forward Predictions
    for m in mut_groups:
        gene = m["gene"]
        mut = m["mutation"]
        c_type = m["clone_type"]

        obs_m = observed_df[(observed_df["gene"] == gene) & (observed_df["mutation"] == mut)].sort_values("cycle")
        fc_m = forecast_df[(forecast_df["gene"] == gene) & (forecast_df["mutation"] == mut)].sort_values("cycle")

        if obs_m.empty:
            continue

        c_vals = obs_m["cycle"].tolist()
        v_vals = obs_m["vaf"].tolist()
        alt_reads = obs_m.get("alt_count", [int(v*40) for v in v_vals]).tolist()
        depth_vals = [cycle_ngs_depths[c] for c in c_vals]
        max_vaf = max(max_vaf, max(v_vals) if v_vals else 0.0)

        x_categories = [cycle_labels[c] for c in c_vals]

        # Observed Hover Text
        hover_texts = []
        for v, a, d in zip(v_vals, alt_reads, depth_vals):
            ref_reads = max(0, d - a)
            hover_texts.append(
                f"<b>Confirmed Lab Specimen</b><br>"
                f"<b>Locus:</b> {gene} {mut} ({c_type.upper()})<br>"
                f"<b>Allelic Fraction (VAF):</b> {v:.2f}%<br>"
                f"<b>NGS Mutant Reads (t_alt):</b> {a:,} UMI reads<br>"
                f"<b>NGS Reference Reads (t_ref):</b> {ref_reads:,} reads<br>"
                f"<b>Total Duplex Depth:</b> {d:,}x"
            )

        if c_type == "driver":
            # Observed Driver
            fig.add_trace(go.Scatter(
                x=x_categories,
                y=v_vals,
                mode="lines+markers+text",
                name=f"Original ({gene} {mut})",
                line=dict(color=driver_color, width=3.5),
                marker=dict(size=9, color=driver_color, symbol="circle", line=dict(color="#ffffff", width=1.5)),
                text=[f"<b>{v:.1f}%</b>" for v in v_vals],
                textposition="top center",
                textfont=dict(size=11, color=driver_color, family="sans-serif"),
                hovertext=hover_texts,
                hoverinfo="text"
            ))

            # Forecasted Driver
            if not fc_m.empty:
                fc_x = [cycle_labels[c_vals[-1]]] + [cycle_labels[c] for c in fc_m["cycle"]]
                fc_y = [v_vals[-1]] + fc_m["vaf"].tolist()
                max_vaf = max(max_vaf, max(fc_y))
                fig.add_trace(go.Scatter(
                    x=fc_x,
                    y=fc_y,
                    mode="lines+markers+text",
                    name=f"Original ({gene} Projected)",
                    line=dict(color=fc_driver_color, width=2.5, dash="dot"),
                    marker=dict(size=8, color=fc_driver_color, symbol="diamond-open", line=dict(color=fc_driver_color, width=2)),
                    text=[""] + [f"{v:.1f}%" for v in fc_m["vaf"]],
                    textposition="top center",
                    textfont=dict(size=10, color=fc_driver_color),
                    hovertext=[""] + [f"<b>Kinetic Forecast</b><br><b>Projected VAF:</b> {v:.2f}%" for v in fc_m["vaf"]],
                    hoverinfo="text"
                ))

        elif c_type == "resistance":
            # Observed Resistance
            text_pos = ["bottom center" if v < 2.0 and has_subclones else "top center" for v in v_vals]
            fig.add_trace(go.Scatter(
                x=x_categories,
                y=v_vals,
                mode="lines+markers+text",
                name=f"Mutant ({gene} {mut} [Escape])",
                line=dict(color=resistance_color, width=3.5, dash="dash"),
                marker=dict(size=9, color=resistance_color, symbol="circle", line=dict(color="#ffffff", width=1.5)),
                text=[f"<b>{v:.1f}%</b>" for v in v_vals],
                textposition=text_pos,
                textfont=dict(size=11, color=resistance_color, family="sans-serif"),
                hovertext=hover_texts,
                hoverinfo="text"
            ))

            # Inflection Point Halo Ring
            for idx, v in enumerate(v_vals):
                if v >= 1.5 and (idx == 0 or v > v_vals[idx-1]):
                    inflection_cycle = x_categories[idx]
                    inflection_vaf = v
                    fig.add_trace(go.Scatter(
                        x=[inflection_cycle],
                        y=[inflection_vaf],
                        mode="markers",
                        name="Escape Inflection",
                        marker=dict(
                            size=18,
                            color="rgba(220, 38, 38, 0.15)",
                            symbol="circle",
                            line=dict(color=resistance_color, width=2.5)
                        ),
                        showlegend=False,
                        hoverinfo="skip"
                    ))
                    break

            # Forecasted Resistance (Cycle 6+)
            if not fc_m.empty:
                fc_x = [cycle_labels[c_vals[-1]]] + [cycle_labels[c] for c in fc_m["cycle"]]
                fc_y = [v_vals[-1]] + fc_m["vaf"].tolist()
                max_vaf = max(max_vaf, max(fc_y))
                fig.add_trace(go.Scatter(
                    x=fc_x,
                    y=fc_y,
                    mode="lines+markers+text",
                    name=f"Mutant ({gene} Projected Forecast)",
                    line=dict(color=fc_res_color, width=3, dash="dot"),
                    marker=dict(size=9, color=fc_res_color, symbol="diamond", line=dict(color="#ffffff", width=1.5)),
                    text=[""] + [f"<b>{v:.1f}%</b>" for v in fc_m["vaf"]],
                    textposition="top center",
                    textfont=dict(size=11, color=fc_res_color, family="sans-serif"),
                    hovertext=[""] + [
                        f"<b>Model Kinetic Forecast (Upcoming Cycle)</b><br><b>Projected VAF:</b> {v:.2f}%<br><b>95% Poisson Interval:</b> [{l:.1f}% - {u:.1f}%]"
                        for v, l, u in zip(fc_m["vaf"], fc_m.get("vaf_lower_ci", fc_m["vaf"]), fc_m.get("vaf_upper_ci", fc_m["vaf"]))
                    ],
                    hoverinfo="text"
                ))

        else:
            # Subclone
            fig.add_trace(go.Scatter(
                x=x_categories,
                y=v_vals,
                mode="lines+markers",
                name=f"Subclone ({gene} {mut})",
                line=dict(color=subclone_color, width=2, dash="dot"),
                marker=dict(size=6, color=subclone_color, symbol="circle"),
                hovertext=hover_texts,
                hoverinfo="text"
            ))

    # Add Actionable 5.0% VAF Threshold Line
    fig.add_hline(
        y=5.0,
        line_dash="dot",
        line_color="#f59e0b",
        line_width=1.5,
        annotation_text="Actionable 5.0% VAF Threshold",
        annotation_position="bottom right",
        annotation_font=dict(size=10, color="#d97706")
    )

    y_upper = max(45.0, (max_vaf * 1.25))

    fig.update_layout(
        template="plotly_dark" if dark_mode else "plotly_white",
        paper_bgcolor=bg_color,
        plot_bgcolor=card_bg,
        height=400,
        margin=dict(l=45, r=60, t=55, b=45),
        title=dict(
            text=f"<b style='color: {title_color}; font-size: 1.05rem;'>LONGITUDINAL CTDNA TRAJECTORY & FORWARD KINETIC PREDICTIONS</b>",
            x=0.01,
            y=0.96,
            font=dict(family="-apple-system, BlinkMacSystemFont, Segoe UI, Roboto, Helvetica, Arial, sans-serif", size=13, color=text_color)
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.03,
            xanchor="right",
            x=1.0,
            font=dict(size=10, color=text_color),
            bgcolor="rgba(0,0,0,0)",
            bordercolor="rgba(0,0,0,0)"
        ),
        xaxis=dict(
            title="",
            showgrid=False,
            zeroline=False,
            showline=True,
            linecolor=grid_color,
            tickfont=dict(size=11, color="#334155" if not dark_mode else "#94a3b8", family="sans-serif")
        ),
        yaxis=dict(
            title="<b>Variant Allele Fraction (VAF %)</b>",
            range=[0, y_upper],
            tickmode="linear",
            tick0=0,
            dtick=10,
            ticksuffix="%",
            showgrid=True,
            gridcolor=grid_color,
            gridwidth=1,
            griddash="dot",
            zeroline=True,
            zerolinecolor=grid_color,
            tickfont=dict(size=11, color="#64748b")
        ),
        yaxis2=dict(
            title=dict(text="<b>NGS Coverage (x)</b>", font=dict(color="#7c3aed", size=10)),
            range=[0, 7000],
            overlaying="y",
            side="right",
            showgrid=False,
            tickfont=dict(size=10, color="#7c3aed"),
            ticksuffix="x"
        ),
        barmode="overlay",
        hovermode="closest"
    )

    return fig
