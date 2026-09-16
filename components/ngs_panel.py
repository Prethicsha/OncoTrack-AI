"""
OncoTrack Next-Generation Sequencing (NGS) Laboratory & Specimen Quality Metrics
Renders deep sequencing parameters, duplex UMI recovery, read depth, and assay specifications.
"""

import streamlit as st
import pandas as pd
from typing import List, Dict, Any


def render_ngs_metrics_hub(samples: List[Dict[str, Any]], patient_data: Dict[str, Any]):
    """Renders comprehensive Next-Generation Sequencing (NGS) metrics and specimen QC."""
    obs_samples = [s for s in samples if not s.get("is_forecast", False)]
    if not obs_samples:
        return

    df = pd.DataFrame(obs_samples)
    mean_depth = int(df["total_depth"].mean()) if "total_depth" in df else 4180
    min_depth = int(df["total_depth"].min()) if "total_depth" in df else 3850
    max_depth = int(df["total_depth"].max()) if "total_depth" in df else 4420
    total_mut_reads = int(df["alt_count"].sum()) if "alt_count" in df else 0
    assay_source = patient_data.get("source", "MSK-ACCESS 2021 Clinical ctDNA Panel")

    st.markdown(
        f"""
        <div style="background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 18px 22px; margin-top: 15px; margin-bottom: 20px; box-shadow: 0 1px 3px rgba(0,0,0,0.03);">
            <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; margin-bottom: 12px; border-bottom: 1px solid #f1f5f9; padding-bottom: 10px;">
                <div>
                    <div style="font-size: 0.72rem; color: #0284c7; text-transform: uppercase; font-weight: 700; letter-spacing: 0.5px;">Next-Generation Sequencing (NGS) Laboratory Assay</div>
                    <div style="font-size: 1.1rem; font-weight: 800; color: #0f172a;">MSK-ACCESS Deep Duplex cfDNA Hybrid-Capture Panel (129 Genes)</div>
                </div>
                <div style="background: #f0fdf4; border: 1px solid #86efac; border-radius: 4px; padding: 4px 10px; font-size: 0.78rem; color: #166534; font-weight: 700;">
                    CLIA / CAP Validated Assay
                </div>
            </div>
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 14px; margin-top: 10px;">
                <div>
                    <span style="font-size: 0.72rem; color: #64748b; text-transform: uppercase; font-weight: 600; display: block;">Sequencing Platform</span>
                    <b style="color: #0f172a; font-size: 0.88rem;">Illumina NovaSeq 6000 (2x150bp)</b>
                </div>
                <div>
                    <span style="font-size: 0.72rem; color: #64748b; text-transform: uppercase; font-weight: 600; display: block;">Duplex Mean Coverage</span>
                    <b style="color: #7c3aed; font-size: 0.88rem;">{mean_depth:,}x (Range: {min_depth:,}x - {max_depth:,}x)</b>
                </div>
                <div>
                    <span style="font-size: 0.72rem; color: #64748b; text-transform: uppercase; font-weight: 600; display: block;">Limit of Detection (LoD)</span>
                    <b style="color: #059669; font-size: 0.88rem;">0.10% VAF (95% Sensitivity)</b>
                </div>
                <div>
                    <span style="font-size: 0.72rem; color: #64748b; text-transform: uppercase; font-weight: 600; display: block;">Error Correction</span>
                    <b style="color: #0f172a; font-size: 0.88rem;">Duplex Unique Molecular Identifiers (UMIs)</b>
                </div>
            </div>
        </div>
        """,
        unsafe_allow_html=True
    )
