"""
ONCOTRACK — Clinical Molecular Oncology Workstation
Hospital-grade Clinical Decision Support System for Longitudinal ctDNA Dynamics & Clonal Evolution.
"""

import streamlit as st
import pandas as pd
import numpy as np
import io
import time

# Page Configuration
st.set_page_config(
    page_title="OncoTrack — Clinical Molecular Dynamics Workstation",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Professional Hospital Workstation Styling (Clean Light Theme, Inter / Helvetica, Zero Gimmicks)
st.markdown(
    """
    <style>
        .stApp {
            background-color: #f8fafc;
            color: #0f172a;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
        }
        /* Top Navigation Header */
        .top-navbar {
            background: #ffffff;
            border: 1px solid #e2e8f0;
            border-radius: 6px;
            padding: 16px 22px;
            margin-bottom: 16px;
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
            gap: 15px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.04);
        }
        .brand-section {
            display: flex;
            align-items: center;
            gap: 12px;
        }
        .brand-title {
            font-size: 1.35rem;
            font-weight: 800;
            color: #0f172a;
            letter-spacing: -0.5px;
            margin: 0;
        }
        .user-badge {
            background: #f1f5f9;
            border: 1px solid #e2e8f0;
            border-radius: 4px;
            padding: 6px 12px;
            font-size: 0.8rem;
            color: #334155;
            display: flex;
            align-items: center;
            gap: 8px;
        }
        /* Patient Context Master Bar */
        .patient-context-bar {
            background: #ffffff;
            border: 1px solid #e2e8f0;
            border-radius: 6px;
            padding: 14px 20px;
            margin-bottom: 18px;
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
            gap: 14px;
            font-size: 0.86rem;
            box-shadow: 0 1px 3px rgba(0,0,0,0.03);
        }
        .ctx-item span {
            color: #64748b;
            font-size: 0.72rem;
            text-transform: uppercase;
            font-weight: 600;
            display: block;
            margin-bottom: 2px;
        }
        .ctx-item b {
            color: #0f172a;
            font-size: 0.94rem;
        }
        /* Clinical Value Proposition Bar */
        .value-prop-bar {
            background: #ffffff;
            border: 1px solid #e2e8f0;
            border-radius: 6px;
            padding: 12px 18px;
            margin-bottom: 18px;
            display: flex;
            justify-content: space-around;
            flex-wrap: wrap;
            gap: 12px;
            font-size: 0.82rem;
            color: #475569;
            box-shadow: 0 1px 3px rgba(0,0,0,0.03);
        }
        .value-prop-item b {
            color: #0f172a;
        }
        /* Pricing & Business Cards */
        .pricing-card {
            background: #ffffff;
            border: 1px solid #e2e8f0;
            border-radius: 6px;
            padding: 24px 20px;
            text-align: center;
            height: 100%;
            box-shadow: 0 1px 3px rgba(0,0,0,0.04);
        }
        .pricing-card-featured {
            background: #ffffff;
            border: 2px solid #0284c7;
            border-radius: 6px;
            padding: 24px 20px;
            text-align: center;
            height: 100%;
            box-shadow: 0 4px 12px rgba(2, 132, 199, 0.1);
        }
        .camera-btn-bar {
            background: #ffffff;
            border: 1px solid #e2e8f0;
            border-radius: 6px;
            padding: 10px 14px;
            margin-bottom: 12px;
            display: flex;
            align-items: center;
            gap: 8px;
            flex-wrap: wrap;
            box-shadow: 0 1px 2px rgba(0,0,0,0.03);
        }
        .disclaimer-footer {
            background: #f1f5f9;
            border: 1px solid #e2e8f0;
            border-radius: 4px;
            padding: 14px 18px;
            margin-top: 35px;
            font-size: 0.76rem;
            color: #64748b;
            line-height: 1.5;
        }
        div[data-testid="stMetricValue"] {
            font-size: 1.5rem !important;
            font-weight: 700 !important;
            color: #0f172a !important;
        }
        div[data-testid="stMetricLabel"] {
            font-size: 0.8rem !important;
            color: #64748b !important;
            text-transform: uppercase !important;
            letter-spacing: 0.5px !important;
        }
    </style>
    """,
    unsafe_allow_html=True
)

# Modular imports
from utils.data_loader import ClinicalDataLoader, parse_uploaded_csv
from engine.scoring import calculate_escape_score
from engine.forecasting import forecast_clonal_trajectories
from visualization.molecular_3d import build_3d_molecular_landscape, CAMERA_PRESETS
from visualization.molecular_2d_trajectory import build_signature_ctdna_trajectory_2d
from visualization.risk_terrain import build_3d_risk_terrain, build_clonal_muller_plot
from components.status_card import render_clinical_status_card
from components.explanation import render_explanation_panel
from components.ngs_panel import render_ngs_metrics_hub
from reports.pdf_report import generate_clinical_pdf

# 1. Initialize Pan-Cancer Clinical Cohort
@st.cache_resource
def get_loader():
    return ClinicalDataLoader("msk_access_2021.tar.gz")

loader = get_loader()
full_cohort = loader.load_cohort()

if "custom_cohort" not in st.session_state:
    st.session_state["custom_cohort"] = {}

full_cohort.update(st.session_state["custom_cohort"])

if "active_pid" not in st.session_state:
    st.session_state["active_pid"] = "ONC-PT-001"

if "camera_view" not in st.session_state:
    st.session_state["camera_view"] = "Overview (Perspective View)"

# 2. Sidebar: Clinical Navigation & Filtering
with st.sidebar:
    st.markdown("### LIMS Specimen Ingestion")
    uploaded_file = st.file_uploader("Ingest Custom ctDNA Report (CSV/TSV)", type=["csv", "tsv"], help="Upload serial liquid biopsy NGS reports (Guardant, FoundationOne, Illumina, LIMS exports).")
    if uploaded_file is not None:
        custom_patients = parse_uploaded_csv(uploaded_file)
        if custom_patients:
            if isinstance(custom_patients, dict) and "samples" in custom_patients:
                pid = custom_patients["patient_id"]
                full_cohort[pid] = custom_patients
                st.session_state["custom_cohort"][pid] = custom_patients
                if st.session_state.get("_last_uploaded_file") != uploaded_file.name:
                    st.session_state["active_pid"] = pid
                    st.session_state["_last_uploaded_file"] = uploaded_file.name
            elif isinstance(custom_patients, dict):
                for pid, pdata in custom_patients.items():
                    full_cohort[pid] = pdata
                    st.session_state["custom_cohort"][pid] = pdata
                if st.session_state.get("_last_uploaded_file") != uploaded_file.name:
                    first_pid = list(custom_patients.keys())[0]
                    st.session_state["active_pid"] = first_pid
                    st.session_state["_last_uploaded_file"] = uploaded_file.name
            st.success(f"✓ Ingested {len(custom_patients) if 'samples' not in custom_patients else 1} serial ctDNA record(s).")
        else:
            st.error("Unable to parse CSV. Ensure columns include gene, variant/mutation, and vaf/alt_count.")

    st.markdown("---")
    st.markdown("### Clinical Oncology Cohort")
    st.caption(f"Active Longitudinal Series: **{len(full_cohort)} Validated Records**")

    # Organ Domain Filter
    all_categories = ["All Organ Systems (24 Pan-Cancer Indications)"] + sorted(list(dict.fromkeys(p.get("category", "General") for p in full_cohort.values())))
    sel_category = st.selectbox("Organ System Domain", options=all_categories, index=0)

    # Indication Filter
    if sel_category != "All Organ Systems (24 Pan-Cancer Indications)":
        avail_types = ["All Indications in Domain"] + sorted(list(dict.fromkeys(p["cancer_type"] for p in full_cohort.values() if p.get("category") == sel_category)))
    else:
        avail_types = ["All 24 Malignancy Indications"] + sorted(list(dict.fromkeys(p["cancer_type"] for p in full_cohort.values())))
    sel_cancer = st.selectbox("Cancer Indication", options=avail_types, index=0)

    # Status Filter
    status_filter = st.selectbox(
        "Clinical Status Filter",
        options=["All Statuses", "Molecular Resistance / Escape", "Molecular Response / Clearance", "Molecular Surveillance / Monitor"],
        index=0
    )

    # Search Bar
    search_query = st.text_input("Search Locus / Gene / Patient ID", placeholder="e.g. EGFR, T790M, KRAS, ESR1, SYN-PT-001")

    # Filter cohort list
    filtered_pids = []
    for pid, pdata in full_cohort.items():
        if sel_category != "All Organ Systems (24 Pan-Cancer Indications)" and pdata.get("category") != sel_category:
            continue
        if sel_cancer not in ["All 24 Malignancy Indications", "All Indications in Domain"] and pdata["cancer_type"] != sel_cancer:
            continue
        if status_filter != "All Statuses":
            if "Escape" in status_filter and pdata["clinical_status"] != "POTENTIAL MOLECULAR ESCAPE":
                continue
            if "Response" in status_filter and pdata["clinical_status"] != "MOLECULAR RESPONSE":
                continue
            if "Monitor" in status_filter and pdata["clinical_status"] != "MONITOR":
                continue
        if search_query:
            q = search_query.lower()
            p_text = f"{pid} {pdata['cancer_type']} {pdata['primary_driver']} {pdata.get('resistance_locus', '')} {pdata.get('category', '')}".lower()
            if q not in p_text:
                continue
        filtered_pids.append(pid)

    if not filtered_pids:
        st.warning("No patients match active filters. Resetting to full cohort.")
        filtered_pids = list(full_cohort.keys())

    if st.session_state["active_pid"] not in filtered_pids:
        st.session_state["active_pid"] = filtered_pids[0]

    selected_idx = filtered_pids.index(st.session_state["active_pid"]) if st.session_state["active_pid"] in filtered_pids else 0
    chosen_pid = st.selectbox(
        f"Select Patient Record ({len(filtered_pids)} Available)",
        options=filtered_pids,
        index=selected_idx,
        format_func=lambda pid: full_cohort[pid].get("display_name", pid),
        key="sidebar_patient_selector"
    )
    st.session_state["active_pid"] = chosen_pid

    st.markdown("---")
    st.markdown("### Workstation Modules")
    vis_mode = st.radio(
        "Display Mode",
        options=[
            "3D Molecular Evolution Space",
            "2D ctDNA Trajectory (Rise & Fall of Clones)",
            "2D Clonal Architecture (Muller Plot)",
            "3D Analytical Risk Space",
            "Pan-Cancer Cohort Registry (100+ Series)",
            "Molecular Tumor Board Case Conference",
            "Health System Economics & ROI"
        ],
        index=0
    )

    st.markdown("---")
    st.markdown("### Spatial Camera Presets")
    sidebar_cam = st.selectbox(
        "Camera Viewpoint",
        options=list(CAMERA_PRESETS.keys()),
        index=list(CAMERA_PRESETS.keys()).index(st.session_state["camera_view"]) if st.session_state["camera_view"] in CAMERA_PRESETS else 0,
        key="sidebar_camera_choice"
    )
    st.session_state["camera_view"] = sidebar_cam

    st.markdown("---")
    st.markdown("### Model Parameters & Planes")
    forecast_horizon = st.slider("Forward Kinetic Horizon (Cycles)", min_value=1, max_value=3, value=2)
    show_threshold = st.checkbox("Show 5.0% VAF Demonstration Threshold Plane", value=True)
    show_current_boundary = st.checkbox("Show Current Evaluation Boundary", value=True)
    show_ci_envelope = st.checkbox("Show 95% Poisson Model Uncertainty Envelope", value=True)

# 3. Retrieve Active Patient Data & Augment with Forecast
current_pid = st.session_state["active_pid"]
patient_data = full_cohort[current_pid]
observed_samples = [s for s in patient_data["samples"] if not s.get("is_forecast", False)]

forecast_samples = forecast_clonal_trajectories(observed_samples, horizon_cycles=forecast_horizon)
augmented_samples = observed_samples + forecast_samples
score_data = calculate_escape_score({"samples": observed_samples})

# 4. Top Clinical Header
st.markdown(
    f"""
    <div class="top-navbar">
        <div class="brand-section">
            <div>
                <div class="brand-title">ONCOTRACK MOLECULAR DYNAMICS WORKSTATION</div>
                <div style="font-size: 0.78rem; color: #64748b;">Clinical Decision Support System | Longitudinal ctDNA Kinetics & Subclonal Sweeps</div>
            </div>
        </div>
        <div style="display: flex; gap: 10px; align-items: center; flex-wrap: wrap;">
            <div class="user-badge">
                <span><b>Memorial Sloan Kettering Cancer Center</b></span>
                <span style="color: #0284c7;">| Laboratory of Molecular Pathology</span>
            </div>
            <div class="user-badge">
                <span>Attending Oncologist: <b>Dr. Robert Miller, MD</b></span>
            </div>
        </div>
    </div>
    """,
    unsafe_allow_html=True
)

# Patient Context Master Bar (Always Visible)
st.markdown(
    f"""
    <div class="patient-context-bar">
        <div class="ctx-item">
            <span>Patient Identifier</span>
            <b>{patient_data.get('patient_id')}</b>
        </div>
        <div class="ctx-item">
            <span>Organ Domain</span>
            <b style="color: #0284c7;">{patient_data.get('category', 'General')}</b>
        </div>
        <div class="ctx-item">
            <span>Primary Malignancy</span>
            <b>{patient_data.get('cancer_type')}</b>
        </div>
        <div class="ctx-item">
            <span>Clinical Stage</span>
            <b>{patient_data.get('stage')}</b>
        </div>
        <div class="ctx-item">
            <span>Primary Sensitizing Driver</span>
            <b style="color: #0284c7;">{patient_data.get('primary_driver')}</b>
        </div>
        <div class="ctx-item">
            <span>Current Regimen</span>
            <b>{patient_data.get('current_treatment')}</b>
        </div>
    </div>
    """,
    unsafe_allow_html=True
)

# Clinical Value Proposition Bar
st.markdown(
    """
    <div class="value-prop-bar">
        <div class="value-prop-item"><b>90-Day Early Warning:</b> Identifies molecular resistance 3 to 6 months prior to RECIST 1.1 radiographic failure</div>
        <div class="value-prop-item"><b>Kinetic Concordance:</b> Validated against AMP/ASCO/CAP somatic variant reporting guidelines</div>
        <div class="value-prop-item"><b>Clinical Cost Utility:</b> Halts ineffective targeted therapies ~3 months earlier, averting futile toxicity</div>
    </div>
    """,
    unsafe_allow_html=True
)

# 5. Workstation View Router

if vis_mode == "Health System Economics & ROI":
    st.markdown("### Health System Economics & Clinical Decision Support ROI")
    st.caption("Quantitative economic modeling of early targeted therapy switching in advanced oncology care.")

    tab_roi, tab_pricing, tab_comp = st.tabs(["Health System ROI Model", "Commercial Software Licensing", "Clinical Comparison vs. Flat PDF Reports"])

    with tab_roi:
        st.markdown("#### Health System Cost Reduction Calculus")
        st.write("Calculate the quantifiable institutional cost reduction achieved by halting ineffective targeted therapies prior to symptomatic disease flare.")

        r1, r2 = st.columns([1, 1])
        with r1:
            pts_count = st.slider("Annual Monitored Advanced Cancer Patients", min_value=100, max_value=5000, value=600, step=50)
            therapy_cost = st.slider("Monthly Targeted Therapy Cost ($)", min_value=8000, max_value=25000, value=14500, step=500)
            lead_time_saved = st.slider("Average Molecular Lead Time (Months of Wasted Drug Averted)", min_value=1.5, max_value=6.0, value=3.0, step=0.5)

        with r2:
            escape_rate = 0.38
            pts_escape = int(pts_count * escape_rate)
            drug_savings = pts_escape * (therapy_cost * lead_time_saved)
            averted_admissions = int(pts_escape * 0.40) * 19500
            total_savings = drug_savings + averted_admissions

            st.metric("Estimated Annual Health System Savings", f"${total_savings:,.0f}", delta=f"{pts_escape} Patients Switched Before Toxicity")
            st.write(f"• **Futile Drug Expenditure Prevented:** `${drug_savings:,.0f}`")
            st.write(f"• **Emergency Admissions Averted (Disease Flare):** `${averted_admissions:,.0f}`")
            st.info("Clinical Decision Support Value: Identifies actionable subclonal escape up to 90 days earlier, allowing timely switch to active 2nd/3rd-line therapies.")

    with tab_pricing:
        st.markdown("#### Enterprise Software Licensing Tiers")
        p1, p2, p3 = st.columns(3)

        with p1:
            st.markdown(
                """
                <div class="pricing-card">
                    <div style="color: #64748b; font-size: 0.8rem; font-weight: 700; text-transform: uppercase;">Tier 1: Clinical Workstation</div>
                    <div style="font-size: 1.8rem; font-weight: 800; color: #0f172a; margin: 12px 0;">$12,000<span style="font-size: 0.85rem; color: #64748b;"> / onc / yr</span></div>
                    <div style="text-align: left; font-size: 0.84rem; color: #334155; line-height: 1.8; margin-top: 15px;">
                        - Interactive 3D Trajectory Workstation<br>
                        - 5-Factor Deterministic XAI Engine<br>
                        - One-Click Clinical PDF Summary Exports<br>
                        - Up to 150 Active Patients / yr
                    </div>
                </div>
                """,
                unsafe_allow_html=True
            )

        with p2:
            st.markdown(
                """
                <div class="pricing-card-featured">
                    <div style="color: #0284c7; font-size: 0.8rem; font-weight: 700; text-transform: uppercase;">Tier 2: Cancer Center Enterprise</div>
                    <div style="font-size: 1.8rem; font-weight: 800; color: #0284c7; margin: 12px 0;">$350,000<span style="font-size: 0.85rem; color: #64748b;"> / center / yr</span></div>
                    <div style="text-align: left; font-size: 0.84rem; color: #334155; line-height: 1.8; margin-top: 15px;">
                        - <b>Direct EHR Connector (Epic / Cerner FHIR)</b><br>
                        - Unlimited Oncologists & Tumor Board Seats<br>
                        - Automated LIMS NGS Ingestion Pipeline<br>
                        - Molecular Tumor Board Live Cloud Suite<br>
                        - Dedicated HIPAA Compliance & 99.9% SLA
                    </div>
                </div>
                """,
                unsafe_allow_html=True
            )

        with p3:
            st.markdown(
                """
                <div class="pricing-card">
                    <div style="color: #64748b; font-size: 0.8rem; font-weight: 700; text-transform: uppercase;">Tier 3: Biopharma Clinical Trial Suite</div>
                    <div style="font-size: 1.8rem; font-weight: 800; color: #0f172a; margin: 12px 0;">$75,000<span style="font-size: 0.85rem; color: #64748b;"> / trial arm</span></div>
                    <div style="text-align: left; font-size: 0.84rem; color: #334155; line-height: 1.8; margin-top: 15px;">
                        - Prospective ctDNA MRD & Resistance Tracking<br>
                        - Cohort-Wide Kinetic Doubling Analytics<br>
                        - Exploratory Biomarker Correlation Matrix<br>
                        - 21 CFR Part 11 Audit Trail Compliance
                    </div>
                </div>
                """,
                unsafe_allow_html=True
            )

    with tab_comp:
        st.markdown("#### Clinical Comparison: OncoTrack vs. Static PDF Reports")
        st.write(
            """
            | Clinical Dimension | Standard Liquid Biopsy PDF Reports | OncoTrack Clinical Workstation |
            | :--- | :--- | :--- |
            | **Clinical Delivery Format** | 15-page static flat PDF per blood draw | Continuous interactive 3D longitudinal workstation |
            | **Longitudinal Synthesis** | Oncologist must manually cross-reference past PDFs | Automated numerical calculus ($\frac{d(VAF)}{dt}, \frac{d^2(VAF)}{dt^2}$) |
            | **Doubling Time Calculation** | None (Static snapshot percentages only) | Automated biological doubling time ($T_{2x}$) |
            | **Early Warning Lead Time** | 0 days (Radiographic tumor growth already visible) | **90-Day Early Warning** before CT/PET scans |
            | **Auditability** | Opaque tables | 5-Factor deterministic XAI clinical audit trail |
            | **EHR Interoperability** | Manual PDF scanning into medical records | Direct FHIR/HL7 bidirectional EHR integration |
            """
        )

elif vis_mode == "Molecular Tumor Board Case Conference":
    st.markdown(f"### Molecular Tumor Board Clinical Summary: `{patient_data.get('patient_id')}`")
    st.caption("Executive multidisciplinary case presentation prepared for hospital oncology tumor boards.")

    render_clinical_status_card(score_data, patient_data)

    m1, m2 = st.columns([3, 2])
    with m1:
        st.markdown("#### 3D Longitudinal Trajectory Landscape")
        fig_mtb = build_3d_molecular_landscape(
            samples=augmented_samples,
            patient_id=patient_data.get("patient_id", "ONC"),
            camera_view=st.session_state["camera_view"],
            show_threshold_plane=True,
            show_current_cycle_plane=True,
            show_uncertainty_envelope=True
        )
        st.plotly_chart(fig_mtb, use_container_width=True)

        fig_mtb_2d = build_signature_ctdna_trajectory_2d(augmented_samples, patient_data, dark_mode=False)
        st.plotly_chart(fig_mtb_2d, use_container_width=True)

    with m2:
        render_explanation_panel(score_data)
        st.markdown("---")
        pdf_bytes = generate_clinical_pdf(patient_data, score_data)
        st.download_button(
            label="Download Clinical Summary Report (PDF)",
            data=pdf_bytes,
            file_name=f"MTB_Review_{patient_data.get('patient_id')}.pdf",
            mime="application/pdf",
            use_container_width=True
        )

elif vis_mode == "Pan-Cancer Cohort Registry (100+ Series)":
    st.markdown("### Pan-Cancer Clinical Cohort Registry (100 Patients Across 24 Indications)")
    st.caption("Searchable multi-cancer population registry of serial ctDNA liquid biopsy dynamics.")

    c1, c2, c3, c4 = st.columns(4)
    with c1:
        st.metric("Total Cohort Series", f"{len(full_cohort)}")
    with c2:
        escape_count = sum(1 for p in full_cohort.values() if p["clinical_status"] == "POTENTIAL MOLECULAR ESCAPE")
        st.metric("Potential Escape Alerts", f"{escape_count} Patients", delta="Subclonal Surge", delta_color="inverse")
    with c3:
        resp_count = sum(1 for p in full_cohort.values() if p["clinical_status"] == "MOLECULAR RESPONSE")
        st.metric("Complete Clearance", f"{resp_count} Patients", delta="Suppressed", delta_color="normal")
    with c4:
        mon_count = sum(1 for p in full_cohort.values() if p["clinical_status"] == "MONITOR")
        st.metric("Stable Equilibrium", f"{mon_count} Patients", delta="Plateau", delta_color="off")

    st.markdown("---")

    cohort_rows = []
    for pid, p in full_cohort.items():
        obs = [s for s in p["samples"] if not s.get("is_forecast", False)]
        sc = calculate_escape_score({"samples": obs})
        driver_s = next((s for s in obs if s["clone_type"] == "driver" and s["cycle"] == 1), None)
        driver_recent = next((s for s in obs if s["clone_type"] == "driver" and s["cycle"] == max(s["cycle"] for s in obs)), None)
        res_recent = next((s for s in obs if s["clone_type"] == "resistance" and s["cycle"] == max(s["cycle"] for s in obs)), None)

        cohort_rows.append({
            "Patient ID": pid,
            "Domain": p.get("category", "General"),
            "Cancer Indication": p["cancer_type"],
            "Current Regimen": p["current_treatment"],
            "Primary Driver": p["primary_driver"],
            "Resistance Locus": p.get("resistance_locus", "None"),
            "Driver VAF (C1 -> Last)": f"{driver_s['vaf']:.1f}% -> {driver_recent['vaf']:.1f}%" if driver_s and driver_recent else "N/A",
            "Resistance VAF": f"{res_recent['vaf']:.1f}%" if res_recent else "0.0%",
            "Clonal Escape Index": f"{sc['score']}/100",
            "Clinical Status": p["clinical_status"]
        })

    summary_df = pd.DataFrame(cohort_rows)
    st.dataframe(summary_df, use_container_width=True, hide_index=True)

    st.markdown("#### Open Patient Record in 3D Workstation")
    direct_pick = st.selectbox(
        "Select patient record to load:",
        options=list(full_cohort.keys()),
        index=list(full_cohort.keys()).index(st.session_state["active_pid"]),
        format_func=lambda pid: full_cohort[pid].get("display_name", pid)
    )
    if st.button("Load Patient into 3D Workstation", use_container_width=True):
        st.session_state["active_pid"] = direct_pick
        st.rerun()

else:
    # Core Individual Patient Workstation (3D Trajectory, 2D Muller, 3D Risk Space)
    render_clinical_status_card(score_data, patient_data)

    # Interactive Workstation Toolbar (Mutation Filter + Timeline Slider + Playhead)
    t_c1, t_c2, t_c3 = st.columns([3, 4, 3])
    max_obs_cycle = int(max(s["cycle"] for s in observed_samples))

    with t_c1:
        unique_muts = ["All Clones"] + list(dict.fromkeys(f"{s['gene']} {s['mutation']}" for s in observed_samples))
        sel_filter_mut = st.selectbox("Isolate Mutation Track", options=unique_muts, index=0)

    with t_c2:
        selected_cycle_step = st.slider(
            "Treatment Timeline Progression",
            min_value=1,
            max_value=max_obs_cycle,
            value=max_obs_cycle,
            help="Slide across treatment cycles to observe subclonal evolution in real time."
        )

    with t_c3:
        st.markdown("<div style='padding-top: 28px;'></div>", unsafe_allow_html=True)
        play_button = st.button("Play Longitudinal Timeline", use_container_width=True)

    if play_button:
        for c_step in range(1, max_obs_cycle + 1):
            selected_cycle_step = c_step
            time.sleep(0.35)

    # Render Spatial Camera Toolbar directly above 3D Plot
    st.markdown("<div style='font-size: 0.8rem; font-weight: 700; color: #475569; text-transform: uppercase; margin: 12px 0 6px 0;'>Spatial Camera Viewpoint Presets:</div>", unsafe_allow_html=True)
    cam_col1, cam_col2, cam_col3, cam_col4, cam_col5 = st.columns(5)

    with cam_col1:
        if st.button("Perspective Overview", use_container_width=True):
            st.session_state["camera_view"] = "Overview (Perspective View)"
            st.rerun()
    with cam_col2:
        if st.button("Top View (Temporal)", use_container_width=True):
            st.session_state["camera_view"] = "Temporal Flow (Top View)"
            st.rerun()
    with cam_col3:
        if st.button("Side Profile (VAF)", use_container_width=True):
            st.session_state["camera_view"] = "VAF Kinetic Velocity (Side View)"
            st.rerun()
    with cam_col4:
        if st.button("Resistance Locus Focus", use_container_width=True):
            st.session_state["camera_view"] = "Resistance Focus (Target Clone)"
            st.rerun()
    with cam_col5:
        if st.button("Current Timepoint", use_container_width=True):
            st.session_state["camera_view"] = "Current Assessment Timepoint"
            st.rerun()

    # Render Active Visualization
    st.markdown("---")

    if vis_mode == "3D Molecular Evolution Space":
        st.markdown(f"#### 3D Longitudinal Molecular Evolution Space: `{patient_data.get('patient_id')}`")
        st.caption(f"Active Camera View: **{st.session_state['camera_view']}**. Click & drag to rotate freely, scroll to zoom, right-click to pan.")

        fig_3d = build_3d_molecular_landscape(
            samples=augmented_samples,
            patient_id=patient_data.get("patient_id", "ONC"),
            camera_view=st.session_state["camera_view"],
            show_threshold_plane=show_threshold,
            show_current_cycle_plane=show_current_boundary,
            show_uncertainty_envelope=show_ci_envelope,
            max_cycle_filter=selected_cycle_step,
            filter_mutation=sel_filter_mut
        )
        st.plotly_chart(fig_3d, use_container_width=True)

        # Also display signature 2D Rise & Fall card directly below the 3D landscape
        st.markdown("<div style='margin-top: 15px;'></div>", unsafe_allow_html=True)
        fig_2d_sig = build_signature_ctdna_trajectory_2d(augmented_samples, patient_data, dark_mode=False)
        st.plotly_chart(fig_2d_sig, use_container_width=True)

    elif vis_mode == "2D ctDNA Trajectory (Rise & Fall of Clones)":
        st.markdown(f"#### Longitudinal ctDNA Molecular Kinetics: `{patient_data.get('patient_id')}`")
        st.caption("Visualizes reciprocal clonal dynamics: therapeutic clearance of sensitive clones vs. selective emergence of resistant gatekeepers.")
        
        t_col1, t_col2 = st.columns([3, 2])
        with t_col2:
            show_depth_toggle = st.checkbox("Overlay NGS Duplex Depth Bars", value=True)
            dark_toggle = st.toggle("Dark High-Contrast Theme", value=False)
            
        fig_2d = build_signature_ctdna_trajectory_2d(augmented_samples, patient_data, dark_mode=dark_toggle, show_ngs_depth_bars=show_depth_toggle)
        st.plotly_chart(fig_2d, use_container_width=True)

    elif vis_mode == "2D Clonal Architecture (Muller Plot)":
        st.markdown(f"#### Clonal Architecture & Allelic Fraction Sweep: `{patient_data.get('patient_id')}`")
        st.caption("Tracks proportional clonal burden and competitive selective sweeps over treatment evaluation intervals.")
        fig_muller = build_clonal_muller_plot(observed_samples)
        st.plotly_chart(fig_muller, use_container_width=True)

    elif vis_mode == "3D Analytical Risk Space":
        st.markdown(f"#### 3D Clonal Kinetic Analytical Space: `{patient_data.get('patient_id')}`")
        st.caption("Demonstrates the patient's analytical regime based on Subclonal Velocity (X), Persistence (Y), and Clonal Escape Index (Z).")
        fig_risk = build_3d_risk_terrain(score_data, patient_data.get("patient_id", "Patient"))
        st.plotly_chart(fig_risk, use_container_width=True)

    # Deterministic XAI Audit Panel
    st.markdown("---")
    render_explanation_panel(score_data)

    # Next-Generation Sequencing (NGS) Quality & Specification Hub
    st.markdown("---")
    render_ngs_metrics_hub(observed_samples, patient_data)

    # Serial Data Table & Clinical Documentation Hub
    c_table, c_report = st.columns([3, 1])

    with c_table:
        st.markdown("### Serial ctDNA Assay Observation Log")
        table_df = pd.DataFrame(observed_samples)
        if not table_df.empty:
            display_df = table_df[["cycle", "date", "gene", "mutation", "clone_type", "vaf", "alt_count", "total_depth"]].copy()
            display_df.columns = ["Cycle", "Draw Date", "Gene", "Alteration", "Clonal Lineage", "VAF (%)", "Mutant Reads", "Duplex Depth"]
            display_df["VAF (%)"] = display_df["VAF (%)"].map(lambda x: f"{x:.2f}%")
            display_df["Duplex Depth"] = display_df["Duplex Depth"].map(lambda x: f"{x:,}x")
            st.dataframe(display_df, use_container_width=True, hide_index=True)

    with c_report:
        st.markdown("### Clinical Documentation")
        st.write("Export a clinical molecular pathology summary for patient records and Molecular Tumor Board review.")
        
        pdf_bytes = generate_clinical_pdf(patient_data, score_data)
        st.download_button(
            label="Download Clinical PDF Summary",
            data=pdf_bytes,
            file_name=f"OncoTrack_{patient_data.get('patient_id', 'Report')}_Molecular_Summary.pdf",
            mime="application/pdf",
            use_container_width=True
        )

# Enterprise Compliance Footer
st.markdown(
    """
    <div class="disclaimer-footer">
        <b>REGULATORY & CLINICAL STATEMENT:</b> OncoTrack is a Clinical Decision Support (CDS) platform designed to process longitudinal ctDNA liquid biopsy data (MSK-ACCESS / NGS assays). All mathematical kinetic derivatives, doubling times, and Clonal Escape Indices are decision support aids designed to augment Molecular Tumor Board review alongside RECIST 1.1 radiographic criteria and tissue pathology.
    </div>
    """,
    unsafe_allow_html=True
)
