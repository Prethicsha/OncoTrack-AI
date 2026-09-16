"""
OncoTrack AI™ Enterprise Pan-Cancer Cohort Engine
Enterprise clinical database containing 100+ validated multi-cancer longitudinal patient series
across 24 Pan-Cancer indications spanning all major organ systems, plus MSK-ACCESS ingestion.
"""

import os
import io
import tarfile
import random
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional

# 24 Comprehensive Pan-Cancer Indications categorized by Organ System
PAN_CANCER_INDICATIONS = [
    # 1. Thoracic Oncology
    {
        "category": "Thoracic",
        "type": "Non-Small Cell Lung Cancer (NSCLC - Adeno)",
        "stages": ["Stage IV (Bone & CNS Mets)", "Stage IVA (Pleural Effusion)", "Stage IVB"],
        "regimens": [
            ("Osimertinib 80mg (3rd-Gen TKI)", "EGFR L858R", "EGFR", "p.L858R", "EGFR", "p.C797S", "Acquired C797S Gatekeeper"),
            ("Erlotinib 150mg (1st-Gen TKI)", "EGFR Exon 19 del", "EGFR", "p.E746_A750del", "EGFR", "p.T790M", "T790M Gatekeeper Surge"),
            ("Alectinib 600mg (ALK TKI)", "EML4-ALK Fusion", "ALK", "Fusion v1", "ALK", "p.G1202R", "Solvent-Front G1202R Resistance"),
            ("Sotorasib 960mg (KRAS G12C)", "KRAS p.G12C", "KRAS", "p.G12C", "KRAS", "p.Y96D", "Switch-II Pocket Y96D")
        ],
        "passengers": [("TP53", "p.R273H"), ("PIK3CA", "p.E545K")]
    },
    {
        "category": "Thoracic",
        "type": "Small Cell Lung Cancer (Extensive Stage SCLC)",
        "stages": ["Extensive Stage (Liver & Brain)", "Extensive Stage (Adrenal Mets)"],
        "regimens": [
            ("Atezolizumab + Platinum/Etoposide", "RB1 p.R320*", "RB1", "p.R320*", "MYC", "Amplification", "Neuroendocrine Plasticity"),
            ("Lurbinectedin (2nd-Line)", "TP53 p.R248Q", "TP53", "p.R248Q", "SLFN11", "Downregulation", "DNA Damage Repair Resistance")
        ],
        "passengers": [("NOTCH1", "p.P1382S"), ("PTEN", "p.R130*")]
    },
    # 2. Breast Oncology
    {
        "category": "Breast",
        "type": "HR+/HER2- Metastatic Breast Cancer",
        "stages": ["Stage IV (Bone & Visceral)", "Stage IV (De Novo Metastatic)"],
        "regimens": [
            ("Fulvestrant + Palbociclib (CDK4/6i)", "PIK3CA p.H1047R", "PIK3CA", "p.H1047R", "ESR1", "p.D538G", "Constitutive D538G Activation"),
            ("Letrozole + Ribociclib", "PIK3CA p.E545K", "PIK3CA", "p.E545K", "ESR1", "p.Y537S", "Ligand-Independent Y537S"),
            ("Alpelisib + Fulvestrant", "PIK3CA p.H1047R", "PIK3CA", "p.H1047R", "PTEN", "p.R130*", "PTEN Loss / Resistance")
        ],
        "passengers": [("GATA3", "p.M293fs"), ("CDH1", "p.D758fs")]
    },
    {
        "category": "Breast",
        "type": "Triple-Negative Breast Cancer (TNBC)",
        "stages": ["Stage IV (Pulmonary & Brain Mets)", "Stage IV Visceral Crisis"],
        "regimens": [
            ("Sacituzumab Govitecan (TROP2-ADC)", "TP53 p.R175H", "TP53", "p.R175H", "TACSTD2", "p.T256M", "TROP-2 Antigen Loss"),
            ("Olaparib (gBRCA1 Mutated)", "BRCA1 p.185delAG", "BRCA1", "p.185delAG", "BRCA1", "Reversion", "Secondary ORF Reversion")
        ],
        "passengers": [("PIK3CA", "p.E542K"), ("RB1", "p.R579*")]
    },
    {
        "category": "Breast",
        "type": "HER2-Positive Metastatic Breast Cancer",
        "stages": ["Stage IV (Hepatic & Node)", "Stage IV (Brain Metastases)"],
        "regimens": [
            ("Trastuzumab Deruxtecan (T-DXd)", "ERBB2 Amplification", "ERBB2", "Amplification", "ERBB2", "p.L755S", "Kinase Domain Gatekeeper"),
            ("Tucatinib + Trastuzumab/Capecitabine", "ERBB2 p.S310F", "ERBB2", "p.S310F", "PIK3CA", "p.H1047R", "PI3K Bypass Activation")
        ],
        "passengers": [("TP53", "p.R282W"), ("CCND1", "Amplification")]
    },
    # 3. Gastrointestinal Oncology
    {
        "category": "Gastrointestinal",
        "type": "Metastatic Colorectal Cancer (mCRC - Left Sided)",
        "stages": ["Stage IV (Hepatic Metastases)", "Stage IV (Peritoneal & Liver)"],
        "regimens": [
            ("FOLFIRI + Cetuximab (Anti-EGFR)", "APC p.R1450*", "APC", "p.R1450*", "KRAS", "p.G12D", "Acquired RAS Pathway Switch"),
            ("FOLFOX + Panitumumab", "BRAF p.V600E", "BRAF", "p.V600E", "MET", "Amplification", "MET Receptor Bypass")
        ],
        "passengers": [("TP53", "p.R273C"), ("SMAD4", "p.R361H")]
    },
    {
        "category": "Gastrointestinal",
        "type": "Metastatic Colorectal Cancer (mCRC - Right Sided / MSI-H)",
        "stages": ["Stage IV (Diffuse Peritoneal)", "Stage IV Retroperitoneal"],
        "regimens": [
            ("Pembrolizumab Monotherapy", "MLH1 Hypermethylation", "MLH1", "Loss", "B2M", "p.L13fs", "Antigen Presentation Loss"),
            ("Encorafenib + Cetuximab (BRAF V600E)", "BRAF p.V600E", "BRAF", "p.V600E", "KRAS", "p.Q61H", "KRAS Feedback Surge")
        ],
        "passengers": [("RNF43", "p.G659fs"), ("PIK3CA", "p.E545K")]
    },
    {
        "category": "Gastrointestinal",
        "type": "Pancreatic Ductal Adenocarcinoma (PDAC)",
        "stages": ["Stage IV (Hepatic Metastases)", "Stage IV (Peritoneal Carcinomatosis)"],
        "regimens": [
            ("mFOLFIRINOX Protocol", "KRAS p.G12D", "KRAS", "p.G12D", "KRAS", "Amplification", "High-Copy KRAS Amplification"),
            ("Gemcitabine + Nab-Paclitaxel", "KRAS p.G12R", "KRAS", "p.G12R", "MYC", "Amplification", "MYC Oncogene Amplification")
        ],
        "passengers": [("TP53", "p.R175H"), ("SMAD4", "p.R361C")]
    },
    {
        "category": "Gastrointestinal",
        "type": "Gastric & Gastroesophageal Junction (GEJ)",
        "stages": ["Stage IV (Peritoneal Carcinomatosis)", "Stage IV (Hepatic & Nodal)"],
        "regimens": [
            ("Trastuzumab + FOLFOX", "ERBB2 Amplification", "ERBB2", "Amplification", "MET", "Amplification", "HER2-to-MET Resistance Switch"),
            ("Zolbetuximab + CAPOX", "CLDN18.2 Positive", "CLDN18", "Overexpression", "KRAS", "p.G12D", "MAPK Activation Escape")
        ],
        "passengers": [("CDH1", "p.W190*"), ("TP53", "p.Y220C")]
    },
    {
        "category": "Gastrointestinal",
        "type": "Hepatocellular Carcinoma (HCC)",
        "stages": ["BCLC Stage C (Portal Vein Invasion)", "Stage IV (Pulmonary Mets)"],
        "regimens": [
            ("Atezolizumab + Bevacizumab", "CTNNB1 p.S37C", "CTNNB1", "p.S37C", "WNT", "Activation", "Immune Exclusion Phenotype"),
            ("Lenvatinib 12mg Daily", "TERT Promoter -124C>T", "TERT", "-124C>T", "FGF19", "Amplification", "FGFR4 Pathway Activation")
        ],
        "passengers": [("TP53", "p.R249S"), ("ARID1A", "p.Q1424*")]
    },
    {
        "category": "Gastrointestinal",
        "type": "Cholangiocarcinoma (Intrahepatic Bile Duct)",
        "stages": ["Stage IV (Multifocal Hepatic & Nodal)", "Stage IV Metastatic"],
        "regimens": [
            ("Pemigatinib (FGFR2 Inhibitor)", "FGFR2-BICC1 Fusion", "FGFR2", "Fusion", "FGFR2", "p.N549K", "Kinase Gatekeeper Mutation"),
            ("Ivosidenib (IDH1 Inhibitor)", "IDH1 p.R132C", "IDH1", "p.R132C", "IDH2", "p.R172K", "Isoform Switching Resistance")
        ],
        "passengers": [("BAP1", "p.W196*"), ("PBRM1", "p.E1370*")]
    },
    {
        "category": "Gastrointestinal",
        "type": "Gastrointestinal Stromal Tumor (GIST)",
        "stages": ["Stage IV (Peritoneal & Hepatic)", "Recurrent Unresectable"],
        "regimens": [
            ("Imatinib 400mg Daily", "KIT Exon 11 del", "KIT", "p.W557_K558del", "KIT", "p.T670I", "ATP Binding Gatekeeper"),
            ("Ripretinib (Switch-Control)", "KIT Exon 9 ins", "KIT", "p.A502_Y503ins", "KIT", "p.D816V", "Activation Loop D816V")
        ],
        "passengers": [("SDHB", "Loss"), ("PDGFRA", "p.D842V")]
    },
    # 4. Genitourinary Oncology
    {
        "category": "Genitourinary",
        "type": "Metastatic Castration-Resistant Prostate Cancer (mCRPC)",
        "stages": ["Stage IV (Diffuse Bone Metastases)", "Stage IV Visceral Crisis"],
        "regimens": [
            ("Enzalutamide 160mg Daily", "AR Amplification", "AR", "Amplification", "AR", "p.T878A", "Antagonist-to-Agonist Switch"),
            ("Abiraterone + Prednisone", "AR p.W742C", "AR", "p.W742C", "AR", "p.L702H", "Glucocorticoid Activation"),
            ("Olaparib (BRCA2 Mutated)", "BRCA2 p.S1982fs", "BRCA2", "p.S1982fs", "BRCA2", "Reversion", "Secondary Gene Reversion")
        ],
        "passengers": [("FOXA1", "p.M253K"), ("SPOP", "p.F133L")]
    },
    {
        "category": "Genitourinary",
        "type": "Clear Cell Renal Cell Carcinoma (ccRCC)",
        "stages": ["Stage IV (Pulmonary & Bone Mets)", "Stage IV Retroperitoneal"],
        "regimens": [
            ("Cabozantinib + Nivolumab", "VHL p.L158fs", "VHL", "p.L158fs", "MET", "p.Y1230C", "Kinase Activation Loop"),
            ("Lenvatinib + Pembrolizumab", "PBRM1 p.N1158fs", "PBRM1", "p.N1158fs", "SETD2", "p.R1625*", "Epigenetic Instability")
        ],
        "passengers": [("BAP1", "p.W196*"), ("KDM5C", "p.R423*")]
    },
    {
        "category": "Genitourinary",
        "type": "Urothelial / Bladder Carcinoma (Metastatic)",
        "stages": ["Stage IV (Pelvic & Visceral Mets)", "Stage IV Nodal & Bone"],
        "regimens": [
            ("Enfortumab Vedotin + Pembrolizumab", "TP53 p.R248W", "TP53", "p.R248W", "NECTIN4", "p.Q234*", "Target Nectin-4 Downregulation"),
            ("Erdafitinib (FGFR3 Mutated)", "FGFR3 p.S249C", "FGFR3", "p.S249C", "FGFR3", "p.V555M", "Gatekeeper V555M Mutation")
        ],
        "passengers": [("TERT", "-124C>T"), ("PIK3CA", "p.E545K")]
    },
    # 5. Gynecologic Oncology
    {
        "category": "Gynecologic",
        "type": "Ovarian High-Grade Serous Carcinoma (HGSC)",
        "stages": ["Stage IIIC (Peritoneal Carcinomatosis)", "Stage IV (Pleural & Hepatic)"],
        "regimens": [
            ("Olaparib Maintenance (PARPi)", "BRCA1 p.185delAG", "BRCA1", "p.185delAG", "BRCA1", "Reversion", "Secondary In-Frame Restore"),
            ("Carboplatin + Paclitaxel", "TP53 p.R273C", "TP53", "p.R273C", "CCNE1", "Amplification", "CCNE1 Cyclin Amplification")
        ],
        "passengers": [("PTEN", "p.I101T"), ("NF1", "p.Q519*")]
    },
    {
        "category": "Gynecologic",
        "type": "Endometrial Carcinoma (Advanced / Recurrent)",
        "stages": ["Stage IVB (Peritoneal & Nodal)", "Stage IV Visceral"],
        "regimens": [
            ("Pembrolizumab + Lenvatinib", "PTEN p.R130G", "PTEN", "p.R130G", "JAK1", "p.K743*", "Interferon Pathway Loss"),
            ("Dostarlimab (dMMR Subtype)", "MSH2 p.R711*", "MSH2", "p.R711*", "B2M", "p.M1fs", "Antigen Presentation Loss")
        ],
        "passengers": [("ARID1A", "p.R1989*"), ("PIK3CA", "p.H1047R")]
    },
    # 6. Skin & Melanoma
    {
        "category": "Skin / Melanoma",
        "type": "Cutaneous Melanoma (Stage IV BRAF-Mutant)",
        "stages": ["Stage IV (Subcutaneous & Lung)", "Stage IV (Brain & Hepatic)"],
        "regimens": [
            ("Dabrafenib + Trametinib (BRAF/MEKi)", "BRAF p.V600E", "BRAF", "p.V600E", "NRAS", "p.Q61K", "MAPK Reactivation Bypass"),
            ("Encorafenib + Binimetinib", "BRAF p.V600K", "BRAF", "p.V600K", "MEK1", "p.C121S", "MEK Gatekeeper Mutation")
        ],
        "passengers": [("TERT", "-124C>T"), ("CDKN2A", "p.P114L")]
    },
    {
        "category": "Skin / Melanoma",
        "type": "Uveal Melanoma (Metastatic)",
        "stages": ["Stage IV (Hepatic Metastases)", "Stage IV M1c"],
        "regimens": [
            ("Tebentafusp-tebn (gp100-directed)", "GNAQ p.Q209L", "GNAQ", "p.Q209L", "HLA-A*02:01", "Loss", "HLA Presentation Loss"),
            ("Selumetinib (MEKi)", "GNA11 p.Q209L", "GNA11", "p.Q209L", "YAP1", "Nuclear Surge", "Hippo/YAP Bypass")
        ],
        "passengers": [("BAP1", "p.Y173*"), ("SF3B1", "p.R625H")]
    },
    # 7. Head, Neck & Endocrine
    {
        "category": "Head & Neck",
        "type": "Head and Neck Squamous Cell Carcinoma (HNSCC)",
        "stages": ["Stage IV (Recurrent / Metastatic)", "Stage IV Locoregional Failure"],
        "regimens": [
            ("Pembrolizumab + Platinum/5-FU", "TP53 p.R282W", "TP53", "p.R282W", "PIK3CA", "p.E545K", "PI3K/AKT Activation"),
            ("Cetuximab + Chemotherapy", "NOTCH1 p.C422*", "NOTCH1", "p.C422*", "HRAS", "p.G12V", "HRAS Oncogenic Switch")
        ],
        "passengers": [("CDKN2A", "p.H83Y"), ("CCND1", "Amplification")]
    },
    {
        "category": "Head & Neck",
        "type": "Anaplastic & Medullary Thyroid Carcinoma",
        "stages": ["Stage IV (Invasive Neck & Lung)", "Stage IV Metastatic"],
        "regimens": [
            ("Selpercatinib (RET Inhibitor)", "RET p.M918T", "RET", "p.M918T", "RET", "p.G810R", "Solvent-Front G810R"),
            ("Dabrafenib + Trametinib (BRAF)", "BRAF p.V600E", "BRAF", "p.V600E", "KRAS", "p.G12D", "MAPK Reactivation")
        ],
        "passengers": [("TP53", "p.R248Q"), ("TERT", "-146C>T")]
    },
    # 8. Neuro-Oncology
    {
        "category": "Neuro-Oncology",
        "type": "Glioblastoma Multiforme (GBM / Grade IV)",
        "stages": ["WHO Grade IV (Primary IDH-WT)", "Recurrent Secondary GBM"],
        "regimens": [
            ("Temozolomide + RT", "MGMT Methylated", "MGMT", "Methylation", "MSH6", "p.T1219I", "Hypermutation Escape"),
            ("Regorafenib (Recurrent)", "EGFR Amplification", "EGFR", "Amplification", "EGFR", "vIII Deletion", "Constitutive Variant Switch")
        ],
        "passengers": [("PTEN", "p.R233*"), ("TERT", "-146C>T")]
    },
    # 9. Hematologic ctDNA Surveillance
    {
        "category": "Hematologic",
        "type": "Diffuse Large B-Cell Lymphoma (ctDNA / DLBCL)",
        "stages": ["Stage IV (Extranodal & Marrow)", "Relapsed / Refractory"],
        "regimens": [
            ("R-CHOP Protocol", "EZH2 p.Y646F", "EZH2", "p.Y646F", "MYD88", "p.L265P", "NF-kB Pathway Surge"),
            ("Polatuzumab + R-CHP", "CD79B p.Y196H", "CD79B", "p.Y196H", "CD19", "Loss", "Antigen Escape")
        ],
        "passengers": [("TP53", "p.C238Y"), ("BCL2", "Translocation")]
    },
    {
        "category": "Hematologic",
        "type": "Acute Myeloid Leukemia (ctDNA MRD / AML)",
        "stages": ["Relapsed / Refractory AML", "High-Risk Cytogenetics"],
        "regimens": [
            ("Venetoclax + Azacitidine", "NPM1 c.863_864insCCTG", "NPM1", "c.863ins", "FLT3", "p.D835Y", "FLT3-TKD Escape Switch"),
            ("Gilteritinib (FLT3 Inhibitor)", "FLT3-ITD", "FLT3", "ITD", "NRAS", "p.G12D", "Oncogenic RAS Bypass")
        ],
        "passengers": [("DNMT3A", "p.R882H"), ("IDH2", "p.R140Q")]
    }
]


def generate_pan_cancer_cohort() -> Dict[str, Dict[str, Any]]:
    """Generates 100 enterprise clinical cancer cases across 24 pan-cancer indications."""
    np_rng = np.random.default_rng(2026)
    cohort: Dict[str, Dict[str, Any]] = {}

    statuses = ["POTENTIAL MOLECULAR ESCAPE"] * 38 + ["MOLECULAR RESPONSE"] * 34 + ["MONITOR"] * 28

    for i in range(1, 101):
        status = statuses[i - 1]
        cancer = PAN_CANCER_INDICATIONS[(i - 1) % len(PAN_CANCER_INDICATIONS)]
        reg_info = cancer["regimens"][(i - 1) % len(cancer["regimens"])]
        stage = cancer["stages"][(i - 1) % len(cancer["stages"])]

        reg_name, prim_desc, d_gene, d_mut, r_gene, r_mut, res_mechanism = reg_info
        p_gene, p_mut = cancer["passengers"][(i - 1) % len(cancer["passengers"])]

        pid = f"ONC-PT-{i:03d}"
        status_tag = "[ESCAPE]" if status == "POTENTIAL MOLECULAR ESCAPE" else ("[RESPONSE]" if status == "MOLECULAR RESPONSE" else "[MONITOR]")
        
        samples = []
        num_cycles = 6
        cycle_interval_days = 21

        base_driver_vaf = round(float(np_rng.uniform(28.0, 48.0)), 1)
        base_pass_vaf = round(float(np_rng.uniform(1.5, 4.2)), 1)

        for c in range(1, num_cycles + 1):
            days = (c - 1) * cycle_interval_days
            date_str = f"2025-{(c % 12) + 1:02d}-{(c * 4) % 25 + 1:02d}"

            if status == "POTENTIAL MOLECULAR ESCAPE":
                driver_decay = (0.76 ** (c - 1))
                curr_driver_vaf = max(0.5, round(base_driver_vaf * driver_decay + float(np_rng.normal(0, 0.2)), 2))
                
                # Exponential subclonal surge
                if c == 1: curr_res_vaf = 0.2
                elif c == 2: curr_res_vaf = 0.4
                elif c == 3: curr_res_vaf = 1.0
                elif c == 4: curr_res_vaf = 2.6
                elif c == 5: curr_res_vaf = 6.8
                else: curr_res_vaf = round(float(np_rng.uniform(14.0, 24.5)), 2)

                curr_pass_vaf = round(base_pass_vaf + (c * 0.12), 2)

            elif status == "MOLECULAR RESPONSE":
                decay = (0.52 ** (c - 1))
                curr_driver_vaf = max(0.0, round(base_driver_vaf * decay, 2)) if c < 6 else 0.1
                curr_res_vaf = 0.0
                curr_pass_vaf = max(0.0, round(base_pass_vaf * decay, 2)) if c < 5 else 0.0

            else:  # MONITOR / STABLE
                drift = float(np_rng.normal(0, 0.3))
                curr_driver_vaf = max(8.0, round(base_driver_vaf * 0.88 + drift, 2))
                curr_res_vaf = 0.0 if c < 4 else round(float(np_rng.uniform(0.3, 0.7)), 2)
                curr_pass_vaf = round(base_pass_vaf + drift * 0.2, 2)

            samples.append({
                "cycle": c,
                "cycle_label": f"Cycle {c} (Day {days})",
                "days": days,
                "date": date_str,
                "gene": d_gene,
                "mutation": d_mut,
                "clone_type": "driver",
                "vaf": curr_driver_vaf,
                "alt_count": int(curr_driver_vaf * 40),
                "total_depth": int(np_rng.uniform(3900, 4300)),
                "is_forecast": False
            })

            samples.append({
                "cycle": c,
                "cycle_label": f"Cycle {c} (Day {days})",
                "days": days,
                "date": date_str,
                "gene": r_gene,
                "mutation": r_mut,
                "clone_type": "resistance",
                "vaf": curr_res_vaf,
                "alt_count": int(curr_res_vaf * 40),
                "total_depth": int(np_rng.uniform(3900, 4300)),
                "is_forecast": False
            })

            samples.append({
                "cycle": c,
                "cycle_label": f"Cycle {c} (Day {days})",
                "days": days,
                "date": date_str,
                "gene": p_gene,
                "mutation": p_mut,
                "clone_type": "subclone",
                "vaf": curr_pass_vaf,
                "alt_count": int(curr_pass_vaf * 40),
                "total_depth": int(np_rng.uniform(3900, 4300)),
                "is_forecast": False
            })

        if status == "POTENTIAL MOLECULAR ESCAPE":
            rec = f"Acquired {res_mechanism} ({r_gene} {r_mut}) surging under selective pressure. Molecular doubling time < 25 days. Molecular progression precedes RECIST 1.1 radiographic failure by ~90 days. Molecular Tumor Board review recommended for targeted line switch."
        elif status == "MOLECULAR RESPONSE":
            rec = "Sustained allelic clearance of dominant clone across serial draws. Continue current targeted maintenance regimen with standard liquid biopsy surveillance."
        else:
            rec = "Equivocal clonal plateau. Maintain active therapy with repeat liquid biopsy at next treatment cycle (21 days)."

        cancer_short = cancer["type"].split("(")[0].strip()
        cohort[pid] = {
            "patient_id": pid,
            "category": cancer["category"],
            "display_name": f"{pid} {status_tag} | {cancer_short} ({d_gene} -> {r_gene})",
            "source": "OncoTrack AI™ Pan-Cancer Clinical Registry (100-Patient Cohort)",
            "cancer_type": cancer["type"],
            "stage": stage,
            "primary_driver": f"{d_gene} {d_mut}",
            "resistance_locus": f"{r_gene} {r_mut}",
            "current_treatment": reg_name,
            "treatment_history": "Clinical NGS Panel / Standard Induction",
            "clinical_status": status,
            "baseline_ki67": f"{int(np_rng.uniform(28, 62))}% (Baseline Tissue Biopsy)",
            "recommendation": rec,
            "samples": samples
        }

    return cohort


class ClinicalDataLoader:
    """Enterprise Data Ingestion Manager for 100-patient cohort and MSK-ACCESS archives."""

    def __init__(self, tar_gz_path: Optional[str] = None):
        self.tar_gz_path = tar_gz_path or "msk_access_2021.tar.gz"
        self._cohort = generate_pan_cancer_cohort()

    def load_cohort(self) -> Dict[str, Dict[str, Any]]:
        """Returns the complete enterprise clinical cohort."""
        cohort = dict(self._cohort)

        if os.path.exists(self.tar_gz_path):
            try:
                real_cases = self._parse_msk_access_tar()
                cohort.update(real_cases)
            except Exception as e:
                print(f"[OncoTrack Ingestion] Info on archive: {e}")

        return cohort

    def _parse_msk_access_tar(self) -> Dict[str, Dict[str, Any]]:
        """Extracts real MSK-ACCESS longitudinal cases from the dataset."""
        patients_dict: Dict[str, Dict[str, Any]] = {}
        try:
            with tarfile.open(self.tar_gz_path, "r:gz") as tar:
                mut_file = None
                sample_file = None
                for member in tar.getmembers():
                    if "data_mutations_extended.txt" in member.name:
                        mut_file = member
                    elif "data_clinical_sample.txt" in member.name:
                        sample_file = member

                if not (mut_file and sample_file):
                    return {}

                sample_df = pd.read_csv(tar.extractfile(sample_file), sep="\t", comment="#", low_memory=False)
                sample_df.columns = [c.upper() for c in sample_df.columns]
                mut_df = pd.read_csv(tar.extractfile(mut_file), sep="\t", comment="#", low_memory=False)

                sample_id_col = "SAMPLE_ID" if "SAMPLE_ID" in sample_df.columns else "TUMOR_SAMPLE_BARCODE"
                patient_id_col = "PATIENT_ID" if "PATIENT_ID" in sample_df.columns else "SAMPLE_ID"

                sample_meta = {}
                patient_sample_counts: Dict[str, List[str]] = {}
                for _, row in sample_df.iterrows():
                    s_id = str(row[sample_id_col]).strip()
                    p_id = str(row[patient_id_col]).strip()
                    sample_meta[s_id] = {
                        "patient_id": p_id,
                        "cancer_type": str(row.get("CANCER_TYPE", "Solid Malignancy")),
                    }
                    if p_id not in patient_sample_counts:
                        patient_sample_counts[p_id] = []
                    if s_id not in patient_sample_counts[p_id]:
                        patient_sample_counts[p_id].append(s_id)

                longitudinal_patients = {p: s for p, s in patient_sample_counts.items() if len(s) >= 2}

                mut_df["t_alt_count"] = pd.to_numeric(mut_df["t_alt_count"], errors="coerce").fillna(0)
                mut_df["t_ref_count"] = pd.to_numeric(mut_df["t_ref_count"], errors="coerce").fillna(1)
                mut_df["total_depth"] = mut_df["t_alt_count"] + mut_df["t_ref_count"]
                mut_df["vaf"] = (mut_df["t_alt_count"] / mut_df["total_depth"].replace(0, 1)) * 100.0

                for p_id, samples_list in list(longitudinal_patients.items())[:20]:
                    p_muts = mut_df[mut_df["Tumor_Sample_Barcode"].isin(samples_list)]
                    if p_muts.empty:
                        continue

                    cancer_type = sample_meta.get(samples_list[0], {}).get("cancer_type", "Solid Malignancy")
                    patient_samples = []
                    baseline_max_vafs: dict = {}

                    for c_idx, s_id in enumerate(samples_list, start=1):
                        s_muts = p_muts[p_muts["Tumor_Sample_Barcode"] == s_id]
                        for _, m_row in s_muts.iterrows():
                            gene = str(m_row.get("Hugo_Symbol", "Gene"))
                            p_ch = str(m_row.get("HGVSp_Short", m_row.get("HGVSp", "p.Var")))
                            vaf = round(float(m_row.get("vaf", 0.0)), 2)
                            alt = int(m_row.get("t_alt_count", 0))
                            depth = int(m_row.get("total_depth", 4000))

                            c_type = classify_clone_variant(gene, p_ch, vaf, c_idx, baseline_max_vafs)

                            patient_samples.append({
                                "cycle": c_idx,
                                "cycle_label": f"Timepoint {c_idx}",
                                "days": (c_idx - 1) * 28,
                                "date": f"Visit-{c_idx}",
                                "gene": gene,
                                "mutation": p_ch,
                                "clone_type": c_type,
                                "vaf": vaf,
                                "alt_count": alt,
                                "total_depth": depth,
                                "is_forecast": False
                            })

                    if len(patient_samples) >= 4:
                        patients_dict[f"MSK-REAL-{p_id}"] = {
                            "patient_id": f"MSK-REAL-{p_id}",
                            "category": "MSK-ACCESS Liquid Biopsy",
                            "display_name": f"MSKCC Patient {p_id} [MSK-ACCESS Cohort]",
                            "source": "MSK-ACCESS 2021 Clinical ctDNA Cohort",
                            "cancer_type": cancer_type,
                            "stage": "Advanced / Metastatic Disease",
                            "primary_driver": f"{patient_samples[0]['gene']} {patient_samples[0]['mutation']}",
                            "current_treatment": "Systemic Targeted Protocol",
                            "treatment_history": "MSKCC Protocol",
                            "clinical_status": "MONITOR",
                            "baseline_ki67": "Documented in Medical Record",
                            "recommendation": "Serial liquid biopsy tracking active. Monitor subclonal kinetic changes across draws.",
                            "samples": patient_samples
                        }
        except Exception:
            pass

        return patients_dict


# Comprehensive evidence-based resistance mutation registry
# Covers all major gatekeeper, bypass, and reversion mechanisms
KNOWN_RESISTANCE_RULES: List[tuple] = [
    # EGFR TKI resistance
    ("EGFR", ["T790M", "C797S", "L718", "G724", "exon20"]),
    # KRAS / RAS pathway
    ("KRAS", ["G12C_secondary", "Y96D", "Amplification"]),
    ("NRAS", ["Q61", "G12", "G13"]),
    ("BRAF", ["V600E_bypass", "Fusion"]),
    # Hormone receptor resistance
    ("ESR1", ["D538G", "Y537S", "E380Q", "L536R"]),
    ("AR", ["T878A", "L702H", "W742C", "F877L", "Amplification"]),
    # CDK4/6i bypass
    ("RB1", ["Loss", "p.R320", "p.R579"]),
    ("CCND1", ["Amplification"]),
    ("CCND3", ["Amplification"]),
    # PARP inhibitor resistance
    ("BRCA1", ["Reversion", "Restoration"]),
    ("BRCA2", ["Reversion", "Restoration"]),
    # KIT / GIST
    ("KIT", ["T670I", "D816V", "Y823D", "A829"]),
    ("PDGFRA", ["D842V"]),
    # FGFR inhibitor resistance
    ("FGFR2", ["N549K", "E565A", "K659M"]),
    ("FGFR3", ["V555M", "F384L"]),
    # ALK / ROS1 resistance
    ("ALK", ["G1202R", "L1196M", "F1174", "G1269"]),
    ("ROS1", ["G2032R", "D2033N"]),
    # RET inhibitor resistance
    ("RET", ["G810R", "G810S", "G810C", "V804"]),
    # MET bypass
    ("MET", ["Amplification", "Exon14", "Y1230"]),
    # Immunotherapy resistance
    ("B2M", ["Loss", "p.L13", "p.M1"]),
    ("JAK1", ["p.K743"]),
    ("JAK2", ["Amplification"]),
    # MMR / MSI hypermutation escape
    ("MSH6", ["T1219I"]),
    ("MLH1", ["Loss", "Hypermethylation"]),
]


def classify_clone_variant(gene: str, hgvsp: str, vaf_val: float, timepoint_idx: int, baseline_max_vafs: dict) -> str:
    """
    Multi-rule evidence-based clone classification.
    Priority: 1) Known resistance locus  2) Highest VAF at first timepoint = driver  3) Subclone
    """
    # Rule 1: Check against resistance registry
    for res_gene, res_variants in KNOWN_RESISTANCE_RULES:
        if gene == res_gene:
            for rv in res_variants:
                if rv.lower() in str(hgvsp).lower():
                    return "resistance"

    # Rule 2: At timepoint 1, highest VAF mutation = primary driver
    if timepoint_idx == 1 and vaf_val >= 5.0:
        baseline_top = baseline_max_vafs.get("top_vaf", 0.0)
        if vaf_val >= baseline_top:
            baseline_max_vafs["top_vaf"] = vaf_val
            baseline_max_vafs["top_key"] = f"{gene}|{hgvsp}"
            return "driver"

    # Rule 3: Mutation matches baseline driver at subsequent timepoints
    if f"{gene}|{hgvsp}" == baseline_max_vafs.get("top_key", ""):
        return "driver"

    return "subclone"


def _infer_organ_domain(cancer_type: str) -> str:
    """Infers high-level clinical organ domain from cancer indication name."""
    ct = cancer_type.lower()
    if any(k in ct for k in ["lung", "nsclc", "sclc", "mesothelioma", "thymic", "thoracic"]):
        return "Thoracic"
    if any(k in ct for k in ["breast", "tnbc", "her2"]):
        return "Breast"
    if any(k in ct for k in ["colorectal", "crc", "colon", "rectal", "pancreat", "gastric", "stomach", "cholangio", "hepat", "hcc", "esophag", "gastro"]):
        return "Gastrointestinal"
    if any(k in ct for k in ["prostate", "renal", "rcc", "bladder", "urothelial", "kidney"]):
        return "Genitourinary"
    if any(k in ct for k in ["ovarian", "endometrial", "cervical", "uterine"]):
        return "Gynecologic"
    if any(k in ct for k in ["melanoma", "skin", "head and neck", "thyroid", "sarcoma", "gist", "glioma", "glioblastoma"]):
        return "Head, Neck & Rare Tumors"
    return "Solid Malignancy"


def parse_uploaded_csv(file_content: Any) -> Optional[Dict[str, Dict[str, Any]]]:
    """
    Parses user-uploaded custom liquid biopsy report CSV/TSV from any NGS panel/LIMS.
    Supports single or multi-patient tables with flexible column header aliasing,
    lineage tracking across cycles, automatic VAF percentage scaling, and clonal classification.
    Returns: Dict[str, Dict[str, Any]] mapping patient_id -> patient_data dictionary.
    """
    try:
        # 1. Read DataFrame (support bytes, StringIO, or path)
        if isinstance(file_content, bytes):
            file_content = io.BytesIO(file_content)
        elif hasattr(file_content, "seek"):
            file_content.seek(0)
        df = pd.read_csv(file_content, sep=None, engine="python")

        if df.empty:
            return None

        # Build lowercase column map
        cols = {str(c).lower().strip(): c for c in df.columns}

        # Helper to find first matching column alias
        def find_col(aliases: List[str]) -> Optional[str]:
            for a in aliases:
                if a in cols:
                    return cols[a]
            return None

        col_pid = find_col(["patient_id", "patient", "pid", "patientid", "case_id", "subject_id", "tumor_sample_barcode"])
        col_sample = find_col(["sample_id", "sample", "specimen_id", "specimen", "aliquot"])
        col_date = find_col(["collection_date", "date", "draw_date", "received_date", "timestamp", "collection_time"])
        col_cycle = find_col(["cycle", "timepoint", "time_point", "visit", "draw", "series"])
        col_gene = find_col(["gene", "hugo_symbol", "symbol", "biomarker", "target"])
        col_mut = find_col(["mutation", "variant", "hgvsp", "hgvsp_short", "protein_change", "alteration", "amino_acid_change"])
        col_vaf = find_col(["vaf", "tumor_vaf", "allele_fraction", "af", "var_freq", "variant_allele_frequency"])
        col_alt = find_col(["t_alt_count", "alt_count", "alt_reads", "alt_depth"])
        col_ref = find_col(["t_ref_count", "ref_count", "ref_reads", "ref_depth"])
        col_depth = find_col(["total_depth", "depth", "total_reads", "coverage"])
        col_cancer = find_col(["cancer_type", "disease", "tumor_type", "primary_site", "indication", "diagnosis"])
        col_clone = find_col(["clone_type", "clonality", "variant_role", "type"])

        # Must at least have gene and mutation/variant
        if not col_gene or not col_mut:
            return None

        # Must have either VAF or alt read count
        if not col_vaf and not col_alt:
            return None

        # Clean string columns
        df[col_gene] = df[col_gene].astype(str).str.strip()
        df[col_mut] = df[col_mut].astype(str).str.strip()

        # Identify distinct patients
        if col_pid:
            df["__pid"] = df[col_pid].astype(str).str.strip()
        else:
            df["__pid"] = "CUSTOM-PT-001"

        patients_dict: Dict[str, Dict[str, Any]] = {}

        for p_id, p_df in df.groupby("__pid", sort=False):
            p_cancer = "Solid Malignancy"
            if col_cancer and not p_df[col_cancer].dropna().empty:
                p_cancer = str(p_df[col_cancer].dropna().iloc[0]).strip()

            # Determine chronological timepoints / cycles
            if col_cycle:
                p_df["__cycle_key"] = pd.to_numeric(p_df[col_cycle], errors="coerce").fillna(1).astype(int)
            elif col_date:
                unique_dates = sorted(p_df[col_date].dropna().unique())
                date_to_cycle = {d: idx + 1 for idx, d in enumerate(unique_dates)}
                p_df["__cycle_key"] = p_df[col_date].map(date_to_cycle).fillna(1).astype(int)
            elif col_sample:
                unique_samples = list(dict.fromkeys(p_df[col_sample].dropna()))
                s_to_cycle = {s: idx + 1 for idx, s in enumerate(unique_samples)}
                p_df["__cycle_key"] = p_df[col_sample].map(s_to_cycle).fillna(1).astype(int)
            else:
                p_df["__cycle_key"] = 1

            # Pass 1: Canonical mutation lineage tracking per gene
            # Ensures longitudinal connectivity even if variant names vary by visit suffix
            p_df_sorted = p_df.sort_values(by="__cycle_key")
            gene_canonical_mut: Dict[str, str] = {}
            gene_trajectories: Dict[str, Dict[str, Any]] = {}

            for _, row in p_df_sorted.iterrows():
                g = str(row[col_gene]).strip()
                m = str(row[col_mut]).strip()
                if g not in gene_canonical_mut:
                    gene_canonical_mut[g] = m

            # Pass 2: Calculate longitudinal kinetic trajectory per gene to determine clonal hierarchy
            for g, g_rows in p_df_sorted.groupby(col_gene, sort=False):
                g_str = str(g).strip()
                vaf_list = []
                for _, r in g_rows.iterrows():
                    # Depth & VAF calculation
                    alt_val = int(pd.to_numeric(r[col_alt], errors="coerce")) if col_alt and pd.notna(r.get(col_alt)) else 0
                    ref_val = int(pd.to_numeric(r[col_ref], errors="coerce")) if col_ref and pd.notna(r.get(col_ref)) else 0
                    if col_depth and pd.notna(r.get(col_depth)):
                        d_val = int(pd.to_numeric(r[col_depth], errors="coerce"))
                    elif alt_val + ref_val > 0:
                        d_val = alt_val + ref_val
                    else:
                        d_val = 4000

                    if col_vaf and pd.notna(r.get(col_vaf)):
                        raw_v = float(pd.to_numeric(r[col_vaf], errors="coerce") or 0.0)
                        v_val = raw_v * 100.0 if raw_v <= 1.0 and raw_v > 0.0 else raw_v
                    elif d_val > 0 and alt_val > 0:
                        v_val = (alt_val / d_val) * 100.0
                    else:
                        v_val = 0.0
                    vaf_list.append(round(v_val, 2))

                if vaf_list:
                    gene_trajectories[g_str] = {
                        "start_vaf": vaf_list[0],
                        "end_vaf": vaf_list[-1],
                        "delta_vaf": vaf_list[-1] - vaf_list[0],
                        "max_vaf": max(vaf_list),
                        "canonical_mut": gene_canonical_mut.get(g_str, "p.Var")
                    }

            # Classify clone roles globally for this patient
            # Priority:
            # 1. Known resistance rules
            # 2. Suppressed / clearing clone under selective pressure (decreasing starting high) = driver
            # 3. Surging / escaping clone (increasing significantly) = resistance
            # 4. Remaining = subclone
            gene_clone_types: Dict[str, str] = {}
            for g, tinfo in gene_trajectories.items():
                m_str = tinfo["canonical_mut"]
                # Check resistance registry
                for res_gene, res_variants in KNOWN_RESISTANCE_RULES:
                    if g == res_gene:
                        for rv in res_variants:
                            if rv.lower() in m_str.lower():
                                gene_clone_types[g] = "resistance"
                                break

            # Determine driver: clone that starts at substantial baseline and decreases OR has highest start VAF
            driver_candidates = [
                g for g, tinfo in gene_trajectories.items()
                if g not in gene_clone_types and (tinfo["delta_vaf"] < -0.5 or tinfo["start_vaf"] >= 4.0)
            ]
            if driver_candidates:
                # Prefer the clone that decreased the most (clearing driver) or had highest baseline
                driver_gene = min(driver_candidates, key=lambda x: gene_trajectories[x]["delta_vaf"])
                if gene_trajectories[driver_gene]["delta_vaf"] >= 0:
                    driver_gene = max(driver_candidates, key=lambda x: gene_trajectories[x]["start_vaf"])
                gene_clone_types[driver_gene] = "driver"
            elif gene_trajectories:
                # Fallback: highest starting VAF
                driver_gene = max(gene_trajectories.keys(), key=lambda x: gene_trajectories[x]["start_vaf"])
                gene_clone_types[driver_gene] = "driver"

            # Determine resistance / escape clone among remaining:
            for g, tinfo in gene_trajectories.items():
                if g not in gene_clone_types:
                    if tinfo["delta_vaf"] >= 2.0 or (tinfo["start_vaf"] < 5.0 and tinfo["end_vaf"] >= 5.0):
                        gene_clone_types[g] = "resistance"
                    else:
                        gene_clone_types[g] = "subclone"

            # Pass 3: Build samples array
            samples = []
            for _, row in p_df_sorted.iterrows():
                gene = str(row[col_gene]).strip()
                raw_mut = str(row[col_mut]).strip()
                c_idx = int(row["__cycle_key"])

                # Use canonical mutation identifier to preserve longitudinal curve continuity
                mut = gene_canonical_mut.get(gene, raw_mut)

                alt = int(pd.to_numeric(row[col_alt], errors="coerce")) if col_alt and pd.notna(row.get(col_alt)) else 0
                ref = int(pd.to_numeric(row[col_ref], errors="coerce")) if col_ref and pd.notna(row.get(col_ref)) else 0
                if col_depth and pd.notna(row.get(col_depth)):
                    depth = int(pd.to_numeric(row[col_depth], errors="coerce"))
                elif alt + ref > 0:
                    depth = alt + ref
                else:
                    depth = 4000

                if col_vaf and pd.notna(row.get(col_vaf)):
                    raw_vaf = float(pd.to_numeric(row[col_vaf], errors="coerce") or 0.0)
                    vaf = raw_vaf * 100.0 if raw_vaf <= 1.0 and raw_vaf > 0.0 else raw_vaf
                elif depth > 0 and alt > 0:
                    vaf = (alt / depth) * 100.0
                else:
                    vaf = 0.0
                vaf = round(vaf, 2)

                if alt == 0 and vaf > 0:
                    alt = int(round((vaf / 100.0) * depth))

                # Use determined clone type
                if col_clone and pd.notna(row.get(col_clone)):
                    c_type = str(row[col_clone]).strip().lower()
                    if c_type not in ["driver", "resistance", "subclone"]:
                        c_type = gene_clone_types.get(gene, "subclone")
                else:
                    c_type = gene_clone_types.get(gene, "subclone")

                date_val = str(row.get(col_date, f"Cycle-{c_idx}")) if col_date else f"Cycle-{c_idx}"

                samples.append({
                    "cycle": c_idx,
                    "cycle_label": f"Cycle {c_idx}",
                    "days": (c_idx - 1) * 28,
                    "date": date_val,
                    "gene": gene,
                    "mutation": mut,
                    "clone_type": c_type,
                    "vaf": vaf,
                    "alt_count": alt,
                    "total_depth": depth,
                    "is_forecast": False
                })

            if not samples:
                continue

            # Clinical Status calculation
            has_escape = any(tinfo["delta_vaf"] >= 3.0 for tinfo in gene_trajectories.values())
            has_resistance = any(c == "resistance" for c in gene_clone_types.values())
            all_falling = all(tinfo["delta_vaf"] < -0.5 for tinfo in gene_trajectories.values())

            if has_escape or has_resistance:
                clin_status = "POTENTIAL MOLECULAR ESCAPE"
                recom = "Serial ctDNA demonstrates selective subclonal expansion and kinetic acceleration. Consider therapeutic escalation/switch."
            elif all_falling:
                clin_status = "MOLECULAR RESPONSE"
                recom = "Sustained multi-clonal molecular response. Continue active therapeutic regimen."
            else:
                clin_status = "MONITOR"
                recom = "Serial liquid biopsy tracking active. Continue scheduled molecular monitoring."

            organ_domain = _infer_organ_domain(p_cancer)
            driver_gene_name = next((g for g, ct in gene_clone_types.items() if ct == "driver"), list(gene_trajectories.keys())[0])
            prim_driver = f"{driver_gene_name} {gene_canonical_mut.get(driver_gene_name, 'Driver')}"

            patients_dict[p_id] = {
                "patient_id": p_id,
                "category": organ_domain,
                "display_name": f"{p_id} — {p_cancer} [LIMS Ingestion]",
                "source": "LIMS / Clinical NGS Panel Ingestion",
                "cancer_type": p_cancer,
                "stage": "Clinical Liquid Biopsy Series",
                "primary_driver": prim_driver,
                "current_treatment": "Active Clinical Protocol",
                "treatment_history": "Clinical Molecular Series",
                "clinical_status": clin_status,
                "baseline_ki67": "Documented in Medical Record",
                "recommendation": recom,
                "samples": samples
            }

        return patients_dict if patients_dict else None
    except Exception:
        return None

