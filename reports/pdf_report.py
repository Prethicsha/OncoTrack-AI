"""
OncoTrack Clinical PDF Summary Report Generator
Generates a downloadable clinical molecular pathology summary report using fpdf2.
Includes strict Latin-1 / ASCII sanitization to prevent Unicode encoding issues.
"""

import io
from typing import Dict, Any
from fpdf import FPDF


def sanitize_text(text: Any) -> str:
    """Sanitizes unicode characters (em-dashes, Greek symbols, curly quotes) for standard PDF fonts."""
    if text is None:
        return ""
    s = str(text)
    replacements = {
        "—": "-",
        "–": "-",
        "’": "'",
        "‘": "'",
        "“": '"',
        "”": '"',
        "≥": ">=",
        "≤": "<=",
        "Δ": "Delta ",
        "•": "*",
        "…": "...",
        "µ": "u",
        "²": "^2",
        "³": "^3",
        "±": "+/-",
        "~": "~",
        "→": "->",
        "←": "<-",
        "↑": "^",
        "↓": "v",
        "🚨": "[ALERT] ",
        "⚠️": "[WARNING] ",
        "✅": "[STABLE] ",
        "⚪": "[NORMAL] ",
        "🔴": "[ALERT] ",
        "💡": "[NOTE] ",
        "📋": "",
        "🧬": "",
        "🏥": "",
        "🎥": "",
        "⚙️": "",
        "📤": "",
        "📄": "",
        "📥": ""
    }
    for orig, rep in replacements.items():
        s = s.replace(orig, rep)
    
    # Fallback encode to latin-1
    return s.encode("latin-1", "replace").decode("latin-1")


class ClinicalReportPDF(FPDF):
    def header(self):
        self.set_fill_color(11, 17, 32)
        self.rect(0, 0, 210, 24, "F")
        self.set_text_color(248, 250, 252)
        self.set_font("Helvetica", "B", 14)
        self.set_xy(10, 6)
        self.cell(0, 8, sanitize_text("ONCOTRACK AI - CLINICAL MOLECULAR DYNAMICS REPORT"), ln=True)
        self.set_font("Helvetica", "", 8)
        self.set_text_color(148, 163, 184)
        self.set_xy(10, 14)
        self.cell(0, 5, sanitize_text("Longitudinal ctDNA Liquid Biopsy & Clonal Evolution Analysis | Research Prototype"), ln=True)
        self.ln(8)

    def footer(self):
        self.set_y(-15)
        self.set_font("Helvetica", "I", 7)
        self.set_text_color(148, 163, 184)
        self.cell(0, 10, sanitize_text("Page " + str(self.page_no()) + " | OncoTrack AI Research Prototype | Not a certified primary diagnostic device."), align="C")


def generate_clinical_pdf(patient_data: Dict[str, Any], score_data: Dict[str, Any]) -> bytes:
    """Generates PDF report for patient chart download."""
    pdf = ClinicalReportPDF(orientation="P", unit="mm", format="A4")
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()

    # 1. Patient & Regimen Demographics Table
    pdf.set_font("Helvetica", "B", 11)
    pdf.set_text_color(30, 41, 59)
    pdf.cell(0, 7, sanitize_text("1. PATIENT DEMOGRAPHICS & CLINICAL REGIMEN"), ln=True)
    pdf.set_draw_color(203, 213, 225)
    pdf.line(10, pdf.get_y(), 200, pdf.get_y())
    pdf.ln(3)

    pdf.set_font("Helvetica", "", 9)
    pdf.set_text_color(51, 65, 85)
    
    col1 = 45
    col2 = 145

    fields = [
        ("Patient Identifier:", patient_data.get("patient_id", "N/A")),
        ("Primary Malignancy:", patient_data.get("cancer_type", "N/A")),
        ("Clinical Stage:", patient_data.get("stage", "N/A")),
        ("Current Treatment Regimen:", patient_data.get("current_treatment", "N/A")),
        ("Treatment History:", patient_data.get("treatment_history", "N/A")),
        ("Baseline Tissue Ki-67:", patient_data.get("baseline_ki67", "N/A")),
        ("Primary Driver Mutation:", patient_data.get("primary_driver", "N/A")),
    ]

    for label, val in fields:
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(col1, 6, sanitize_text(label), border=0)
        pdf.set_font("Helvetica", "", 9)
        pdf.cell(col2, 6, sanitize_text(str(val)), border=0, ln=True)

    pdf.ln(4)

    # 2. Molecular Status & Escape Score
    pdf.set_font("Helvetica", "B", 11)
    pdf.set_text_color(30, 41, 59)
    pdf.cell(0, 7, sanitize_text("2. LONGITUDINAL MOLECULAR STATUS & CLONAL ESCAPE SCORE"), ln=True)
    pdf.line(10, pdf.get_y(), 200, pdf.get_y())
    pdf.ln(3)

    status = score_data.get("status", "MONITOR")
    score = score_data.get("score", 0)

    # Highlight box for score
    if score >= 70:
        pdf.set_fill_color(254, 226, 226)
        pdf.set_text_color(185, 28, 28)
    elif score >= 40:
        pdf.set_fill_color(254, 243, 199)
        pdf.set_text_color(180, 83, 9)
    else:
        pdf.set_fill_color(209, 250, 229)
        pdf.set_text_color(4, 120, 87)

    pdf.set_font("Helvetica", "B", 12)
    pdf.cell(0, 10, sanitize_text(f"STATUS: {status}  |  CLONAL ESCAPE SCORE: {score}/100"), border=1, ln=True, fill=True, align="C")
    pdf.ln(2)

    pdf.set_font("Helvetica", "", 9)
    pdf.set_text_color(51, 65, 85)
    pdf.multi_cell(0, 5, sanitize_text(f"Clinical Interpretation: {score_data.get('summary', '')}"))
    pdf.ln(2)
    pdf.set_font("Helvetica", "I", 8.5)
    pdf.multi_cell(0, 4.5, sanitize_text(f"Guidance: {patient_data.get('recommendation', '')}"))
    pdf.ln(4)

    # 3. Deterministic Clinical Audit Factors
    pdf.set_font("Helvetica", "B", 11)
    pdf.set_text_color(30, 41, 59)
    pdf.cell(0, 7, sanitize_text("3. DETERMINISTIC CLINICAL AUDIT FACTORS (XAI)"), ln=True)
    pdf.line(10, pdf.get_y(), 200, pdf.get_y())
    pdf.ln(3)

    for idx, f in enumerate(score_data.get("factors", []), start=1):
        is_pos = f.get("is_positive", False)
        prefix = "[ALERT]" if is_pos else "[STABLE]"
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(0, 5, sanitize_text(f"{idx}. {f.get('title')} ({prefix})"), ln=True)
        pdf.set_font("Helvetica", "", 8.5)
        pdf.multi_cell(0, 4.5, sanitize_text(f"   - {f.get('description')}"))
        pdf.ln(1)

    pdf.ln(4)

    # 4. Serial ctDNA Observation Log Table
    pdf.set_font("Helvetica", "B", 11)
    pdf.set_text_color(30, 41, 59)
    pdf.cell(0, 7, sanitize_text("4. SERIAL ctDNA OBSERVATION LOG"), ln=True)
    pdf.line(10, pdf.get_y(), 200, pdf.get_y())
    pdf.ln(3)

    pdf.set_font("Helvetica", "B", 8)
    pdf.set_fill_color(241, 245, 249)
    pdf.cell(20, 6, sanitize_text("Cycle"), border=1, fill=True)
    pdf.cell(30, 6, sanitize_text("Date"), border=1, fill=True)
    pdf.cell(25, 6, sanitize_text("Gene"), border=1, fill=True)
    pdf.cell(35, 6, sanitize_text("Mutation"), border=1, fill=True)
    pdf.cell(30, 6, sanitize_text("Lineage"), border=1, fill=True)
    pdf.cell(25, 6, sanitize_text("VAF (%)"), border=1, fill=True)
    pdf.cell(25, 6, sanitize_text("Depth"), border=1, fill=True, ln=True)

    pdf.set_font("Helvetica", "", 8)
    for s in patient_data.get("samples", []):
        if not s.get("is_forecast", False):
            pdf.cell(20, 5, sanitize_text(f"C{s['cycle']}"), border=1)
            pdf.cell(30, 5, sanitize_text(str(s.get('date', f"C{s['cycle']}"))), border=1)
            pdf.cell(25, 5, sanitize_text(str(s['gene'])), border=1)
            pdf.cell(35, 5, sanitize_text(str(s['mutation'])), border=1)
            pdf.cell(30, 5, sanitize_text(str(s['clone_type']).capitalize()), border=1)
            pdf.cell(25, 5, sanitize_text(f"{s['vaf']:.2f}%"), border=1)
            pdf.cell(25, 5, sanitize_text(f"{s.get('total_depth', 4000):,}x"), border=1, ln=True)

    pdf.ln(6)
    pdf.set_font("Helvetica", "I", 7.5)
    pdf.set_text_color(100, 116, 139)
    pdf.multi_cell(0, 3.5, sanitize_text("SCIENTIFIC DISCLAIMER: OncoTrack AI is a research and clinical hackathon prototype using retrospective/synthetic longitudinal datasets. All trajectories and kinetic indices are non-diagnostic decision support aids and must be validated in prospective clinical trials."))

    return bytes(pdf.output())
