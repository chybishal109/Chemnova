import io
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Table, TableStyle, Spacer, Image
from reportlab.lib import colors
import chem


def pdf(name, smiles, m, d):
    buf = io.BytesIO()
    st = getSampleStyleSheet()
    img = io.BytesIO(); chem.image(m).save(img, format="PNG"); img.seek(0)
    rows = [["Property", "Value"], ["Formula", d["formula"]], ["Mol. weight", d["mw"]], ["LogP", d["logp"]],
            ["TPSA", d["tpsa"]], ["H-bond donors", d["hbd"]], ["H-bond acceptors", d["hba"]],
            ["Rotatable bonds", d["rotb"]], ["Rings", d["rings"]], ["Lipinski violations", d["violations"]]]
    t = Table(rows, colWidths=[200, 200])
    t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#16271e")),
                           ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                           ("GRID", (0, 0), (-1, -1), .4, colors.grey)]))
    SimpleDocTemplate(buf, pagesize=A4).build([
        Paragraph(f"MoleculeVerse 2.0 - {name}", st["Title"]), Paragraph(f"SMILES: {smiles}", st["Normal"]),
        Spacer(1, 12), Image(img, 280, 200), Spacer(1, 12), t])
    return buf.getvalue()
