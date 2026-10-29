import json
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
import streamlit.components.v1 as components

import chem, db, report

st.set_page_config(page_title="MoleculeVerse 2.0", page_icon="🧪", layout="wide")
db.init()

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Fraunces:wght@600;800&family=DM+Sans:wght@400;500;700&family=DM+Mono&display=swap');
html, body, [class*="css"], .stMarkdown p { font-family:'DM Sans',sans-serif; letter-spacing:.2px; }
h1,h2,h3 { font-family:'Fraunces',serif !important; letter-spacing:-.5px; color:#f5b335; }
code, .stDataFrame { font-family:'DM Mono',monospace !important; }
.stApp { background:radial-gradient(circle at 15% 0%, #1b3326 0%, #0e1a14 55%); }
[data-testid=stSidebar] { background:#0a130e; border-right:1px solid #24402f; }
.card { background:#16271e; border:1px solid #2b4a38; border-radius:14px; padding:16px 20px; }
.card b { font-size:1.9rem; color:#b6f36b; font-family:'Fraunces',serif; }
.card span { display:block; color:#9db8a7; font-size:.8rem; text-transform:uppercase; letter-spacing:1.5px; }
</style>""", unsafe_allow_html=True)

rows = db.all_rows()
df = pd.DataFrame(rows)
PROPS = {"mw": "Mol. weight", "logp": "LogP", "tpsa": "TPSA", "hbd": "HBD", "hba": "HBA", "rotb": "Rot. bonds", "rings": "Rings"}


def viewer(block, style="stick"):
    sty = {"stick": {"stick": {}}, "sphere": {"sphere": {"scale": .3}}, "ball & stick": {"stick": {"radius": .15}, "sphere": {"scale": .25}}}[style]
    components.html(f"""<script src="https://3dmol.org/build/3Dmol-min.js"></script>
    <div id="v" style="height:380px;width:100%;position:relative;border-radius:12px;overflow:hidden"></div>
    <script>let v=$3Dmol.createViewer("v",{{backgroundColor:"#0a130e"}});
    v.addModel(`{block}`,"sdf");v.setStyle({{}},{json.dumps(sty)});v.zoomTo();v.spin(true);v.render();</script>""", height=390)


def show_mol(r):
    m = chem.parse(r["smiles"])
    c1, c2 = st.columns(2)
    c1.image(chem.image(m), caption="2D structure")
    with c2:
        style = st.radio("3D style", ["ball & stick", "stick", "sphere"], horizontal=True, key=f"s{r['id']}")
        blk = chem.mol_block(m)
        viewer(blk, style) if blk else st.warning("3D embedding failed.")


page = st.sidebar.radio("Navigate", ["Dashboard", "Molecular Library", "Add & Analyze", "Compare", "Assessments", "About"])
st.sidebar.caption("MoleculeVerse 2.0 · RDKit · 3Dmol.js")

if page == "Dashboard":
    st.title("🧪 MoleculeVerse 2.0")
    k = st.columns(4)
    for col, (lab, val) in zip(k, [("Molecules", len(df)), ("Avg MW", round(df.mw.mean(), 1)),
                                    ("Avg LogP", round(df.logp.mean(), 2)), ("Lipinski-pass", int((df.violations == 0).sum()))]):
        col.markdown(f"<div class='card'><span>{lab}</span><b>{val}</b></div>", unsafe_allow_html=True)
    a, b = st.columns(2)
    a.plotly_chart(px.scatter(df, x="mw", y="logp", color="category", hover_name="name", size="tpsa", title="MW vs LogP",
                              template="plotly_dark", color_discrete_sequence=px.colors.qualitative.Bold), use_container_width=True)
    b.plotly_chart(px.histogram(df, x="category", color="category", title="Molecules by category", template="plotly_dark"), use_container_width=True)

elif page == "Molecular Library":
    st.title("Molecular Library")
    q = st.text_input("Search name / formula / SMILES")
    cats = st.multiselect("Category", sorted(df.category.unique()))
    f = df[df.apply(lambda r: q.lower() in f"{r['name']} {r['formula']} {r['smiles']}".lower(), axis=1)]
    if cats:
        f = f[f.category.isin(cats)]
    st.dataframe(f[["name", "category", "formula", "mw", "logp", "tpsa", "violations"]], use_container_width=True, hide_index=True)
    if len(f):
        pick = st.selectbox("Inspect molecule", f.name)
        r = f[f.name == pick].iloc[0].to_dict()
        st.subheader(r["name"]); st.caption(r["notes"]); st.code(r["smiles"])
        show_mol(r)
        st.download_button("⬇ PDF report", report.pdf(r["name"], r["smiles"], chem.parse(r["smiles"]), r), f"{r['name']}.pdf")
        if st.button("🗑 Delete from library"):
            db.delete(int(r["id"])); st.rerun()

elif page == "Add & Analyze":
    st.title("Add & Analyze")
    name = st.text_input("Name"); cat = st.text_input("Category", "Custom"); notes = st.text_area("Notes")
    smi = st.text_input("SMILES", "CC(=O)Oc1ccccc1C(=O)O")
    m = chem.parse(smi)
    if m is None:
        st.error("Invalid SMILES")
    else:
        d = chem.analyze(m)
        c1, c2 = st.columns(2)
        c1.image(chem.image(m))
        c2.json(d)
        st.success("Passes Lipinski's rule of five") if d["violations"] == 0 else st.warning(f"{d['violations']} Lipinski violation(s)")
        blk = chem.mol_block(m)
        if blk: viewer(blk)
        st.download_button("⬇ PDF report", report.pdf(name or "Unnamed", smi, m, d), "report.pdf")
        if st.button("💾 Save to library") and name:
            db.add(name, cat, smi, notes); st.success("Saved"); st.rerun()

elif page == "Compare":
    st.title("Compare")
    sel = st.multiselect("Pick 2–4 molecules", list(df.name), default=list(df.name[:2]), max_selections=4)
    if len(sel) >= 2:
        sub = df[df.name.isin(sel)]
        st.dataframe(sub.set_index("name")[list(PROPS)].T, use_container_width=True)
        fig = go.Figure()
        for _, r in sub.iterrows():
            vals = [r[p] / df[p].max() if df[p].max() else 0 for p in PROPS]
            fig.add_trace(go.Scatterpolar(r=vals + vals[:1], theta=list(PROPS.values()) + [list(PROPS.values())[0]], fill="toself", name=r["name"]))
        st.plotly_chart(fig.update_layout(template="plotly_dark", title="Normalised property radar"), use_container_width=True)
        cols = st.columns(len(sel))
        for col, (_, r) in zip(cols, sub.iterrows()):
            col.image(chem.image(chem.parse(r["smiles"]), (300, 220)), caption=r["name"])

elif page == "Assessments":
    st.title("Assessments")
    qs = [("Lipinski: maximum molecular weight?", ["300", "500", "800"], "500"),
          ("Lipinski: maximum H-bond donors?", ["2", "5", "10"], "5"),
          ("Which has the higher LogP (more lipophilic)?", ["Citric acid", "Limonene"], "Limonene"),
          ("SMILES stands for…", ["Simplified Molecular-Input Line-Entry System", "Standard Molecular Index of Linked Elements"], "Simplified Molecular-Input Line-Entry System"),
          ("TPSA predicts mainly…", ["Membrane permeability", "Boiling point", "Colour"], "Membrane permeability")]
    player = st.text_input("Your name")
    ans = [st.radio(q, o, index=None, key=f"q{i}") for i, (q, o, _) in enumerate(qs)]
    if st.button("Submit"):
        sc = sum(a == c for a, (_, _, c) in zip(ans, qs))
        st.metric("Score", f"{sc}/{len(qs)}"); db.save_score(player, sc, len(qs))
    st.subheader("Leaderboard"); st.dataframe(pd.DataFrame(db.scores()), hide_index=True, use_container_width=True)

else:
    st.title("About")
    st.markdown("**MoleculeVerse 2.0** is a cheminformatics explorer built with Streamlit, RDKit, Plotly, 3Dmol.js, SQLite and ReportLab. "
                "Browse a curated library, analyse any SMILES, compare molecules, export PDF reports and test yourself.")
