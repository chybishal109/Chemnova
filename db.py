import sqlite3
import chem

DB = "moleculeverse.db"

SEED = [
    ("Vanillin", "Flavour", "COc1cc(C=O)ccc1O", "Main component of vanilla flavour."),
    ("Menthol", "Flavour", "CC(C)C1CCC(C)CC1O", "Cooling terpene alcohol from mint."),
    ("Capsaicin", "Natural product", "COc1cc(CNC(=O)CCCCC=CC(C)C)ccc1O", "Pungent compound of chilli peppers."),
    ("Limonene", "Terpene", "CC1=CCC(CC1)C(C)=C", "Citrus-scented cyclic terpene."),
    ("Citric acid", "Acid", "OC(=O)CC(O)(CC(O)=O)C(O)=O", "Tricarboxylic acid of the Krebs cycle."),
    ("Dopamine", "Neurotransmitter", "NCCc1ccc(O)c(O)c1", "Reward and movement neurotransmitter."),
    ("Serotonin", "Neurotransmitter", "NCCc1c[nH]c2ccc(O)cc12", "Mood and sleep regulating amine."),
    ("Melatonin", "Hormone", "COc1ccc2[nH]cc(CCNC(C)=O)c2c1", "Sleep-cycle hormone."),
    ("Ascorbic acid", "Vitamin", "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O", "Vitamin C, antioxidant."),
    ("Naproxen", "Drug", "COc1ccc2cc(ccc2c1)[C@H](C)C(O)=O", "NSAID pain reliever."),
    ("Lidocaine", "Drug", "CCN(CC)CC(=O)Nc1c(C)cccc1C", "Local anaesthetic."),
    ("Resveratrol", "Natural product", "Oc1ccc(C=Cc2cc(O)cc(O)c2)cc1", "Polyphenol found in grape skins."),
]


def conn():
    c = sqlite3.connect(DB)
    c.row_factory = sqlite3.Row
    return c


def init():
    with conn() as c:
        c.execute("""CREATE TABLE IF NOT EXISTS molecules(
            id INTEGER PRIMARY KEY, name TEXT, category TEXT, smiles TEXT UNIQUE, notes TEXT,
            formula TEXT, mw REAL, logp REAL, tpsa REAL, hbd INT, hba INT, rotb INT, rings INT, violations INT)""")
        c.execute("CREATE TABLE IF NOT EXISTS scores(id INTEGER PRIMARY KEY, player TEXT, score INT, total INT, ts DEFAULT CURRENT_TIMESTAMP)")
        if c.execute("SELECT COUNT(*) FROM molecules").fetchone()[0] == 0:
            for n, cat, s, note in SEED:
                add(n, cat, s, note, c)


def add(name, category, smiles, notes="", c=None):
    m = chem.parse(smiles)
    if m is None:
        raise ValueError("Invalid SMILES")
    d = chem.analyze(m)
    own = c is None
    c = c or conn()
    c.execute("INSERT OR IGNORE INTO molecules(name,category,smiles,notes,formula,mw,logp,tpsa,hbd,hba,rotb,rings,violations) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)",
              (name, category, smiles, notes, d["formula"], d["mw"], d["logp"], d["tpsa"], d["hbd"], d["hba"], d["rotb"], d["rings"], d["violations"]))
    if own:
        c.commit(); c.close()


def all_rows():
    with conn() as c:
        return [dict(r) for r in c.execute("SELECT * FROM molecules ORDER BY name")]


def delete(mid):
    with conn() as c:
        c.execute("DELETE FROM molecules WHERE id=?", (mid,))


def save_score(player, score, total):
    with conn() as c:
        c.execute("INSERT INTO scores(player,score,total) VALUES(?,?,?)", (player or "Guest", score, total))


def scores():
    with conn() as c:
        return [dict(r) for r in c.execute("SELECT player,score,total,ts FROM scores ORDER BY score DESC, ts DESC LIMIT 10")]
