from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, AllChem, Draw


def parse(smiles):
    return Chem.MolFromSmiles((smiles or "").strip()) if (smiles or "").strip() else None


def analyze(m):
    d = {
        "formula": rdMolDescriptors.CalcMolFormula(m),
        "mw": round(Descriptors.MolWt(m), 2),
        "logp": round(Crippen.MolLogP(m), 2),
        "tpsa": round(rdMolDescriptors.CalcTPSA(m), 1),
        "hbd": rdMolDescriptors.CalcNumHBD(m),
        "hba": rdMolDescriptors.CalcNumHBA(m),
        "rotb": rdMolDescriptors.CalcNumRotatableBonds(m),
        "rings": rdMolDescriptors.CalcNumRings(m),
    }
    d["violations"] = int(d["mw"] > 500) + int(d["logp"] > 5) + int(d["hbd"] > 5) + int(d["hba"] > 10)
    return d


def image(m, size=(380, 280)):
    return Draw.MolToImage(m, size=size)


def mol_block(m):
    h = Chem.AddHs(m)
    if AllChem.EmbedMolecule(h, randomSeed=42) != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(h)
    except Exception:
        pass
    return Chem.MolToMolBlock(h)
