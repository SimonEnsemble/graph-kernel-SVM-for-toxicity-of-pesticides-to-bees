from rdkit import Chem
from rdkit.Chem import MACCSkeys
import numpy as np
import pandas as pd
from rdkit.Chem import Draw
from rdkit.Chem.MACCSkeys import smartsPatts
from rdkit.Chem.Draw import IPythonConsole
import copy

df = pd.read_csv("new_smiles.csv")
mols = [Chem.MolFromSmiles(smile) for smile in df["SMILES"]]
fps = [MACCSkeys.GenMACCSKeys(mol) for mol in mols]

id_mol = 347

def viz_maccs(id_mol):
    onbits = list(fps[id_mol].GetOnBits())
    m = copy.deepcopy(mols[id_mol])
    mMols = [m] * len(onbits)
    highlight_list = []
    for onbit in onbits:
        smart = smartsPatts[onbit][0]
        sub = Chem.MolFromSmarts(smart)
        highlight = m.GetSubstructMatches(sub)
        hit_ats = []
        for tp in highlight:
            for i in tp:
                hit_ats.append(i)
        highlight_list.append(list(hit_ats))

    img = Draw.MolsToGridImage(mMols, molsPerRow=4, highlightAtomLists=highlight_list, subImgSize=(250, 250), useSVG=False)
    img.save("viz_MACCS.png")

viz_maccs(id_mol)