from rdkit import Chem
from rdkit.Chem import MACCSkeys
import numpy as np
import pandas as pd
from rdkit.Chem import Draw
from rdkit.Chem.MACCSkeys import smartsPatts
import copy

# load SMILES and compute MACCS fingerprints
df = pd.read_csv("new_smiles.csv")
smiles = [smile for smile in df["SMILES"]]
mols = [Chem.MolFromSmiles(smile) for smile in df["SMILES"]]

def viz_maccs(id_mol):
    # get mol
    this_smiles = smiles[id_mol]
    mol = copy.deepcopy(Chem.MolFromSmiles(this_smiles))

    # compute maccs fingerprint
    fp = MACCSkeys.GenMACCSKeys(mol)

    # find which bits are on
    onbits = list(fp.GetOnBits())
    
    # create list of the molecules for different highlights
    mols = [mol] * len(onbits)
    highlight_list = []
    for onbit in onbits:
        # get smarts pattern and substructure corresponding to this on bit
        smart = smartsPatts[onbit][0]
        sub = Chem.MolFromSmarts(smart)

        # build highlights in molecule by getting these substructure matches
        highlight = mol.GetSubstructMatches(sub)
        hit_ats = []
        for tp in highlight:
            for i in tp:
                hit_ats.append(i)

        # append these highlights to list of highlights for each molecule
        highlight_list.append(list(hit_ats))
    
    img = Draw.MolsToGridImage(mols, molsPerRow=4, highlightAtomLists=highlight_list, subImgSize=(250, 250), useSVG=True)
    # save to SVG
    f = open("mol_w_highlights_{}.svg".format(id_mol), "w")
    f.write(img)
    f.close()

id_mol = 347
viz_maccs(id_mol)
