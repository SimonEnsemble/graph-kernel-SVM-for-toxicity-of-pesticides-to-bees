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

def get_mol(id_mol):
    # get mol
    this_smiles = smiles[id_mol]
    mol = copy.deepcopy(Chem.MolFromSmiles(this_smiles))
    return mol

def save_to_svg(svg_string, svg_filename):
    f = open(svg_filename, "w")
    f.write(svg_string)
    f.close()

# build highlights in molecule by getting these substructure matches
def get_highlight(mol, substructure):
    highlight = mol.GetSubstructMatches(substructure)
    hit_ats = []
    for tp in highlight:
        for i in tp:
            hit_ats.append(i)
    return list(hit_ats)

def viz_maccs(id_mol):
    mol = get_mol(id_mol)

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
        substructure = Chem.MolFromSmarts(smart)

        # get highlights
        highlight = get_highlight(mol, substructure)

        # append these highlights to list of highlights for each molecule
        highlight_list.append(highlight)
    
    svg_string = Draw.MolsToGridImage(mols, molsPerRow=4, highlightAtomLists=highlight_list, subImgSize=(250, 250), useSVG=True)
    save_to_svg(svg_string, "mol_w_highlights_{}.svg".format(id_mol))
    # save to SVG

def viz_one_maccs(id_mol, id_maccs):
    # get mol
    mol = get_mol(id_mol)

    # compute maccs fingerprint; assert this is on
    fp = MACCSkeys.GenMACCSKeys(mol)
#    assert id_maccs in list(fp.GetOnBits())



id_mol = 347
viz_maccs(id_mol)
