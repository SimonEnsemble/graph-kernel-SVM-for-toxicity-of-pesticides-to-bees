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

hl_color = (0.5, 0.97, 0.73, 0.0) # first one is alpha
grid_image_kwargs = {"molsPerRow": 12, "useSVG": True, "subImgSize": (250, 250)}

def build_highlight_colors():
    atom_colors = {}
    bond_colors = {}
    for i in range(100):
        atom_colors[i] = hl_color
        bond_colors[i] = hl_color
    return atom_colors, bond_colors

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
def get_highlight(mol, substructure, just_one):
    if just_one:
        highlight = mol.GetSubstructMatch(substructure)
        hit_ats = []
        for tp in highlight:
            hit_ats.append(tp)
    else:
        highlight = mol.GetSubstructMatches(substructure)
        hit_ats = []
        for tp in highlight:
            for i in tp:
                hit_ats.append(i)
    return list(hit_ats)

def viz_molecule(id_mol):
    mol = get_mol(id_mol)
    svg_string = Draw.MolsToGridImage([mol], **grid_image_kwargs)
    #molsPerRow=4, subImgSize=(250, 250), useSVG=True)
    save_to_svg(svg_string, "mol_{}.svg".format(id_mol))

def viz_maccs(id_mol, just_one=False):
    mol = get_mol(id_mol)

    # compute maccs fingerprint
    fp = MACCSkeys.GenMACCSKeys(mol)

    # find which bits are on
    onbits = list(fp.GetOnBits())
    
    # create list of the molecules for different highlights
    mols = [mol] * len(onbits)
    highlight_list = []
    atom_colors = []
    bond_colors = []
    for onbit in onbits:
        # get smarts pattern and substructure corresponding to this on bit
        smart = smartsPatts[onbit][0]
        substructure = Chem.MolFromSmarts(smart)
        print("\ton bit: ", onbit)
        print("\tsmart: ", smart)

        # get highlights
        highlight = get_highlight(mol, substructure, just_one)

        # append these highlights to list of highlights for each molecule
        highlight_list.append(highlight)
    
        atom_color, bond_color = build_highlight_colors()
        atom_colors.append(atom_color)
        bond_colors.append(bond_color)
    svg_string = Draw.MolsToGridImage(mols, highlightAtomLists=highlight_list, **grid_image_kwargs, 
                                      highlightAtomColors=atom_colors, highlightBondColors=bond_colors)
    svg_filename = "mol_w_highlights_{}_{}.svg".format(id_mol, just_one)
    save_to_svg(svg_string, svg_filename)
    print("see ", svg_filename)
    # save to SVG

def viz_one_maccs(id_mol, onbit, just_one=False):
    # get mol
    mol = get_mol(id_mol)

    # compute maccs fingerprint; assert this is on
    fp = MACCSkeys.GenMACCSKeys(mol)
    assert onbit in list(fp.GetOnBits())

    smart = smartsPatts[onbit][0]
    
    # get smarts pattern and substructure corresponding to this on bit
    smart = smartsPatts[onbit][0]
    substructure = Chem.MolFromSmarts(smart)

    # get highlights
    highlight = get_highlight(mol, substructure, just_one)
    
    svg_string = Draw.MolsToGridImage([mol], highlightAtomLists=[highlight], **grid_image_kwargs)
    svg_filename = "mol_w_highlights_{}_onbit_{}_{}.svg".format(id_mol, onbit, just_one)
    save_to_svg(svg_string, svg_filename)
    print("see ", svg_filename)

# rsvg-convert -f pdf -o mygraph.pdf mol_w_highlights_347.svg to convert svg to pdf

id_mol = 347
viz_maccs(id_mol, False)
viz_maccs(id_mol, True)
viz_molecule(id_mol)
