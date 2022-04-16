from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
import pandas as pd
import numpy as np

df = pd.read_csv("new_smiles.csv")

fps_morgan = []
fps_MACCS  = []
for smile in df["SMILES"]:
    m = Chem.MolFromSmiles(smile)
    fp_morgan = AllChem.GetMorganFingerprint(m, 4)
    fp_MACCS  = MACCSkeys.GenMACCSKeys(m)
    fps_morgan.append(fp_morgan)
    fps_MACCS.append(fp_MACCS)

n = len(fps_morgan)
K_morgan = np.eye(n)
K_maccs = np.eye(n)
for i in range(n):
    for j in range(n):
        K_morgan[i, j] = DataStructs.TanimotoSimilarity(fps_morgan[i], fps_morgan[j])
        K_maccs[i, j]  = DataStructs.TanimotoSimilarity(fps_MACCS[i],fps_MACCS[j])

np.save('MFTSK', K_morgan)
np.save('MACCSTSK', K_maccs)