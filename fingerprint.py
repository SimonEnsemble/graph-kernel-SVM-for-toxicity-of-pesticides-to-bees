from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np

df = pd.read_csv("new_smiles.csv")

fps = []
for smile in df["SMILES"]:
    m = Chem.MolFromSmiles(smile)
    fp = AllChem.GetMorganFingerprint(m, 4)
    fps.append(fp)

n = len(fps)
K = np.eye(n)
for i in range(n):
    for j in range(n):
        K[i, j] = DataStructs.DiceSimilarity(fps[i], fps[j])
        K[j, i] = K[i, j]


np.save('MFDSK', K)