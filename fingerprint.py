from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
import pandas as pd
import numpy as np

# load bee tox data
df = pd.read_csv("new_smiles.csv")

# compute and store fingerprints here.
fps  = []
for smile in df["SMILES"]:
    m = Chem.MolFromSmiles(smile)
    fps.append(
        MACCSkeys.GenMACCSKeys(m)
    )

# compute pairwise similarity matrix
n = len(fps)
K = np.eye(n)
for i in range(n):
    for j in range(n):
        K[i, j]  = DataStructs.TanimotoSimilarity(fps[i],fps[j])

np.save('MACCS_TS_matrix', K)
