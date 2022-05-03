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
K_TS = np.eye(n) # TanimotoSimilarity
K_DP = np.eye(n) # dot product
for i in range(n):
    x_i = np.array(fps[i].ToList())
    for j in range(n):
        x_j = np.array(fps[j].ToList())
        K_DP[i, j] = np.dot(x_i, x_j)
        K_TS[i, j] = DataStructs.TanimotoSimilarity(fps[i],fps[j])

np.save('MACCS_TS_matrix', K_TS)
np.save('MACCS_DP_matrix', K_DP)
