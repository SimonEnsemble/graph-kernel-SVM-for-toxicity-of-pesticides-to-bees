from rdkit import Chem

from pandas import pandas

df = pandas.read_csv("BeeToxAI Data/File S1 Acute contact toxicity dataset for classification.csv")

new_smiles = []
for smile in df["SMILES"]:
    m = Chem.MolFromSmiles(smile)
    new_smile = Chem.MolToSmiles(m, kekuleSmiles = True)
    new_smiles.append(new_smile)

df["SMILES"] = new_smiles
df.to_csv("new_smiles.csv", encoding = 'utf-8', index = False)