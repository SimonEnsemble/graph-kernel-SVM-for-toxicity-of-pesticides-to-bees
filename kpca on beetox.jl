using RWK, LinearAlgebra, JLD2, ProgressMeter, UnicodePlots, MolecularGraph, CSV, DataFrames

# load BeeToxAI dataset

df_contact = CSV.read("BeeToxAI Data/File S1 Acute contact toxicity dataset for classification.csv", DataFrame)

df_contact[!, "SMILES"]

contact_tox_mol = [smilestomol(smiles) for smiles in df_contact[30:383, "SMILES"]]

smilestomol("COc1nn(CSP(=S)(OC)OC)c(=O)s1")

df_contact[29, "SMILES"]