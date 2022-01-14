using RWK, LinearAlgebra, JLD2, ProgressMeter, UnicodePlots, MolecularGraph, CSV, DataFrames

# load BeeToxAI dataset

df_contact = CSV.read("new_smiles.csv", DataFrame)

df_contact[:, "SMILES"]

contact_tox_mol = [smilestomol(smiles) for smiles in df_contact[:, "SMILES"]]

beetox_mol = atomsymbol.(contact_tox_mol)