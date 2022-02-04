push!(LOAD_PATH, joinpath(pwd(), "src"))
using RWK, JLD2, MolecularGraph, CSV, DataFrames

# settings
addhydrogens_flag = false
γs = [0.01]

# load BeeToxAI dataset
data = CSV.read("new_smiles.csv", DataFrame)
toxicity = data[:, "Outcome"]

# convert SMILES to graphs
mols = [smilestomol(smiles) for smiles in data[:, "SMILES"]]
if addhydrogens_flag
	mols = addhydrogens.(mols)
end

for γ in γs
    K = compute_Gram_matrix(mols, γ) # uncentered
    Kcentered = centered_Gram_matrix(K)

    savename = "BeeTox_Gram_matrix_γ_$γ"
    if addhydrogens_flag
        savename *= "w_Hs"
    end
    savename *= ".jld2"

	jldsave(savename; K, Kcentered, mols, toxicity)
end
