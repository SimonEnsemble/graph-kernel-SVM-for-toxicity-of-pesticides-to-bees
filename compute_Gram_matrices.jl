push!(LOAD_PATH, joinpath(pwd(), "src"))
using RWK, JLD2, MolecularGraph, CSV, DataFrames

#=
settings
=#
include_hydrogens = false
γs = [0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001] # for geometric rwk
γs = [0.01]

#=
load BeeToxAI dataset
convert SMILES to graphs
=#
data = CSV.read("new_smiles.csv", DataFrame)
toxicity = data[:, "Outcome"] # targets

mols = [smilestomol(smiles) for smiles in data[:, "SMILES"]]
if include_hydrogens
	mols = addhydrogens.(mols)
end

#=
compute Gram matrix for geometric rwk
=#
for γ in γs
    rwk = dpg -> grw_kernel(dpg, γ)
    K = compute_Gram_matrix(mols, rwk) # uncentered
    Kcentered = centered_Gram_matrix(K)

    if include_hydrogens
        savename = "with_hydrogens_file/"
    else
        savename = "none_hydrogens_file/"
    end
    savename *= "BeeTox_Gram_matrix_γ_$γ"
    if include_hydrogens
        savename *= "w_Hs"
    end
    savename *= ".jld2"

	jldsave(savename; K, Kcentered, mols, toxicity)
end

