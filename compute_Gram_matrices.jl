push!(LOAD_PATH, joinpath(pwd(), "src"))
using RWK, JLD2, MolecularGraph, CSV, DataFrames

#=
settings
=#
include_hydrogens = false
kernels = [fixed_length_rw_kernel, grw_kernel]
params = Dict(
     #"geometric" => [0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001], # for geometric rwk
     grw_kernel            => [0.05, 0.04, 0.03, 0.02, 0.01], # γ's
     fixed_length_rw_kernel => [1, 2, 3, 4, 5] # l's
)

# set up directory for storing data
savedir = "include_hydrogens_$include_hydrogens"
if ! isdir(savedir)
   mkdir(savedir)
end

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
compute Gram matrix
=#
for kernel in kernels
    println("kernel = ", kernel)
    for p in params[kernel]
        println("\tparam = ", p)
        # set up kernel
        rwk = dpg -> kernel(dpg, p)
        
        # compute gram matrix
        K = compute_Gram_matrix(mols, rwk) # uncentered
        
        # center the gram matrix
        Kcentered = centered_Gram_matrix(K)
        
        # save to file
        savename = "BeeTox_$(kernel)_$p.jld2"
        jldsave(joinpath(savedir, savename); K, Kcentered, mols, toxicity)
    end
end
