push!(LOAD_PATH, joinpath(pwd(), "src"))
using Distributed, SharedArrays, JLD2, MolecularGraph, CSV, DataFrames, ProgressMeter
@everywhere using RWK

#=
settings
=#
include_hydrogens = false
kernels = [fixed_length_rw_kernel, grw_kernel]
params = Dict(
     # grw_kernel            => [0.05, 0.04, 0.03, 0.02, 0.01], # γ's
     grw_kernel            => [],
     # fixed_length_rw_kernel => [1, 2, 3, 4, 5] # l's
     fixed_length_rw_kernel => [1] # l's
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

@everywhere mols = [smilestomol(smiles) for smiles in data[1:10, "SMILES"]]
if include_hydrogens
	mols = addhydrogens.(mols)
end

# set up kernel
@everywhere rwk = dpg -> fixed_length_rw_kernel(dpg, 2)

# compute gram matrix
@everywhere n_mol = length(mols) # number of molecules
pbar = Progress(n_mol, 1)
K = SharedMatrix{Float64}(n_mol, n_mol) # Gram matrix
@sync @distributed for m = 1:n_mol
    for n = m:n_mol
        dpg = direct_product_graph(mols[m], mols[n], verbose=false)

        K[m, n] = rwk(dpg)
        K[n, m] = K[m, n]
    end
end
println(Matrix(K))
#=
compute Gram matrix
=#
#for kernel in kernels
#    println("kernel = ", kernel)
#    for p in params[kernel]
#        println("\tparam = ", p)
#        # set up kernel
#        rwk = dpg -> kernel(dpg, p)
#        
#        # compute gram matrix
#        @everywhere n_mol = length(mols) # number of molecules
#        pbar = Progress(n_mol, 1)
#        K = SharedMatrix{Float64}(n_mol, n_mol) # Gram matrix
#        @distributed for m = 1:n_mol
#            #next!(pbar)
#            for n = m:n_mol
#                dpg = direct_product_graph(mols[m], mols[n], verbose=false)
#
#                K[m, n] = rwk(dpg)
#                K[n, m] = K[m, n]
#            end
#        end
#
#        K = Matrix(K)
#        
#        # center the gram matrix
#        Kcentered = centered_Gram_matrix(K)
#        
#        # save to file
#        savename = "BeeTox_$(kernel)_$p.jld2"
#        jldsave(joinpath(savedir, savename); K, Kcentered, mols, toxicity)
#    end
#end
