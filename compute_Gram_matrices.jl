using Distributed, SharedArrays, JLD2, MolecularGraph, CSV, DataFrames, ProgressMeter
@everywhere push!(LOAD_PATH, joinpath(pwd(), "src"))
@everywhere using RWK

#=
settings
=#
# kernel = fixed_length_rw_kernel
# params = [1, 2, 3, 4, 5, 6, 7, 8] # l's
kernel = grw_kernel
params = [0.06, 0.05, 0.04, 0.03, 0.02, 0.01] # γ's

println("settings:\n\tkernel = ", kernel)

# set up directory for storing data
savedir = "gram_matrices"
if ! isdir(savedir)
   mkdir(savedir)
end

#=
load BeeToxAI dataset
convert SMILES to graphs
=#
data = CSV.read("new_smiles.csv", DataFrame)
toxicity = data[:, "Outcome"] # targets

mols = [smilestomol(smiles) for smiles in data[:, "SMILES"]] # no hydrogen's included.
#if include_hydrogens
#	mols = addhydrogens.(mols)
#end

#=
 compute Gram matrix for the different kernel params
=#
n_mol = length(mols) # number of molecules
n_jobs = Int(n_mol * (n_mol - 1) / 2 + n_mol) # for progress bar
pbar = Progress(n_jobs, 1)
println("\t# molecules = ", n_mol)
for p in params
    println("\t\tcomputing Gram matrix for kernel param = ", p)

    # set up kernel
    rwk = dpg -> kernel(dpg, p)

    # compute gram matrix
    K = SharedMatrix{Float64}(n_mol, n_mol) # Gram matrix
    @sync @distributed for m = 1:n_mol
        for n = m:n_mol
            dpg = direct_product_graph(mols[m], mols[n], verbose=false)

            K[m, n] = rwk(dpg)
            K[n, m] = K[m, n]
            
            #next!(pbar)
        end
    end
    K = Matrix(K)

    # center the gram matrix
    Kcentered = centered_Gram_matrix(K)
    
    # save to file
    savename = "BeeTox_$(kernel)_$p.jld2"
    jldsave(joinpath(savedir, savename); K, Kcentered, mols, toxicity)
end
