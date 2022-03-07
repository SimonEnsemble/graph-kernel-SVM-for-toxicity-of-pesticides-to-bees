@everywhere using Distributed, SharedArrays, MolecularGraph, CSV, DataFrames, BenchmarkTools

@everywhere push!(LOAD_PATH, joinpath(pwd(), "src"))
@everywhere using RWK

data = CSV.read("new_smiles.csv", DataFrame)

mols = [smilestomol(smiles) for smiles in data[:, "SMILES"]]

n_mol = length(mols) # number of molecules

@everywhere function parallel_compute_gram_matrix(kernel::Function, p::Int64, K::SharedMatrix{Float64})
    rwk = dpg -> kernel(dpg, p)

    @sync @distributed for m = 1:n_mol
        for n = m:n_mol
            dpg = direct_product_graph(mols[m], mols[n], verbose = false)

            K[m, n] = rwk(dpg)
            K[n, m] = K[m, n]

        end
    end
    return K
end

kernel = fixed_length_rw_kernel
p = 12
K = SharedMatrix{Float64}(n_mol, n_mol)

@btime parallel_compute_gram_matrix(kernel, p, K)