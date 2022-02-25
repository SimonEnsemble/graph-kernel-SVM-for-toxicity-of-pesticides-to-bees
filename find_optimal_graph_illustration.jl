push!(LOAD_PATH, joinpath(pwd(), "src"))
using RWK, LinearAlgebra, JLD2, ProgressMeter, MolecularGraph, CSV, DataFrames, Graphs

data = CSV.read("new_smiles.csv", DataFrame)

mols = [smilestomol(smiles) for smiles in data[:, "SMILES"]]

n_mol = length(mols)

n_jobs = Int(n_mol * (n_mol - 1) / 2 + n_mol)

pbar = Progress(n_jobs, 1)

for m = 1:n_mol
	for n = m:n_mol

		dpg = direct_product_graph(mols[m], mols[n], verbose=false)

        if 3<ne(dpg)<10 && nv(dpg)<=24
            println("$m \t $n")
        end

		next!(pbar)
	end
end
