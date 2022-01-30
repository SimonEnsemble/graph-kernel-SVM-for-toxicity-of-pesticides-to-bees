push!(LOAD_PATH, joinpath(pwd(), "src"))
using RWK, LinearAlgebra, JLD2, ProgressMeter, MolecularGraph, CSV, DataFrames

# load BeeToxAI dataset
data = CSV.read("new_smiles.csv", DataFrame)

mols = [smilestomol(smiles) for smiles in data[:, "SMILES"]]

addhydrogens_flag = false
if addhydrogens_flag
	mols = addhydrogens.(mols)
end

n_mol = length(mols)

n_jobs = Int(n_mol * (n_mol - 1) / 2 + n_mol)

pbar = Progress(n_jobs, 1)

K = zeros(n_mol, n_mol) # Gram matrix
γ = 0.01
for m = 1:n_mol
	for n = m:n_mol
		dpg = direct_product_graph(mols[m], mols[n], verbose=false)

		K[m, n] = grw_kernel(dpg, γ)
		K[n, m] = K[m, n]

		next!(pbar)
	end
end

# center the Gram matrix
C = I(n_mol) - 1 / n_mol * ones(n_mol, n_mol)
Kcenter = C * K * C

if addhydrogens_flag
	jldsave("BeeToxKwh.jld2"; Kcenter)
else
	jldsave("BeeToxK.jld2"; Kcenter)
end
