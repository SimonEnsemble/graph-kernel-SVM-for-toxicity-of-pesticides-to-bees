using RWK, LinearAlgebra, JLD2, ProgressMeter, MolecularGraph, CSV, DataFrames

# load BeeToxAI dataset

df_contact = CSV.read("new_smiles.csv", DataFrame)

contact_tox_mol = [smilestomol(smiles) for smiles in df_contact[:, "SMILES"]]

addhydrogens_flag = true

if addhydrogens_flag
	contact_tox_mol = addhydrogens.(contact_tox_mol)
end

n_mol = length(contact_tox_mol)

n_jobs = Int(n_mol*(n_mol-1)/2 + n_mol)

pbar = Progress(n_jobs, 1)

K = zeros(n_mol, n_mol)
γ = 0.01
for m = 1:n_mol
	for n = m:n_mol
		# compute grw_kernel between every two different cages
		dpg = direct_product_graph(contact_tox_mol[m], contact_tox_mol[n], verbose = false)
		K[m, n] = grw_kernel(dpg, γ)
		K[n, m] = K[m, n]
		next!(pbar)
	end
end

N = size(K)[1]

C = I(N) - 1/N * ones(N, N)

Kcenter = C*K*C

if addhydrogens_flag
	jldsave("BeeToxKwh.jld2"; Kcenter)
else
	jldsave("BeeToxK.jld2"; Kcenter)
end