using RWK, LinearAlgebra, JLD2, ProgressMeter, UnicodePlots, MolecularGraph, CSV, DataFrames

# load BeeToxAI dataset

df_contact = CSV.read("new_smiles.csv", DataFrame)

contact_tox_mol = [addhydrogens(smilestomol(smiles)) for smiles in df_contact[:, "SMILES"]]

n_mol = length(contact_tox_mol)

n_jobs = Int(n_mol*(n_mol-1)/2 + n_mol)

pbar = Progress(n_jobs, 1)

id_pair = 0
K = zeros(n_mol, n_mol)
for m = 1:n_mol
	for n = m:n_mol
		# compute grw_kernel between every two different cages
		dpg = direct_product_graph(contact_tox_mol[m], contact_tox_mol[n], verbose = true)
		K[m, n] = fixed_point_grw_kernel(dpg, 0.1, ϵ = 0.01)
		K[n, m] = K[m, n]
		next!(pbar)
	end
end

jldsave("BeeToxK.jld2"; K)

all_ones = ones(size(K)) / size(K)[1]

tildeK = K - all_ones * K - K * all_ones + all_ones * K * all_ones

λₛ, Qₛ = eigen(tildeK)

barplot(["λ$i" for i = 1:length(λₛ)], abs.(λₛ), title = "eigenvalues of ̃K")