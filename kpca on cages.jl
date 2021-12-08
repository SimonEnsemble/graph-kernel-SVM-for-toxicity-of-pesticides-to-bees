using RWK, Graphs, MetaGraphs, LinearAlgebra, Xtals, Distributed, JLD2, ProgressMeter

cages = read_xyz.(readdir("all_cages\\normal type", join = true))

cages_bonds = [infer_bonds(cage) for cage in cages]

n_cages = length(cages)

n_jobs = Int(n_cages*(n_cages-1)/2 + n_cages)
pbar = Progress(n_jobs, 1)

K = zeros(n_cages, n_cages)
dpg_sizes = zeros(Int, n_jobs)
id_pair = 0
for m = 1:n_cages
	for n = m:n_cages
		# compute grw_kernel between every two different cages
		dpg = direct_product_graph(cages_bonds[m], cages[m].species, cages_bonds[n], cages[n].species, verbose = true)
		id_pair += 1
		dpg_sizes[id_pair] = nv(dpg)
		# K[m, n] = fixed_point_grw_kernel(dpg, 0.1, ϵ = 0.01, verbose = true)
		# K[n, m] = K[m, n]
		next!(pbar)
	end
end
print(maximum(dpg_sizes))

K

jldsave("K.jld2"; K)

all_ones = ones(n_cages, n_cages) / n_cages

K̃ = K - all_ones * K - K * all_ones + all_ones * K * all_ones

## goal: write a function for pmap.

