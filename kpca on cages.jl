using RWK, Graphs, MetaGraphs, LinearAlgebra, Xtals, Distributed, JLD2

@everywhere cages = read_xyz.(readdir("all_cages\\normal type", join = true))

@everywhere cages_bonds = [infer_bonds(cage) for cage in cages]

@everywhere n_cages = length(cages)

@everywhere K = zeros(n_cages, n_cages)
@distributed for m = 1:n_cages
	for n = m:n_cages
		# compute grw_kernel between every two different cages
		K[m, n] = fixed_point_grw_kernel(cages_bonds[m], cages[m].species, cages_bonds[n], cages[n].species, 0.1, ϵ = 0.01)
		K[n, m] = K[m, n]
	end
end

K = r

jldsave("K.jld2"; K)

all_ones = ones(n_cages, n_cages) / n_cages

K̃ = K - all_ones * K - K * all_ones + all_ones * K * all_ones

## goal: write a function for pmap.

