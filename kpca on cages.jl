using RWK, Graphs, MetaGraphs, LinearAlgebra, Xtals, Distributed

@everywhere cages = read_xyz.(readdir("all_cages\\normal type", join = true))

@everywhere cages_bonds = [infer_bonds(cage) for cage in cages]

@everywhere K = zeros(length(cages), length(cages))
@distributed for m = 1:length(cages)
	for n = m:length(cages)
		# compute grw_kernel between every two different cages
		K[m, n] = grw_kernel(cages_bonds[m], cages[m].species, cages_bonds[n], cages[n].species, 0.01)
		K[n, m] = K[m, n]
	end
end

K̃ = K - ones(size(K)[1], size(K)[2]) * K / size(K)[1] - K * ones(size(K)[1], size(K)[2]) / size(K)[1] + ones(size(K)[1], size(K)[2]) * K * ones(size(K)[1], size(K)[2]) / size(K)[1]^2

## goal: write a function for pmap.