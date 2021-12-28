using RWK, Graphs, MetaGraphs, LinearAlgebra, Xtals, Distributed, JLD2, ProgressMeter, CairoMakie

cages = read_xyz.(readdir("all_cages/normal type", join = true))

cages_bonds = [infer_bonds(cage) for cage in cages]

n_cages = length(cages)

n_jobs = Int(n_cages*(n_cages-1)/2 + n_cages)
n_jobs_dpg = Int(74*74)
pbar = Progress(n_jobs_dpg, 1)

K = zeros(n_cages, n_cages)
dpg_sizes = zeros(Int, n_jobs_dpg)
id_pair = 0
for m = 1:n_cages
	for n = 1:n_cages
		# compute grw_kernel between every two different cages
		dpg = direct_product_graph(cages_bonds[m], cages[m].species, cages_bonds[n], cages[n].species, verbose = true)
		id_pair += 1
		dpg_sizes[id_pair] = nv(dpg)
		# K[m, n] = fixed_point_grw_kernel(dpg, 0.1, ϵ = 0.01)
		# K[n, m] = K[m, n]
		# jldsave("K.jld2"; K)
		next!(pbar)
	end
end

dpg_sizes

barplot(dpg_sizes, strokewidth = 1, strokecolor = :blue)

largevector = findall(size -> size > 100000, dpg_sizes)

function findlocation(dpglargevectors::Vector)
	location_vector = Array{Tuple{Int64, Int64}, 1}(undef, length(dpglargevectors))
	for i = 1:length(dpglargevectors)
		location_vector[i] = (div(dpglargevectors[i], 74), dpglargevectors[i]%74)
	end
	return location_vector
end

# print(maximum(dpg_sizes))

location_vector = findlocation(largevector)

remove_cages = [10, 11, 26, 27, 28, 51, 52, 60, 61, 63, 67]

cages_bonds_new = deleteat!(cages_bonds, remove_cages)

dpg_sizes = zeros(Int, n_jobs)
id_pair = 0
for m = 1:n_cages
	for n = m:n_cages
		# compute grw_kernel between every two different cages
		dpg = direct_product_graph(cages_bonds[m], cages[m].species, cages_bonds[n], cages[n].species, verbose = true)
		id_pair += 1
		dpg_sizes[id_pair] = nv(dpg)
		K[m, n] = fixed_point_grw_kernel(dpg, 0.1, ϵ = 0.01)
		K[n, m] = K[m, n]
		jldsave("K.jld2"; K)
		next!(pbar)
	end
end 

all_ones = ones(n_cages, n_cages) / n_cages

K̃ = K - all_ones * K - K * all_ones + all_ones * K * all_ones

## goal: write a function for pmap.

