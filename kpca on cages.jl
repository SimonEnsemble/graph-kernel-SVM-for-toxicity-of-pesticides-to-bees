### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 814d5ded-773b-48ca-8854-0dbf294ae98d
begin
	import Pkg
	Pkg.develop("RWK")
end

# ╔═╡ b7c3d633-0577-4cb7-9a80-4a2de8368557
using RWK, Graphs, MetaGraphs, LinearAlgebra, Xtals

# ╔═╡ 884f2941-537a-4a3b-a90b-3fcddb9db4fc
cages = read_xyz.(readdir("all_cages\\normal type", join = true))

# ╔═╡ 11e970ca-fa5f-49f9-ba16-bd049baa044b
begin
	cages_bonds = []
	for i = 1:length(cages)
		push!(cages_bonds, infer_bonds(cages[i]))
	end
	cages_bonds
end

# ╔═╡ 4148b293-fad9-4d70-accf-dcd058977684
direct_product_graph(cages_bonds[1], cages[1].species, cages_bonds[2], cages[2].species)

# ╔═╡ 1f2413a5-04c6-4730-9207-a4cb96a6dd6a
direct_product_graph(cages_bonds[1], cages[1].species, cages_bonds[3], cages[3].species)

# ╔═╡ 8a90d8c2-e10f-4238-b228-5aa02dd8ecf5
# begin
# 	K̃ = Matrix{Float64}(undef, length(cages), length(cages))
# 	for m = 1:length(cages)
# 		for n = m:length(cages)
# 			# compute grw_kernel between every two different cages
# 			K̃[m, n] = grw_kernel(cages_bonds[m], cages[m].species, cages_bonds[n], cages[n].species, 0.01)
#           K̃[n, m] = K̃[m, n]
# 		end
# 	end
# 	K̃
# end

# ╔═╡ f6fba150-7ead-4b1d-a8e7-0050092eaef3
#compute the eigenvalues for kernel matrix K
λₛ, Qₛ = eigen(K̃)

# ╔═╡ Cell order:
# ╠═814d5ded-773b-48ca-8854-0dbf294ae98d
# ╠═b7c3d633-0577-4cb7-9a80-4a2de8368557
# ╠═884f2941-537a-4a3b-a90b-3fcddb9db4fc
# ╠═11e970ca-fa5f-49f9-ba16-bd049baa044b
# ╠═4148b293-fad9-4d70-accf-dcd058977684
# ╠═1f2413a5-04c6-4730-9207-a4cb96a6dd6a
# ╠═8a90d8c2-e10f-4238-b228-5aa02dd8ecf5
# ╠═f6fba150-7ead-4b1d-a8e7-0050092eaef3
