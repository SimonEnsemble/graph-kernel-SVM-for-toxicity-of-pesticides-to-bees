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

# ╔═╡ 8a90d8c2-e10f-4238-b228-5aa02dd8ecf5
begin
	K̃ = Matrix{Float64}(undef, length(cages), length(cages))
	for m = 1:length(cages)
		for n = 1:length(cages)
			# compute grw_kernel between every two different cages
			K̃[m, n] = grw_kernel(cages_bonds[m], cages[m].species, cages_bonds[n], cages[n].species, 0.01)
		end
	end
	K̃
end

# ╔═╡ f6fba150-7ead-4b1d-a8e7-0050092eaef3
#compute the eigenvalues for kernel matrix K
λₛ, Qₛ = eigen(K̃)

# ╔═╡ Cell order:
# ╠═814d5ded-773b-48ca-8854-0dbf294ae98d
# ╠═b7c3d633-0577-4cb7-9a80-4a2de8368557
# ╠═884f2941-537a-4a3b-a90b-3fcddb9db4fc
# ╠═11e970ca-fa5f-49f9-ba16-bd049baa044b
# ╠═8a90d8c2-e10f-4238-b228-5aa02dd8ecf5
# ╠═f6fba150-7ead-4b1d-a8e7-0050092eaef3
