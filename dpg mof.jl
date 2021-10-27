### A Pluto.jl notebook ###
# v0.16.3

using Markdown
using InteractiveUtils

# ╔═╡ 6380053e-3103-11ec-1fb0-e77e1b25c9d6
using Xtals,GraphPlot, MetaGraphs, Colors, LinearAlgebra , PlutoUI, Graphs

# ╔═╡ 9ca0bdec-0d4a-42bc-9cdf-bf3edcb0e053
rc[:paths][:crystals]

# ╔═╡ 68b40fac-2830-42f7-b886-3f18b47c0bb0
begin
	xtal1 = Crystal("IRMOF-1.cif")
	infer_bonds!(xtal1, true)
end

# ╔═╡ a694954a-63ac-4a6c-8c22-17d5e6ee7682
begin
	xtal2 = Crystal("SBMOF-1.cif")
	infer_bonds!(xtal2, true)
end

# ╔═╡ c9a4e9dd-9a99-42ea-baec-647fc2370f99
xtal1.bonds

# ╔═╡ 97085c0c-984d-4e88-8f72-c4b825c4b6ba
xtal2.bonds

# ╔═╡ b8d86f2c-1c5c-4355-a6eb-cdc315dea337
toy_graphs = MetaGraph(SimpleGraph(3))

# ╔═╡ baf83f29-68fe-45ec-bc2a-d91607d8dfce
gplot(xtal1.bonds)

# ╔═╡ 2eccdc1f-9caa-44fd-93fe-3560b84ddf75
function viz_graph(xtal::Crystal)
	gplot(SimpleGraph(xtal.bonds), 
	  nodefillc=[RGB((Xtals.DEFAULT_CPK_COLORS[a] ./ 255)...) for a in xtal.atoms.species], 
	  nodestrokec=colorant"black",
      nodestrokelw=3)
end

# ╔═╡ 9517ea54-f96a-4e35-b2ee-427bccb7f996
viz_graph(xtal1)

# ╔═╡ d77f0972-7370-40c5-9801-2f044f81938d
viz_graph(xtal2)

# ╔═╡ 350639df-fe11-4742-abae-4df3fab4e004
function direct_product_graph(molecule_a, molecule_b)
	# find how many vertices in the product graph
	nb_vertex = 0
	for i = 1:nv(molecule_a.bonds)
		for j = 1:nv(molecule_b.bonds)
			if molecule_a.atoms.species[i] == molecule_b.atoms.species[j]
				nb_vertex += 1
			end
		end
	end
	# New Graph, prepare to add new edges
	axb = MetaGraph(SimpleGraph(nb_vertex))
	# add vertex properties
	k = 1
	if k <= nb_vertex 
		for i = 1:nv(molecule_a.bonds)
			for j = 1:nv(molecule_b.bonds)
				if molecule_a.atoms.species[i] == molecule_b.atoms.species[j]
					set_props!(axb, k, Dict(:uv => (i, j)))
					k += 1
				end
			end
		end
	end
	# add edges
	for i = 1:nb_vertex
		for j = 1:nb_vertex
			# every prerequisite should be satisfied
			if i !=j && has_edge(molecule_a.bonds, get_prop(axb, i, :uv)[1], get_prop(axb, j, :uv)[1]) && has_edge(molecule_b.bonds, get_prop(axb, i, :uv)[2], get_prop(axb, j, :uv)[2])
				add_edge!(axb, i, j)
			end
		end
	end
	# # add color information
	# for k = 1:nb_vertex
	# 	if molecule_a.species[get_prop(axb, k, :uv)[1]] == :C
	# 		set_props!(axb, k, Dict(:color => colorant"gray"))
	# 	else
	# 		set_props!(axb, k, Dict(:color => colorant"red"))
	# 	end
	# end
	return axb
end

# ╔═╡ 932c5e16-d6fb-485e-b14c-b0241fcdc8a0
begin
	function dpg_fast(graph_a::SimpleGraph, species_a::Vector{Symbol}, graph_b::SimpleGraph, species_b::Vector{Symbol}; verbose::Bool=false)
		axb = MetaGraph(SimpleGraph(0))
		ab_vertex_pair_to_axb_vertex = Dict{Tuple{Int, Int}, Int}()
		for i = 1:nv(graph_a)
			for j = 1:nv(graph_b)
				if species_a[i] == species_b[j]
					add_vertex!(axb)
					set_props!(axb, nv(axb), Dict(:ij => (i, j)))
					ab_vertex_pair_to_axb_vertex[(i, j)] = nv(axb)
				end
			end
		end
		if verbose
			println("# nodes in dpg: ", nv(axb))
		end
		
		# TODO to make faster.
		for ed_a in edges(graph_a)
			i1, i2 = Tuple(ed_a)
			for ed_b in edges(graph_b)
				j1, j2 = Tuple(ed_b)
				if haskey(ab_vertex_pair_to_axb_vertex, (i1, j1)) &&
					haskey(ab_vertex_pair_to_axb_vertex, (i2, j2))
					add_edge!(axb, ab_vertex_pair_to_axb_vertex[(i1, j1)],
					          ab_vertex_pair_to_axb_vertex[(i2, j2)])
				end
				if haskey(ab_vertex_pair_to_axb_vertex, (i1, j2)) &&
					haskey(ab_vertex_pair_to_axb_vertex, (i2, j1))
					add_edge!(axb, ab_vertex_pair_to_axb_vertex[(i1, j2)],
					          ab_vertex_pair_to_axb_vertex[(i2, j1)])
				end
			end
		end
		return axb
	end
		
		
		# TODO change this add edges. this way is SLOW b/c u look over pairs of vertices in dpg.
	# 	for i = 1:nv(axb)
	# 		vᵢ, v′ᵢ = get_prop(axb, i, :ij)
	# 		for j = 1:nv(axb)
	# 			vⱼ, v′ⱼ = get_prop(axb, j, :ij)
	# 			# every prerequisite should be satisfied
	# 			if i != j && has_edge(molecule_a.bonds, vᵢ, vⱼ) && has_edge(
	# 					molecule_b.bonds, v′ᵢ, v′ⱼ)
	# 				add_edge!(axb, i, j)
	# 			end
	# 		end
	# 	end
	# 	if verbose
	# 		println("# edges in dpg: ", ne(axb))
	# 	end
	# 	return axb
	# end
	
	function dpg_fast(molecule_a::Crystal, molecule_b::Crystal; verbose::Bool=false)
		return dpg_fast(molecula_a.bonds, molecule_a.atoms.species, molecule_b.bonds, molecule_b.atoms.species, verbose = verbose)
	end
end

# ╔═╡ 07467c7f-4f42-40e7-9095-0d2d8936115a
with_terminal() do
	@time axb = dpg_fast(xtal1, xtal2)
end

# ╔═╡ 75ff87c8-115e-4a1d-8d21-2f830f8c6e77
axb = dpg_fast(xtal1, xtal2)

# ╔═╡ 9b8ed5ba-c748-4ebe-bc20-6563fdd6af47
md"## Random Walk Kernel"

# ╔═╡ b1596f42-7863-4b34-9198-961a4e2c9c85
I(3)

# ╔═╡ cd15d883-3794-4f63-8f4a-9a13706e1613
function dpg_kernel(dpg::MetaGraph, γ::Float64)
	## γ need to be < 1/a, a >= Δ
	A = Matrix(adjacency_matrix(dpg))
	B = I(size(A)[1]) - γ * A
	invB = inv(B)
	return sum(invB)
end

# ╔═╡ cc9c9c62-62dc-4787-adf0-44b9c3486f4c
md" $\Delta(g)$ will Return the maximum degree of vertices in $g$."

# ╔═╡ 8c33e755-b6a2-43b6-95ee-6413d95fe22d
Δ(axb)

# ╔═╡ a2920e9e-5866-4585-b108-933d7db6f771
# dpg_kernel(axb, 0.1)
# takes about 624 s

# ╔═╡ a5ad2868-ac7f-4263-b602-819d1de5763a
md"## Try k-PCA on cages
"

# ╔═╡ f21b12e5-cfb3-44d9-ad8e-0d09449ab2bf
B9 = read_xyz("all_cages/B9.xyz")

# ╔═╡ c4d6f3c1-37be-4167-9348-20c44ec7a863
A11 = read_xyz("all_cages/A11.xyz")

# ╔═╡ d511a2f3-bdc9-454a-839e-10a050f0af3f
begin
	B9_bonds = infer_bonds(B9)
	A11_bonds = infer_bonds(A11)
end

# ╔═╡ 419cf666-8a18-4e3d-8ea4-804b26ea10ff


# ╔═╡ ec49084c-5a87-41d5-a3f6-6bc86f387fca
function viz_atoms_graph(atoms::Atoms)
	gplot(atoms.bonds, 
	  nodefillc=[RGB((Xtals.DEFAULT_CPK_COLORS[a] ./ 255)...) for a in atoms.species], 
	  nodestrokec=colorant"black",
      nodestrokelw=3)
end

# ╔═╡ 650d7090-54a6-49af-97d1-39a093e3dab7
dpg_fast(A11_bonds, B9_bonds)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
GraphPlot = "a2cc645c-3eea-5389-862e-a155d0052231"
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
MetaGraphs = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Xtals = "ede5f01d-793e-4c47-9885-c447d1f18d6d"

[compat]
Colors = "~0.12.8"
GraphPlot = "~0.5.0"
Graphs = "~1.4.1"
MetaGraphs = "~0.6.8"
PlutoUI = "~0.7.16"
Xtals = "~0.3.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "f87e559f87a45bece9c9ed97458d3afe98b1ebb9"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.1.0"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bio3DView]]
deps = ["Requires"]
git-tree-sha1 = "7f472efd9b6af772307dd017f9deeff2a243754f"
uuid = "99c8bb3a-9d13-5280-9740-b4880ed9c598"
version = "0.1.3"

[[CSV]]
deps = ["Dates", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode"]
git-tree-sha1 = "b83aa3f513be680454437a0eee21001607e5d983"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.8.5"

[[CategoricalArrays]]
deps = ["DataAPI", "Future", "JSON", "Missings", "Printf", "Statistics", "StructTypes", "Unicode"]
git-tree-sha1 = "18d7f3e82c1a80dd38c16453b8fd3f0a7db92f23"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.9.7"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "32a2b8af383f11cbb65803883837a149d10dfe8a"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.10.12"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "31d0151f5716b655421d9d75b7fa74cc4e744df2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.39.0"

[[Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "c6461fc7c35a4bb8d00905df7adafcff1fe3a6bc"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.2"

[[Conda]]
deps = ["JSON", "VersionParsing"]
git-tree-sha1 = "299304989a5e6473d985212c28928899c74e9421"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.5.2"

[[Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataFrames]]
deps = ["CategoricalArrays", "Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d50972453ef464ddcebdf489d11885468b7b83a3"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "0.22.7"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "3c041d2ac0a52a12a27af2782b34900d9c3ee68c"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.11.1"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GraphPlot]]
deps = ["ArnoldiMethod", "ColorTypes", "Colors", "Compose", "DelimitedFiles", "Graphs", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "5e51d9d9134ebcfc556b82428521fe92f709e512"
uuid = "a2cc645c-3eea-5389-862e-a155d0052231"
version = "0.5.0"

[[Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "92243c07e786ea3458532e199eb3feee0e7e08eb"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.4.1"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "5efcf53d798efede8fee5b2c8b09284be359bf24"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.2"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLD2]]
deps = ["DataStructures", "FileIO", "MacroTools", "Mmap", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "46b7834ec8165c541b0b5d1c8ba63ec940723ffb"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.15"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LightGraphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "432428df5f360964040ed60418dd5601ecd240b6"
uuid = "093fc24a-ae57-5d10-9952-331d41423f4d"
version = "1.3.5"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "5a5bc6bf062f0f95e62d0fe0a2d99699fed82dd9"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.8"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[MetaGraphs]]
deps = ["JLD2", "LightGraphs", "Random"]
git-tree-sha1 = "81c0488104fb0dc977f38b4acaff81e6e79efc4d"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.6.8"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f8c673ccc215eb50fcadb285f522420e29e69e1c"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "0.4.5"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "bfd7d8c7fd87f04543810d9cbd3995972236ba1b"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.2"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlutoUI]]
deps = ["Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "4c8a7d080daca18545c56f1cac28710c362478f3"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.16"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a193d6ad9c45ada72c14b731a318bedd3c2f00cf"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.3.0"

[[PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "574a6b3ea95f04e8757c0280bb9c29f1a5e35138"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "0.11.1"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "4ba3651d33ef76e24fef6a598b63ffd1c5e1cd17"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.92.5"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "f45b34656397a1f6e729901dc9ef679610bd12b5"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.8"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures", "Random", "Test"]
git-tree-sha1 = "03f5898c9959f8115e30bc7226ada7d0df554ddd"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "0.3.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "d24a825a95a6d98c385001212dc9020d609f2d4f"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.8.1"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "fed34d0e71b91734bf0a7e10eb1bb05296ddbcd0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[VersionParsing]]
git-tree-sha1 = "e575cf85535c7c3292b4d89d89cc29e8c3098e47"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.1"

[[Xtals]]
deps = ["Bio3DView", "CSV", "DataFrames", "JLD2", "LightGraphs", "LinearAlgebra", "Logging", "MetaGraphs", "Printf", "PyCall", "UUIDs"]
git-tree-sha1 = "28da83a0d4b2691e26b2af51a8bf58a20f738b1a"
uuid = "ede5f01d-793e-4c47-9885-c447d1f18d6d"
version = "0.3.7"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═6380053e-3103-11ec-1fb0-e77e1b25c9d6
# ╠═9ca0bdec-0d4a-42bc-9cdf-bf3edcb0e053
# ╠═68b40fac-2830-42f7-b886-3f18b47c0bb0
# ╠═a694954a-63ac-4a6c-8c22-17d5e6ee7682
# ╠═c9a4e9dd-9a99-42ea-baec-647fc2370f99
# ╠═97085c0c-984d-4e88-8f72-c4b825c4b6ba
# ╠═b8d86f2c-1c5c-4355-a6eb-cdc315dea337
# ╠═baf83f29-68fe-45ec-bc2a-d91607d8dfce
# ╠═2eccdc1f-9caa-44fd-93fe-3560b84ddf75
# ╠═9517ea54-f96a-4e35-b2ee-427bccb7f996
# ╠═d77f0972-7370-40c5-9801-2f044f81938d
# ╠═350639df-fe11-4742-abae-4df3fab4e004
# ╠═932c5e16-d6fb-485e-b14c-b0241fcdc8a0
# ╠═07467c7f-4f42-40e7-9095-0d2d8936115a
# ╠═75ff87c8-115e-4a1d-8d21-2f830f8c6e77
# ╟─9b8ed5ba-c748-4ebe-bc20-6563fdd6af47
# ╠═b1596f42-7863-4b34-9198-961a4e2c9c85
# ╠═cd15d883-3794-4f63-8f4a-9a13706e1613
# ╟─cc9c9c62-62dc-4787-adf0-44b9c3486f4c
# ╠═8c33e755-b6a2-43b6-95ee-6413d95fe22d
# ╠═a2920e9e-5866-4585-b108-933d7db6f771
# ╟─a5ad2868-ac7f-4263-b602-819d1de5763a
# ╠═f21b12e5-cfb3-44d9-ad8e-0d09449ab2bf
# ╠═c4d6f3c1-37be-4167-9348-20c44ec7a863
# ╠═d511a2f3-bdc9-454a-839e-10a050f0af3f
# ╠═419cf666-8a18-4e3d-8ea4-804b26ea10ff
# ╠═ec49084c-5a87-41d5-a3f6-6bc86f387fca
# ╠═650d7090-54a6-49af-97d1-39a093e3dab7
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
