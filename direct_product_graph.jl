### A Pluto.jl notebook ###
# v0.16.0

using Markdown
using InteractiveUtils

# ╔═╡ d7f84110-2228-11ec-0dbf-7983fe8880c9
using PlutoUI, Statistics, LightGraphs, GraphPlot, PyPlot, LinearAlgebra, Colors, MetaGraphs, Markdown

# ╔═╡ 738ecf75-c936-4f21-965c-424352784cf0
md"# Random Walk Kernel
_by Ping Yang_

Last time we have walked through the basic concept of kernels. This time I will focus on **Graph Kernels**. 

> **Graph Kernels:** a graph kernel is a kernel function that computes an inner product on graphs. Graph kernels can be intuitively understood as functions measuring the similarity of pairs of graphs [^1]. 

Specifically speaking, I will focus on **Random Walk Kernel**.

!!! note \"Random Walk\"

	- Walks: A _**walk**_ $w$ in a graph is an ordered sequence of vertices $w=(a, ...,b)$ such that any two subsequent vertices are connected by an edge.
	- Paths: A $(a,b)$ _**path**_ is a walk that starts in $a$ and ends in $b$ with no repeated vertices.
	- Random Walks: From a vertex $i$ randomly jump to any adjacent vertex $j$

## Principle of random walk kernel

- Count the number of common walks in two input graphs $G$ and $H$

- Walks are sequences of nodes that allow repetitions of nodes

> Gärtner et al. [^2] define the feature space of their kernel as the label sequences derived from walks, and propose a method of computation based on the direct (tensor) product graph of two labeled input graphs.

"

# ╔═╡ 93ef8a13-ae88-4d0b-872f-7a41d922dd47
md" 
## Direct Product Graph

Suppose $G$ and $H$ are graphs with
 -  $\mathcal{V}(G)=\{v_1, v_2, ..., v_m\}$
 -  $\mathcal{V}(H)=\{u_1, u_2, ..., u_n\}$
where $n$ and $m$ can be different number.

Then $G\times H$ is the graph with vertex set

```math
\begin{aligned}
\mathcal{V}(G\times H)=&\mathcal{V}(G) \times \mathcal{V}(H)\\
                      =&\{(v_i, u_j)|v_i \in \mathcal{V}(G), u_j \in \mathcal{V}(H)\}
\end{aligned}
```

and $e$ is an edge of $G\times H$ iff $e=((v_i, u_j),(v_k, u_e))$

where: 
$(v_i, v_k)\in \mathcal{E}(G) \land (u_j, u_e)\in \mathcal{E}(H)$

### Toy example for simple situation
"

# ╔═╡ 89e877bf-b281-484d-b4b9-c48f121e291b
begin
	G = SimpleGraph(3)
	H = SimpleGraph(2)
	G_edges = [(1,2), (2,3), (1,3)]
	H_edges = [(1,2), (2,3)]
	for (vi, vj) in G_edges
		add_edge!(G, vi, vj)
	end
	for (vi, vj) in H_edges
		add_edge!(H, vi, vj)
	end
end

# ╔═╡ 7b7b4cc2-0500-4fd9-b798-0c0a5b5045ef
begin
	gplot(G,layout=spring_layout,nodefillc=RGB(1.0, 1.0, 1.0), 
		  nodestrokec=RGB(0.0, 0.0, 0.0),
		  nodestrokelw=0.5,
		  edgestrokec=colorant"gray",
		  nodelabel=["v<sub>$i</sub>" for i = 1:nv(G)])
end

# ╔═╡ 0be35826-826a-43bc-800d-167081b029eb
begin
	gplot(H,layout=spring_layout,nodefillc=RGB(1.0, 1.0, 1.0), 
		  nodestrokec=RGB(0.0, 0.0, 0.0),
		  nodestrokelw=0.5,
		  edgestrokec=colorant"gray",
		  nodelabel=["u<sub>$i</sub>" for i = 1:nv(H)])
end

# ╔═╡ 0050b1dc-9cd0-49cf-ab77-ed6326f56a74
md" Let's see how it looks like in Graph $G\times H$
"

# ╔═╡ ac0a9d22-d7c2-4b4b-9254-082e9454f122
begin
	GxH = SimpleGraph(3*2)
	GxH_edges = [(1,4), (1,6), (2,3), (2,5), (3,6), (4,5)]
	for (vi, vj) in GxH_edges
		add_edge!(GxH, vi, vj)
	end
	GxH
end

# ╔═╡ 75fedce3-50e7-4574-94a1-a315185c6363
begin
	gplot(GxH,layout=spring_layout,nodefillc=RGB(1.0, 1.0, 1.0), 
		  nodestrokec=RGB(0.0, 0.0, 0.0),
		  nodestrokelw=0.5,
		  edgestrokec=colorant"gray",
		  nodelabel=["(v<sub>$i</sub>,u<sub>$j</sub>)" for j=1:nv(H), i=1:nv(G)])
end

# ╔═╡ 6942a6fc-4e7f-412e-96a1-7aa31b3900e7
begin
	GXH = tensor_product(G,H)
	gplot(GXH, layout = spring_layout)
end

# ╔═╡ df04c74e-79e6-46e1-8dd3-c453b372fe79
md"
But most time in the field of cheminformatics, the graph should be labeled.

For two labeled graphs $G = (\mathcal{V}(G), \mathcal{E}(G))$ and $H = (\mathcal{V}(H), \mathcal{E}(H))$ the _**direct product graph is denoted by**_ $G\times H = (\mathcal{V}(G\times H), \mathcal{E}(G \times H))$ _**and defined as**_:

$\mathcal{V}(G\times H) = \{(v_i, u_j)\in\mathcal{V}(G) \times \mathcal{V}(H)|l(v_i)=l(u_j)\}$

```math
\begin{aligned}
\mathcal{E}(G\times H) = \{((v_i, u_j),(v_k, u_e))\in\mathcal{V}(G\times H)|(v_i,v_k)\in\mathcal{E}(G)\\ \land (u_j,u_e)\in\mathcal{E}(H)\land l(v_i,v_k)=l(u_j,u_e)\}
\end{aligned}
```

"

# ╔═╡ e3053f87-ac53-4559-a7a8-f8f9a295fc64
md"### example for labeled graph"

# ╔═╡ 062a2597-b5f6-4c94-b1ab-39239e3e1a78
# struct molecule
struct Molecule
	species::Vector{Symbol}
	graph::SimpleGraph
end

# ╔═╡ 4ba7884f-d149-4c1b-8bc1-0b5c6ebf4a16
begin
	molecule_g = Molecule([:C,:O,:O,:C], SimpleGraph(4))
	molecule_h = Molecule([:C,:C,:O,:O], SimpleGraph(4))
end

# ╔═╡ 81c51d2f-3f99-4716-9394-2b72bdcaf73c
# assign each atom with a color
atom_to_color = Dict(:C => colorant"gray", :O => colorant"red")

# ╔═╡ 87e356dd-7720-481f-8193-28ff77ea54f0
# draw molecule g and h
begin
	# add edge
	g_edges = [(1,2), (1,3), (1,4), (2,4), (3,4)]
	h_edges = [(1,2), (2,3), (2,4)]
	for (vi,vj) in g_edges
		add_edge!(molecule_g.graph, vi, vj)
	end
	for (ui,uj) in h_edges
		add_edge!(molecule_h.graph, ui, uj)
	end
end

# ╔═╡ b0d42085-1bf0-40ca-aff4-432af6da0f0e
begin
	nodelabels_g = [String(species)*"$i" for (i, species) in enumerate(molecule_g.species)]
	vertex_colors_g = [atom_to_color[molecule_g.species[v_i]] for v_i =                                      1:nv(molecule_g.graph)]
	loc_x, loc_y = spring_layout(molecule_g.graph, 4)
	gplot(molecule_g.graph, loc_x, loc_y, nodefillc = vertex_colors_g, nodelabel = nodelabels_g)
end

# ╔═╡ 6eab1489-2e25-48a5-adc9-ad154b5c3418
begin
	nodelabels_h = [String(species)*"$i" for (i, species) in enumerate(molecule_h.species)]
	vertex_colors_h = [atom_to_color[molecule_h.species[v_i]] for v_i =                                      1:nv(molecule_h.graph)]
	loc_x2, loc_y2 = spring_layout(molecule_h.graph, 2)
	gplot(molecule_h.graph, loc_x2, loc_y2, nodefillc = vertex_colors_h, nodelabel = nodelabels_h)
end

# ╔═╡ 3546f9b7-ede7-437f-b0b4-a87dc5eeabf6
md"Now the most important part is to write a function to make this process automatically.
"

# ╔═╡ b4430262-0f57-46fb-8411-241201075d5f
function direct_product_graph(molecule_a::Molecule, molecule_b::Molecule)
	# find how many vertices in the product graph
	nb_vertex = 0
	for i = 1:nv(molecule_a.graph)
		for j = 1:nv(molecule_b.graph)
			if molecule_a.species[i] == molecule_b.species[j]
				nb_vertex += 1
			end
		end
	end
	# New Graph, prepare to add new edges
	axb = MetaGraph(SimpleGraph(nb_vertex))
	# add vertex properties
	k = 1
	if k <= nb_vertex 
		for i = 1:nv(molecule_a.graph)
			for j = 1:nv(molecule_b.graph)
				if molecule_a.species[i] == molecule_b.species[j]
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
			if i !=j && has_edge(molecule_a.graph, get_prop(axb, i, :uv)[1], get_prop(axb, j, :uv)[1]) && has_edge(molecule_b.graph, get_prop(axb, i, :uv)[2], get_prop(axb, j, :uv)[2])
				add_edge!(axb, i, j)
			end
		end
	end
	# add color information
	for k = 1:nb_vertex
		if molecule_a.species[get_prop(axb, k, :uv)[1]] == :C
			set_props!(axb, k, Dict(:color => colorant"gray"))
		else
			set_props!(axb, k, Dict(:color => colorant"red"))
		end
	end
	return axb
end

# ╔═╡ df459344-b8e1-4e70-869e-b41337d835ce
begin
	gxh = direct_product_graph(molecule_g, molecule_h)
	nodelabel = [get_prop(gxh, i, :uv) for i = 1:nv(gxh)]
	combine_colors = [get_prop(gxh, i, :color) for i=1:nv(gxh)]
	loc_xx, loc_yx = spring_layout(gxh, 1)
	gplot(gxh, loc_xx, loc_yx, nodefillc = combine_colors, nodelabel = nodelabel)
end

# ╔═╡ 9cb9ffce-37e2-4335-a22e-153c3e9bbbbb
md"## Direct Product Kernel

> There is a one-to-one correspondence between walks in $G\times H$ and walks in the graphs $G$ and $H$ with the same label sequence [^3].

The _**direct product kernel**_ is then defined as

```math
K_{RW}(G,H)=\sum_{i,j=1}^{|\mathcal{V}|}[\sum_{l=0}^\infty\lambda_lA_\times^l]
```

Where:
 
 $A_\times$ is the adjacency matrix of $G\times H$

 $\lambda = (\lambda_0,\lambda_1,...)$ is a sequence of weights.
 
 $l$ is the length of walks.

!!! note \"fun fact\"

	Powers of Adjacency matrix involved because they count the random walks.

	Element $(i,j)$ of $A^n$ is the number of walks of length $n$ from vertex $i$ to $j$. 
"

# ╔═╡ ed1c371f-9d78-4e16-bdb2-55646ded7c8d
begin
	M = SimpleGraph(5)
	M_edges = [(1,2), (1,3), (2,3), (3,4), (3,5), (1,5)]
	for (vi, vj) in M_edges
		add_edge!(M, vi, vj)
	end
	loc_xm, loc_ym = spring_layout(M, 2)
	gplot(M,loc_xm, loc_ym,
		  edgestrokec=colorant"gray",
		  nodelabel=[i for i = 1:nv(M)])
end

# ╔═╡ 22cf20e8-3a9b-4b9a-bf96-dcc6c4d93765
A_M = Matrix(adjacency_matrix(M))

# ╔═╡ 4b05aed7-69ab-4420-b0c3-5e5c690cc7f0
A_M^2

# ╔═╡ 2b32e0b7-3b8a-41b2-bb97-531d7599dd0a
A_M^3

# ╔═╡ a6d8da12-23b5-42b5-a478-4117e28c9d73
md"## References 
[^1]. [Graph Kernel Wiki](https://en.wikipedia.org/wiki/Graph_kernel)

[^2]. [Gärtner T, Flach P, Wrobel S (2003) On graph kernels: Hardness results and efficient alternatives. In: Learning Theory and
Kernel Machines. Springer. pp 129–143](https://doi.org/10.1007/978-3-540-45167-9_11)

[^3]. [Kriege, N.M., Johansson, F.D. & Morris, C. A survey on graph kernels. Appl Netw Sci 5, 6 (2020). https://doi.org/10.1007/s41109-019-0195-3](https://doi.org/10.1007/s41109-019-0195-3)
"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
GraphPlot = "a2cc645c-3eea-5389-862e-a155d0052231"
LightGraphs = "093fc24a-ae57-5d10-9952-331d41423f4d"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Markdown = "d6f4376e-aef5-505a-96c1-9c027394607a"
MetaGraphs = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
Colors = "~0.12.8"
GraphPlot = "~0.4.4"
LightGraphs = "~1.3.5"
MetaGraphs = "~0.6.8"
PlutoUI = "~0.7.14"
PyPlot = "~2.10.0"
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

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

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

[[GraphPlot]]
deps = ["ArnoldiMethod", "ColorTypes", "Colors", "Compose", "DelimitedFiles", "LightGraphs", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "dd8f15128a91b0079dfe3f4a4a1e190e54ac7164"
uuid = "a2cc645c-3eea-5389-862e-a155d0052231"
version = "0.4.4"

[[HypertextLiteral]]
git-tree-sha1 = "72053798e1be56026b81d4e2682dbe58922e5ec9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.0"

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

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[JLD2]]
deps = ["DataStructures", "FileIO", "MacroTools", "Mmap", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "192934b3e2a94e897ce177423fd6cf7bdf464bce"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.14"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

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
git-tree-sha1 = "a8709b968a1ea6abc2dc1967cb1db6ac9a00dfb6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.5"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlutoUI]]
deps = ["Base64", "Dates", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "d1fb76655a95bf6ea4348d7197b22e889a4375f4"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.14"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "169bb8ea6b1b143c5cf57df6d34d022a7b60c6db"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.92.3"

[[PyPlot]]
deps = ["Colors", "LaTeXStrings", "PyCall", "Sockets", "Test", "VersionParsing"]
git-tree-sha1 = "14c1b795b9d764e1784713941e787e1384268103"
uuid = "d330b81b-6aea-500a-939a-2ce795aea3ee"
version = "2.10.0"

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

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

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
git-tree-sha1 = "80229be1f670524750d905f8fc8148e5a8c4537f"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.0"

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
# ╠═d7f84110-2228-11ec-0dbf-7983fe8880c9
# ╟─738ecf75-c936-4f21-965c-424352784cf0
# ╟─93ef8a13-ae88-4d0b-872f-7a41d922dd47
# ╠═89e877bf-b281-484d-b4b9-c48f121e291b
# ╠═7b7b4cc2-0500-4fd9-b798-0c0a5b5045ef
# ╠═0be35826-826a-43bc-800d-167081b029eb
# ╟─0050b1dc-9cd0-49cf-ab77-ed6326f56a74
# ╠═ac0a9d22-d7c2-4b4b-9254-082e9454f122
# ╠═75fedce3-50e7-4574-94a1-a315185c6363
# ╠═6942a6fc-4e7f-412e-96a1-7aa31b3900e7
# ╟─df04c74e-79e6-46e1-8dd3-c453b372fe79
# ╟─e3053f87-ac53-4559-a7a8-f8f9a295fc64
# ╠═062a2597-b5f6-4c94-b1ab-39239e3e1a78
# ╠═4ba7884f-d149-4c1b-8bc1-0b5c6ebf4a16
# ╠═81c51d2f-3f99-4716-9394-2b72bdcaf73c
# ╠═87e356dd-7720-481f-8193-28ff77ea54f0
# ╠═b0d42085-1bf0-40ca-aff4-432af6da0f0e
# ╠═6eab1489-2e25-48a5-adc9-ad154b5c3418
# ╟─3546f9b7-ede7-437f-b0b4-a87dc5eeabf6
# ╠═b4430262-0f57-46fb-8411-241201075d5f
# ╠═df459344-b8e1-4e70-869e-b41337d835ce
# ╟─9cb9ffce-37e2-4335-a22e-153c3e9bbbbb
# ╠═ed1c371f-9d78-4e16-bdb2-55646ded7c8d
# ╠═22cf20e8-3a9b-4b9a-bf96-dcc6c4d93765
# ╠═4b05aed7-69ab-4420-b0c3-5e5c690cc7f0
# ╠═2b32e0b7-3b8a-41b2-bb97-531d7599dd0a
# ╟─a6d8da12-23b5-42b5-a478-4117e28c9d73
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
