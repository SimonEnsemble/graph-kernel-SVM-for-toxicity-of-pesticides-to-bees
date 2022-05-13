### A Pluto.jl notebook ###
# v0.19.3

using Markdown
using InteractiveUtils

# ╔═╡ e85331d0-9732-11ec-350c-e7f3c76bca59
begin
	push!(LOAD_PATH, joinpath(pwd(), "src"))
	import Pkg
	# Pkg.rm("MolecularGraph")
	Pkg.add(url="https://github.com/SimonEnsemble/MolecularGraph.jl")
	using MolecularGraph, CSV, DataFrames, Graphs, MetaGraphs, Colors, GraphPlot, Cairo, Compose, RWK
end

# ╔═╡ 36af2a96-2edc-4163-9b90-c31dc72feb39
DRAW_SETTING[:C_visible] = true

# ╔═╡ b73e7f22-cd3e-4e47-8561-220f9151f948
data = CSV.read("new_smiles.csv", DataFrame)

# ╔═╡ 079b7e01-25f7-4e3e-90c3-4783c70b4790
begin
	unique_elements = []
	for i = 1:nrow(data)
		mol = smilestomol(data[i, "SMILES"])
		for elem in unique(atomsymbol(mol))
			if ! (elem in unique_elements)
				push!(unique_elements, elem)
			end
		end
	end
end

# ╔═╡ daea02b2-2f26-4c78-b96f-82a524aff8c4
begin
	# good pairs: #7 and #134, #35 and #134
	id_a = 35
	id_b = 134
end

# ╔═╡ 5d67ebae-d960-46d8-9443-10029fcdbee7
smiles_a = data[id_a, "SMILES"]

# ╔═╡ 95a5ed24-e1f7-45c2-95f9-ac3bcfd3b208
data[id_a, "Outcome"]

# ╔═╡ 2db6218b-b74b-4dc4-8579-c6b70f6d3f3c
data[id_a, "Outcome"]

# ╔═╡ 5cd90fd8-abb2-4d8c-8606-250c45811d2b
smiles_b = data[id_b, "SMILES"]

# ╔═╡ 03c431ed-d285-4013-8a0e-02e9a7f789c3
data[348, "SMILES"]

# ╔═╡ ce93e1c7-cccd-46bb-bdbe-7a8f0f75be70
begin
	mol_a = smilestomol(smiles_a)
	mol_b = smilestomol(smiles_b)
	other_mol = smilestomol(data[348, "SMILES"])
end

# ╔═╡ ea275a9a-a7df-49ff-9c71-ff88033ed28a
atomcolor(mol_b)

# ╔═╡ a6e669b2-bdf8-4786-bf30-6b73302cca27
function viz_molecule(mol, savename)
	canvas = SvgCanvas()
	draw2d!(canvas, mol)
	drawatomindex!(canvas, mol, 
		bgcolor=MolecularGraph.Color(255, 255, 255), opacity=0.6)
	svg_string = tosvg(canvas, 300, 300)
	savesvg(svg_string, savename)
	return svg_string
end

# ╔═╡ f27b6ba8-72b4-4ff2-bb8e-5bd8892771da
begin
	d = viz_molecule(other_mol, "mol_347.svg")
	HTML(d)
end

# ╔═╡ b907dfea-cab3-4350-a51f-00f6e8cd7670
begin
	da = viz_molecule(mol_a, "molecule_a.svg")
	HTML(da)
end

# ╔═╡ c6b5616b-d1a8-4cf1-b97b-cf0c142b1fa3
begin
	db = viz_molecule(mol_b, "molecule_b.svg")
	HTML(db)
end

# ╔═╡ a521f480-c5b5-4a60-8f9f-52b45ad39dac
run(`inkscape molecule_a.svg --export-pdf=Fig2a_molecule_a.pdf`)

# ╔═╡ 179dc2b3-59be-4719-9cec-059dfa447301
run(`inkscape molecule_b.svg --export-pdf=Fig2b_molecule_b.pdf`)

# ╔═╡ eac14468-06d5-42ce-ac86-890ef38c0b00
function viz_graph(mol::GraphMol)
    graph = SimpleGraph(length(atomsymbol(mol)))
    for (vᵢ, vⱼ) in mol.edges
        add_edge!(graph, vᵢ, vⱼ)
    end
	
    edgelabels = ["$b" for b in bondorder(mol)]
    edgelabels[isaromaticbond(mol)] .= "a"

	locs_x, locs_y = spring_layout(graph, C=0.25)
	# locs_x, locs_y = circular_layout(graph)
	gplot(graph, locs_x, locs_y,
	      nodestrokec=[RGB(rgb.r/255, rgb.g/255, rgb.b/255) for rgb in atomcolor(mol)],
		  nodefillc=RGBA(1.0, 1.0, 1.0, 0.0),
		NODELABELSIZE=5.0,
		EDGELABELSIZE=5.0,
          # nodefillc=RGB(1.0,1.0,1.0),
		EDGELINEWIDTH=20.0/nv(graph),
	      # nodestrokec = colorant"black",
	      nodestrokelw=3,
          nodelabel=["$i" for i = 1:nv(graph)],
          edgelinewidth=2,
          edgelabel = edgelabels)
end

# ╔═╡ a3f80472-8d74-4bad-b0a8-6a32cef149ed
viz_graph(other_mol)

# ╔═╡ 26ac6c3c-6f60-46b7-af9d-b4c3c35295a8
begin
	graph_a = SimpleGraph(length(atomsymbol(mol_a)))
    for (vᵢ, vⱼ) in mol_a.edges
        add_edge!(graph_a, vᵢ, vⱼ)
    end
	
    edgelabels = ["$b" for b in bondorder(mol_a)]
    edgelabels[isaromaticbond(mol_a)] .= "a"

	nodefill_a = []
	for i = 1:length(atomsymbol(mol_a))
		if 1<=i<=2 || 5<=i<=6 || 9<=i<=13
			push!(nodefill_a, RGBA(1.0,1.0,1.0,0.0))
		else
			push!(nodefill_a, RGB(1.0,0.988,0.0))
		end
	end

	locs_x_a, locs_y_a = spring_layout(graph_a, C=0.25)
	# locs_x, locs_y = circular_layout(graph)
	gga = gplot(graph_a, locs_x_a, locs_y_a,
	      nodestrokec=[RGB(rgb.r/255, rgb.g/255, rgb.b/255) for rgb in atomcolor(mol_a)],
		  nodefillc=nodefill_a,
		NODELABELSIZE=5.0,
		EDGELABELSIZE=5.0,
          # nodefillc=RGB(1.0,1.0,1.0),
		EDGELINEWIDTH=20.0/nv(graph_a),
	      # nodestrokec = colorant"black",
	      nodestrokelw=3,
          nodelabel=["$i" for i = 1:nv(graph_a)],
          edgelinewidth=2,
          edgelabel = edgelabels)
end

# ╔═╡ b89d2a05-b5d7-4aba-8a62-4c01e8ca59b0
draw(PDF("Fig2c_graph_a_hl.pdf", 20cm, 20cm), gga)

# ╔═╡ d7912d52-f141-4ed7-ac18-39a588891de6
begin
	graph_b = SimpleGraph(length(atomsymbol(mol_b)))
    for (vᵢ, vⱼ) in mol_b.edges
        add_edge!(graph_b, vᵢ, vⱼ)
    end
	
    edgelabels_b = ["$b" for b in bondorder(mol_b)]
    edgelabels_b[isaromaticbond(mol_b)] .= "a"

	nodefill_b = []
	for i = 1:length(atomsymbol(mol_b))
		if 1<=i<=3 || i == 5
			push!(nodefill_b, RGB(1.0,0.988,0.0))
		else
			push!(nodefill_b, RGB(1.0,1.0,1.0))
		end
	end

	locs_x_b, locs_y_b = spring_layout(graph_b, C=0.25)
	# locs_x, locs_y = circular_layout(graph)
	ggb = gplot(graph_b, locs_x_b, locs_y_b,
	      nodestrokec=[RGB(rgb.r/255, rgb.g/255, rgb.b/255) for rgb in atomcolor(mol_b)],
		  nodefillc=nodefill_b,
		NODELABELSIZE=5.0,
		EDGELABELSIZE=5.0,
          # nodefillc=RGB(1.0,1.0,1.0),
		EDGELINEWIDTH=20.0/nv(graph_b),
	      # nodestrokec = colorant"black",
	      nodestrokelw=3,
          nodelabel=["$i" for i = 1:nv(graph_b)],
          edgelinewidth=2,
          edgelabel = edgelabels_b)
end

# ╔═╡ 82d4d058-1ac5-4c2a-b85a-a45ded8fe082
draw(PDF("Fig2d_graph_b_hl.pdf", 15cm, 15cm), ggb)

# ╔═╡ 3569a3dc-a1d4-4574-9e60-a254c53377f1
axb = direct_product_graph(mol_a, mol_b, 
	store_vertex_pair=true, store_nodelabel=true, store_edgelabel=true)

# ╔═╡ 64c7263e-e959-4068-a4d7-9f1c806ee893
w = [(4, 5), (3, 3), (2, 2), (1, 1)]

# ╔═╡ 0b7d99ac-a23d-4ced-8b27-f57741924bc9
begin
	axb_node = [get_prop(axb, i, :label) for i in vertices(axb)]
	axb_colors = []
	for v in 1:nv(axb)
		rgb = atomcolor(mol_a)[get_prop(axb, v, :vertex_pair)[1]]
		push!(axb_colors, RGB(rgb.r/255, rgb.g/255, rgb.b/255))
	end
	axb_nodepair = [get_prop(axb, i, :vertex_pair) for i in vertices(axb)]
	axb_edgelabel = [get_prop(axb, i, :bondorder) for i in 1:ne(axb)]
	order_list = ["$i" for i = 1:nv(axb)]
	nlist = Vector{Vector{Int}}(undef, 2) # two shells
	nlist[1] = 1:6 # first shell
	nlist[2] = 7:17 # second shell
	@assert length(unique(vcat(nlist[1], nlist[2]))) == nv(axb)
	_locs_x, _locs_y = shell_layout(axb, nlist)
	locs_x = deepcopy(_locs_x)
	locs_y = deepcopy(_locs_y)
	# swap
	locs_x[12] = _locs_x[13]
	locs_x[13] = _locs_x[12]
	locs_y[12] = _locs_y[13]
	locs_y[13] = _locs_y[12]
	
	locs_x[9] = _locs_x[10]
	locs_x[10] = _locs_x[9]
	locs_y[9] = _locs_y[10]
	locs_y[10] = _locs_y[9]

	locs_x[1] = _locs_x[2]
	locs_x[2] = _locs_x[1]
	locs_y[1] = _locs_y[2]
	locs_y[2] = _locs_y[1]
	
	locs_y[6] -= 0.1
	locs_y[7] -= 0.1
	
	# for highlight a walk example
	nodefill = []
	for i = 1:17
		if i == 14|| i == 12|| i == 5|| i == 7
			push!(nodefill, RGB(1.0, 0.988, 0.0))
		else
			push!(nodefill, RGBA(1.0,1.0,1.0,0.0))
		end
	end

	g = gplot(axb, locs_x, locs_y, 
	      nodefillc = nodefill, # change to nodefill to highlight
		# linetype="curve",
	      nodestrokec=axb_colors,
	      nodestrokelw=1,
			# NODELABELSIZE=5.0,
		  EDGELABELSIZE=5.0,
	      NODESIZE=0.45 / sqrt(nv(axb)),
	      nodelabel=axb_nodepair, # change to orderlist to debug
	      EDGELINEWIDTH=25.0/nv(axb),
	      edgelabel=axb_edgelabel,
	      edgelabelsize=2.0)
end

# ╔═╡ ee75cb9c-6cae-4a38-bb1c-57af21f6dfd7
draw(PDF("Fig2e_dpg_example.pdf", 16cm, 16cm), g)

# ╔═╡ 4dbcad70-c777-44b9-b838-810d0f258eaf
fixed_length_rw_kernel(mol_a, mol_b, 0)

# ╔═╡ 593f6fc2-8494-480e-ac04-c4ccb548795c
nv(axb)

# ╔═╡ e49af556-b443-43ea-b319-7d47b74c9e98
fixed_length_rw_kernel(mol_a, mol_b, 1)

# ╔═╡ 2f88108a-e954-49ed-88fe-015249dd8bdd
ne(axb)

# ╔═╡ Cell order:
# ╠═e85331d0-9732-11ec-350c-e7f3c76bca59
# ╠═36af2a96-2edc-4163-9b90-c31dc72feb39
# ╠═b73e7f22-cd3e-4e47-8561-220f9151f948
# ╠═079b7e01-25f7-4e3e-90c3-4783c70b4790
# ╠═daea02b2-2f26-4c78-b96f-82a524aff8c4
# ╠═5d67ebae-d960-46d8-9443-10029fcdbee7
# ╠═95a5ed24-e1f7-45c2-95f9-ac3bcfd3b208
# ╠═2db6218b-b74b-4dc4-8579-c6b70f6d3f3c
# ╠═5cd90fd8-abb2-4d8c-8606-250c45811d2b
# ╠═03c431ed-d285-4013-8a0e-02e9a7f789c3
# ╠═ce93e1c7-cccd-46bb-bdbe-7a8f0f75be70
# ╠═ea275a9a-a7df-49ff-9c71-ff88033ed28a
# ╠═a6e669b2-bdf8-4786-bf30-6b73302cca27
# ╠═f27b6ba8-72b4-4ff2-bb8e-5bd8892771da
# ╠═a3f80472-8d74-4bad-b0a8-6a32cef149ed
# ╠═b907dfea-cab3-4350-a51f-00f6e8cd7670
# ╠═c6b5616b-d1a8-4cf1-b97b-cf0c142b1fa3
# ╠═a521f480-c5b5-4a60-8f9f-52b45ad39dac
# ╠═179dc2b3-59be-4719-9cec-059dfa447301
# ╠═eac14468-06d5-42ce-ac86-890ef38c0b00
# ╠═26ac6c3c-6f60-46b7-af9d-b4c3c35295a8
# ╠═b89d2a05-b5d7-4aba-8a62-4c01e8ca59b0
# ╠═d7912d52-f141-4ed7-ac18-39a588891de6
# ╠═82d4d058-1ac5-4c2a-b85a-a45ded8fe082
# ╠═3569a3dc-a1d4-4574-9e60-a254c53377f1
# ╠═64c7263e-e959-4068-a4d7-9f1c806ee893
# ╠═0b7d99ac-a23d-4ced-8b27-f57741924bc9
# ╠═ee75cb9c-6cae-4a38-bb1c-57af21f6dfd7
# ╠═4dbcad70-c777-44b9-b838-810d0f258eaf
# ╠═593f6fc2-8494-480e-ac04-c4ccb548795c
# ╠═e49af556-b443-43ea-b319-7d47b74c9e98
# ╠═2f88108a-e954-49ed-88fe-015249dd8bdd
