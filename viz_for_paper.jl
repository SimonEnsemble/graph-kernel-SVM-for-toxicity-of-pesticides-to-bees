### A Pluto.jl notebook ###
# v0.18.1

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

# ╔═╡ ce93e1c7-cccd-46bb-bdbe-7a8f0f75be70
begin
	mol_a = smilestomol(smiles_a)
	mol_b = smilestomol(smiles_b)
end

# ╔═╡ ea275a9a-a7df-49ff-9c71-ff88033ed28a
atomcolor(mol_b)

# ╔═╡ a6e669b2-bdf8-4786-bf30-6b73302cca27
function viz_molecule(mol, savename)
	canvas = SvgCanvas()
	draw2d!(canvas, mol)
	drawatomindex!(canvas, mol, 
		bgcolor=MolecularGraph.Color(255, 255, 255))
	d = tosvg(canvas, 300, 300)

	f = open(savename, "w")
	write(f, d)
	close(f)
	return d
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
run(`inkscape molecule_a.svg --export-pdf=molecule_a.pdf`)

# ╔═╡ 179dc2b3-59be-4719-9cec-059dfa447301
run(`inkscape molecule_b.svg --export-pdf=molecule_b.pdf`)

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
		NODELABELSIZE=5.0,
		EDGELABELSIZE=5.0,
          nodefillc=RGB(1.0,1.0,1.0),
	      # nodestrokec = colorant"black",
	      nodestrokelw=3,
          nodelabel=["$i" for i = 1:nv(graph)],
          edgelinewidth=2,
          edgelabel = edgelabels)
end

# ╔═╡ 74b9375d-b5ee-4288-980a-300e8c204da3
ga = viz_graph(mol_a)

# ╔═╡ 0703740d-c695-4792-ae0a-395a1da7e8f3
draw(PDF("graph_a.pdf", 15cm, 15cm), ga)

# ╔═╡ d134ca8b-f413-45f4-ac48-b3bb1decf089
gb = viz_graph(mol_b)

# ╔═╡ 305d0c8c-e8c1-4a9e-a51c-8630d8689ea8
draw(PDF("graph_b.pdf", 15cm, 15cm), gb)

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
	nlist[1] = 1:7 # first shell
	nlist[2] = 8:17 # second shell
	locs_x, locs_y = shell_layout(axb, nlist)
	g = gplot(axb, locs_x, locs_y, 
	      nodefillc = RGB(1.0,1.0,1.0),
	      nodestrokec=axb_colors,
	      nodestrokelw=1,
			# NODELABELSIZE=5.0,
		  EDGELABELSIZE=5.0,
	      NODESIZE=0.45 / sqrt(nv(axb)),
	      nodelabel=axb_nodepair,
	      EDGELINEWIDTH=15.0/nv(axb),
	      edgelabel=axb_edgelabel,
	      edgelabelsize=2.0)
end

# ╔═╡ ee75cb9c-6cae-4a38-bb1c-57af21f6dfd7
draw(PDF("dpg_example.pdf", 16cm, 16cm), g)

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
# ╠═b73e7f22-cd3e-4e47-8561-220f9151f948
# ╠═079b7e01-25f7-4e3e-90c3-4783c70b4790
# ╠═daea02b2-2f26-4c78-b96f-82a524aff8c4
# ╠═5d67ebae-d960-46d8-9443-10029fcdbee7
# ╠═95a5ed24-e1f7-45c2-95f9-ac3bcfd3b208
# ╠═2db6218b-b74b-4dc4-8579-c6b70f6d3f3c
# ╠═5cd90fd8-abb2-4d8c-8606-250c45811d2b
# ╠═ce93e1c7-cccd-46bb-bdbe-7a8f0f75be70
# ╠═ea275a9a-a7df-49ff-9c71-ff88033ed28a
# ╠═a6e669b2-bdf8-4786-bf30-6b73302cca27
# ╠═b907dfea-cab3-4350-a51f-00f6e8cd7670
# ╠═c6b5616b-d1a8-4cf1-b97b-cf0c142b1fa3
# ╠═a521f480-c5b5-4a60-8f9f-52b45ad39dac
# ╠═179dc2b3-59be-4719-9cec-059dfa447301
# ╠═eac14468-06d5-42ce-ac86-890ef38c0b00
# ╠═74b9375d-b5ee-4288-980a-300e8c204da3
# ╠═0703740d-c695-4792-ae0a-395a1da7e8f3
# ╠═d134ca8b-f413-45f4-ac48-b3bb1decf089
# ╠═305d0c8c-e8c1-4a9e-a51c-8630d8689ea8
# ╠═3569a3dc-a1d4-4574-9e60-a254c53377f1
# ╠═64c7263e-e959-4068-a4d7-9f1c806ee893
# ╠═0b7d99ac-a23d-4ced-8b27-f57741924bc9
# ╠═ee75cb9c-6cae-4a38-bb1c-57af21f6dfd7
# ╠═4dbcad70-c777-44b9-b838-810d0f258eaf
# ╠═593f6fc2-8494-480e-ac04-c4ccb548795c
# ╠═e49af556-b443-43ea-b319-7d47b74c9e98
# ╠═2f88108a-e954-49ed-88fe-015249dd8bdd
