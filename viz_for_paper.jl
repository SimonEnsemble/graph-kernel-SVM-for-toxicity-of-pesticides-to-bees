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

# ╔═╡ 0d035b37-a9e6-4081-a95a-229cc3b5befa
md"find unique elements for paper"

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
	unique_elements
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
function viz_graph(mol::GraphMol; savename=nothing, append_prime::Bool=false)
    graph = SimpleGraph(length(atomsymbol(mol)))
    for (vᵢ, vⱼ) in mol.edges
        add_edge!(graph, vᵢ, vⱼ) # warning: does not preserve edge order.
    end
	
	# for (i, (vᵢ, vⱼ)) in enumerate(mol.edges)
	# 	println("mol edge ", i, " is ", (vᵢ, vⱼ))
	# 	println("\tspecies: ", atomsymbol(mol)[vᵢ], "-", atomsymbol(mol)[vⱼ])
	# 	println("\taromatic? ", isaromaticbond(mol)[i])
	# end
	
	edgelabels = ["" for b in bondorder(mol)] # in order of graph
	for (i_graph, ed) in enumerate(edges(graph))
		# find which edge this is in the mol.
		for (i_mol, (vᵢ, vⱼ)) in enumerate(mol.edges)
			if ((vᵢ, vⱼ) == (ed.src, ed.dst)) || ((vᵢ, vⱼ) == (ed.dst, ed.src))
				# println("graph edge ", i_graph, " is mol edge ", i_mol)
				# edgelabels[i_graph] = "$((vᵢ, vⱼ))"
				# println("\taromatic?", isaromaticbond(mol)[i_mol])
				edgelabels[i_graph] = "$(bondorder(mol)[i_mol])"
				if isaromaticbond(mol)[i_mol]
					edgelabels[i_graph] = "a"
				end
				break # found it!
			end
		end
	end

	locs_x, locs_y = spring_layout(graph, C=0.25)
	# locs_x, locs_y = circular_layout(graph)

	nodelabels = ["v<sub>$i</sub>" for i = 1:nv(graph)]
	if append_prime
		nodelabels = ["v′<sub>$i</sub>" for i = 1:nv(graph)]
	end
	
	gp = gplot(graph, locs_x, locs_y,
	      nodestrokec=[RGB(rgb.r/255, rgb.g/255, rgb.b/255) for rgb in atomcolor(mol)],
		nodesize=1.0,
		NODESIZE=0.3 / sqrt(nv(graph)),
		  nodefillc=[RGBA(rgb.r/255, rgb.g/255, rgb.b/255, 0.075) for rgb in atomcolor(mol)], #RGBA(1.0, 1.0, 1.0, 0.0),
		  # nodelabelc=RGB(1.0, 1.0, 1.0),
		  NODELABELSIZE=5.0,
		  EDGELABELSIZE=5.0,
          # nodefillc=RGB(1.0,1.0,1.0),
		  EDGELINEWIDTH=15.0/nv(graph),
	      # nodestrokec = colorant"black",
	      nodestrokelw=1,
          nodelabel=nodelabels,
          edgelinewidth=1,
          edgelabel=edgelabels
	)

	if ! isnothing(savename)
		draw(PDF(savename * ".pdf", 13cm, 13cm), gp)
	end
	return gp
end

# ╔═╡ afda4de4-909e-4b82-bc4c-3c3695321ffd
viz_graph(other_mol, savename="mol_347_graph")

# ╔═╡ 22bbe372-62c9-4982-a02c-12081e84c984
viz_graph(mol_a, savename="mol_a_graph")

# ╔═╡ bf7555f6-7fbb-4bcf-800a-558c93ff8d48
viz_graph(mol_b, savename="mol_b_graph", append_prime=true)

# ╔═╡ 3569a3dc-a1d4-4574-9e60-a254c53377f1
axb = direct_product_graph(mol_a, mol_b, 
	store_vertex_pair=true, store_nodelabel=true, store_edgelabel=true)

# ╔═╡ 0b7d99ac-a23d-4ced-8b27-f57741924bc9
begin
	axb_node = [get_prop(axb, i, :label) for i in vertices(axb)]
	axb_colors = []
	axb_fcolors = []
	for v in 1:nv(axb)
		rgb = atomcolor(mol_a)[get_prop(axb, v, :vertex_pair)[1]]
		push!(axb_colors, RGB(rgb.r/255, rgb.g/255, rgb.b/255))
		push!(axb_fcolors, RGBA(rgb.r/255, rgb.g/255, rgb.b/255, 0.075))
	end
	axb_nodepair = [get_prop(axb, i, :vertex_pair) for i in vertices(axb)]
	axb_nodelabel = ["( v<sub>$i</sub>, v′<sub>$j</sub> )" for (i, j) in axb_nodepair]
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
	
	g = gplot(axb, locs_x, locs_y, 
	      nodefillc =axb_fcolors, # change to nodefill to highlight
		# linetype="curve",
	      nodestrokec=axb_colors,
	      nodestrokelw=1,
			# NODELABELSIZE=5.0,
		  EDGELABELSIZE=5.0,
	      NODESIZE=0.7 / sqrt(nv(axb)),
	      nodelabel=axb_nodelabel, # change to orderlist to debug
	      EDGELINEWIDTH=25.0/nv(axb),
	      edgelabel=axb_edgelabel,
	      edgelabelsize=2.0)

	draw(PDF("dpg_example.pdf", 16cm, 16cm), g)
	g
end

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
# ╟─0d035b37-a9e6-4081-a95a-229cc3b5befa
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
# ╠═b907dfea-cab3-4350-a51f-00f6e8cd7670
# ╠═c6b5616b-d1a8-4cf1-b97b-cf0c142b1fa3
# ╠═a521f480-c5b5-4a60-8f9f-52b45ad39dac
# ╠═179dc2b3-59be-4719-9cec-059dfa447301
# ╠═afda4de4-909e-4b82-bc4c-3c3695321ffd
# ╠═22bbe372-62c9-4982-a02c-12081e84c984
# ╠═bf7555f6-7fbb-4bcf-800a-558c93ff8d48
# ╠═eac14468-06d5-42ce-ac86-890ef38c0b00
# ╠═3569a3dc-a1d4-4574-9e60-a254c53377f1
# ╠═0b7d99ac-a23d-4ced-8b27-f57741924bc9
# ╠═4dbcad70-c777-44b9-b838-810d0f258eaf
# ╠═593f6fc2-8494-480e-ac04-c4ccb548795c
# ╠═e49af556-b443-43ea-b319-7d47b74c9e98
# ╠═2f88108a-e954-49ed-88fe-015249dd8bdd
