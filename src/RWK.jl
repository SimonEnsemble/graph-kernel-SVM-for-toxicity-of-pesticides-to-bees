module RWK

using Graphs, MetaGraphs, Xtals, MolecularGraph, LinearAlgebra

export direct_product_graph, grw_kernel, fixed_point_grw_kernel

# the fastest method, store the vertex pair as a matrix
function direct_product_graph(graph_a::Union{SimpleGraph, MetaGraph}, 
	                          species_a::Vector{Symbol},
							  graph_b::Union{SimpleGraph, MetaGraph},
							  species_b::Vector{Symbol}; 
							  verbose::Bool=false)
	axb = SimpleGraph(0)
	# use a matrix to store the vertex no. from both graph_a & graph_b, Aᵢⱼ
	ab_vertex_pair_to_axb_vertex = zeros(Int, nv(graph_a), nv(graph_b))
	for a = 1:nv(graph_a)
		for b = 1:nv(graph_b)
			if species_a[a] == species_b[b]
				add_vertex!(axb)
				ab_vertex_pair_to_axb_vertex[a, b] = nv(axb)
			end
		end
	end
	if verbose
		println("# nodes in dpg: ", nv(axb))
	end

	# loop over every edge in graph a&b
	for ed_a in edges(graph_a)
		a_1, a_2 = Tuple(ed_a)
		for ed_b in edges(graph_b)
			b_1, b_2 = Tuple(ed_b)
			# check if a_1&b_1, a_2&b_2 are paired
			if ab_vertex_pair_to_axb_vertex[a_1, b_1] * ab_vertex_pair_to_axb_vertex[a_2, b_2] !== 0
				add_edge!(axb, ab_vertex_pair_to_axb_vertex[a_1, b_1],
				          ab_vertex_pair_to_axb_vertex[a_2, b_2])
			end
			# check if a_1&b_2, a_2&b_1 are paired
			if ab_vertex_pair_to_axb_vertex[a_1, b_2] * ab_vertex_pair_to_axb_vertex[a_2, b_1] !== 0
				add_edge!(axb, ab_vertex_pair_to_axb_vertex[a_1, b_2],
				          ab_vertex_pair_to_axb_vertex[a_2, b_1])
			end
		end
	end
	if verbose
		println("# edges in dpg: ", ne(axb))
	end
	
	return axb
end

# use the same function when we compute the dpg of two crystal
function direct_product_graph(crystal_a::Crystal, 
	                          crystal_b::Crystal; 
							  verbose::Bool=false)
	return direct_product_graph(crystal_a.bonds, 
								crystal_a.atoms.species, 
								crystal_b.bonds, 
								crystal_b.atoms.species, 
								verbose = verbose)
end

# for GraphMol, need to confirm if the bond has the same order
function direct_product_graph(molecule_a::GraphMol, 
	                          molecule_b::GraphMol; 
							  verbose::Bool=false)
	species_a = atomsymbol(molecule_a)
	species_b = atomsymbol(molecule_b)
	n_a = length(species_a)
	n_b = length(species_b)
	edge_a_bond = zeros(Int, n_a, n_a)
	edge_b_bond = zeros(Int, n_b, n_b)
	for (k, (aᵢ, aⱼ)) in enumerate(molecule_a.edges)
		edge_a_bond[aᵢ, aⱼ] = bondorder(molecule_a)[k]
	end
	for (k, (bᵢ, bⱼ)) in enumerate(molecule_b.edges)
		edge_b_bond[bᵢ, bⱼ] = bondorder(molecule_b)[k]
	end

	axb = SimpleGraph(0)
	ab_vertex_pair_to_axb_vertex = zeros(Int, n_a, n_b)
	for a = 1:n_a
		for b = 1:n_b
			if species_a[a] == species_b[b]
				add_vertex!(axb)
				ab_vertex_pair_to_axb_vertex[a, b] = nv(axb)
			end
		end
	end
	if verbose
		println("# nodes in dpg: ", nv(axb))
	end

	for ed_a in molecule_a.edges
		a_1, a_2 = Tuple(ed_a)
		for ed_b in molecule_b.edges
			b_1, b_2 = Tuple(ed_b)
			if ab_vertex_pair_to_axb_vertex[a_1, b_1] * ab_vertex_pair_to_axb_vertex[a_2, b_2] !== 0 && 
				edge_a_bond[a_1, a_2] == edge_b_bond[b_1, b_2]
				add_edge!(axb, ab_vertex_pair_to_axb_vertex[a_1, b_1],
				          ab_vertex_pair_to_axb_vertex[a_2, b_2])
			end
			if ab_vertex_pair_to_axb_vertex[a_1, b_2] * ab_vertex_pair_to_axb_vertex[a_2, b_1] !== 0 && 
				edge_a_bond[a_1, a_2] == edge_b_bond[b_1, b_2]
				add_edge!(axb, ab_vertex_pair_to_axb_vertex[a_1, b_2],
				          ab_vertex_pair_to_axb_vertex[a_2, b_1])
			end
		end
	end
	if verbose
		println("# edges in dpg: ", ne(axb))
	end

	return axb
end

function grw_kernel(dpg::SimpleGraph, γ::Float64)
	if γ >= 1 / Δ(dpg)
		error("γ is greater than 1 / Δ(dpg)")
	end
	A_x = Matrix(adjacency_matrix(dpg))
	B = I(size(A_x)[1]) - γ * A_x
	invB = inv(B)
	return sum(invB)
end

function grw_kernel(crystal_a::Crystal, crystal_b::Crystal, γ::Float64)
	dpg = direct_product_graph(crystal_a, crystal_b)
	return grw_kernel(dpg, γ)
end

function grw_kernel(molecule_a::GraphMol, molecule_b::GraphMol, γ::Float64)
	dpg = direct_product_graph(molecule_a, molecule_b)
	return grw_kernel(dpg, γ)
end

function grw_kernel(graph_a::Union{SimpleGraph, MetaGraph},
					species_a::Vector{Symbol},
					graph_b::Union{SimpleGraph, MetaGraph},
					species_b::Vector{Symbol},
					γ::Float64)
	dpg = direct_product_graph(graph_a, species_a, graph_b, species_b)
	return grw_kernel(dpg, γ)
end

function fixed_point_grw_kernel(A_x::Matrix, γ::Float64; ϵ::Float64=0.001)
    # B = I - γ*A
    # to compute Inverse B 
    n_x = size(A_x)[1]
    y = rand(n_x)
	y_old = rand(n_x)
	one_vector = ones(n_x)
    while true
        y_old .= y
        y .= one_vector + γ*A_x*y
        if norm(y - y_old) < ϵ
            return sum(y)
        end
    end
end

function fixed_point_grw_kernel(dpg::SimpleGraph, γ::Float64; ϵ::Float64=0.001)
	if γ >= 1 / Δ(dpg)
		error("γ is greater than 1 / Δ(dpg)")
	end
	A_x = Matrix(adjacency_matrix(dpg))
	return fixed_point_grw_kernel(A_x, γ, ϵ = ϵ)
end

function fixed_point_grw_kernel(graph_a::Union{SimpleGraph, MetaGraph},
	                            species_a::Vector{Symbol},
								graph_b::Union{SimpleGraph, MetaGraph},
								species_b::Vector{Symbol},
								γ::Float64;
								ϵ::Float64=0.001,
								verbose::Bool=false)
	dpg = direct_product_graph(graph_a, species_a, graph_b, species_b, verbose = verbose)
	return fixed_point_grw_kernel(dpg, γ, ϵ = ϵ)
end

function fixed_point_grw_kernel(crystal_a::Crystal, crystal_b::Crystal, γ::Float64; ϵ::Float64=0.001)
	dpg = direct_product_graph(crystal_a, crystal_b)
	return fixed_point_grw_kernel(dpg, γ, ϵ = ϵ)
end

function fixed_point_grw_kernel(molecule_a::GraphMol, molecule_b::GraphMol, γ::Float64; ϵ::Float64=0.001)
	dpg = direct_product_graph(molecule_a, molecule_b)
	return fixed_point_grw_kernel(dpg, γ, ϵ = ϵ)
end

end