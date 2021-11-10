module RWK

using LinearAlgebra, Graphs, MetaGraphs, Xtals

export direct_product_graph, grw_kernel

begin
	# the fastest method, store the vertex pair as a matrix
	function direct_product_graph(graph_a::Union{SimpleGraph, MetaGraph}, 
		                          species_a::Vector{Symbol},
		                          graph_b::Union{SimpleGraph, MetaGraph},
		                          species_b::Vector{Symbol}; 
		                          verbose::Bool=false)
		axb = SimpleGraph(0)
		# use a matrix to store the vertex no. from both graph_a & graph_b, Aᵢⱼ
		ab_vertex_pair_to_axb_vertex = zeros(Int, nv(graph_a), nv(graph_b))
		for i = 1:nv(graph_a)
			for j = 1:nv(graph_b)
				if species_a[i] == species_b[j]
					add_vertex!(axb)
					ab_vertex_pair_to_axb_vertex[i, j] = nv(axb)
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
				if ab_vertex_pair_to_axb_vertex[a_1, b_1] !== 0 &&
					ab_vertex_pair_to_axb_vertex[a_2, b_2] !== 0
					add_edge!(axb, ab_vertex_pair_to_axb_vertex[a_1, b_1],
					          ab_vertex_pair_to_axb_vertex[a_2, b_2])
				end
				# check if a_1&b_2, a_2&b_1 are paired
				if ab_vertex_pair_to_axb_vertex[a_1, b_2] !== 0 &&
					ab_vertex_pair_to_axb_vertex[a_2, b_1] !== 0
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
end

begin
	function grw_kernel(dpg::SimpleGraph, γ::Float64)
		if γ >= 1 / Δ(dpg)
			error("γ is greater than 1 / Δ(dpg)")
		end
		A = Matrix(adjacency_matrix(dpg))
		B = I(size(A)[1]) - γ * A
		invB = inv(B)
		return sum(invB)
	end
	
	function grw_kernel(crystal_a::Crystal, crystal_b::Crystal, γ::Float64)
		dpg = direct_product_graph(crystal_a, crystal_b)
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
end

end