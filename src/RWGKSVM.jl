module RWGKSVM

using FIGlet, Graphs, LinearAlgebra, MetaGraphs, MolecularGraph

include("direct_product_graph.jl")
include("graph_kernels.jl")
include("misc.jl")

export direct_product_graph, grw_kernel, fixed_point_grw_kernel, centered_Gram_matrix, fixed_length_rw_kernel

end
