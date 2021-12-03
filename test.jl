using Graphs, MetaGraphs, Xtals, Test, LinearAlgebra, RWK

struct Molecule
	species::Vector{Symbol}
	graph::SimpleGraph
end

molecule_a = Molecule([:C, :O, :O, :C], SimpleGraph(4))
molecule_b = Molecule([:C, :C, :O, :O], SimpleGraph(4))
a_edges = [(1,2), (1,3), (1,4), (2,4), (3,4)]
b_edges = [(1,2), (2,3), (2,4)]
for (vᵢ, vⱼ) in a_edges
    add_edge!(molecule_a.graph, vᵢ, vⱼ)
end
for (uᵢ, uⱼ) in b_edges
    add_edge!(molecule_b.graph, uᵢ, uⱼ)
end

axb = SimpleGraph(8)
axb_edges = [(1,8), (2,3), (2,4), (2,5), (2,6), (2,7), (3,8), (4,8), (5,8), (6,8)]
for (i, j) in axb_edges
    add_edge!(axb, i, j)
end
axb

@test axb == direct_product_graph(molecule_a.graph, molecule_a.species, molecule_b.graph, molecule_b.species)

begin
	xtal1 = Crystal("IRMOF-1.cif")
	infer_bonds!(xtal1, true)
	xtal2 = Crystal("SBMOF-1.cif")
	infer_bonds!(xtal2, true)
end

@time x12 = direct_product_graph(xtal1, xtal2)

@time grw_kernel(x12, 0.1)

typeof(x12)

@time fixed_point_rwk(x12, 0.1, ϵ = 0.0001)