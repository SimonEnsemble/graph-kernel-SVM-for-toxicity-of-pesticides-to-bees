module RWK

using Graphs, MolecularGraph, LinearAlgebra, MetaGraphs, ProgressMeter

export direct_product_graph, grw_kernel, fixed_point_grw_kernel, compute_Gram_matrix, centered_Gram_matrix

# for GraphMol, need to confirm if the bond has the same order
function direct_product_graph(mol_a::GraphMol, mol_b::GraphMol; verbose::Bool=false)
    # unpack species, bond orders
    vertex_labels_a, vertex_labels_b = atomsymbol(mol_a), atomsymbol(mol_b)
    n_a, n_b = length(vertex_labels_a), length(vertex_labels_b)
    edge_labels_a, edge_labels_b = bondorder(mol_a), bondorder(mol_b)

    #= 
    construct vertices of the direct product graph, axb.
        each vertex v=(a,b) of axb represents a pair of vertices
          in graph a and graph b with the same species.

        ab_vertex_pair_to_axb_vertex stores the mapping from
           (a, b) -> v
    =#
    axb = MetaGraph(SimpleGraph(0))
    ab_vertex_pair_to_axb_vertex = zeros(Int, n_a, n_b)
    for a = 1:n_a
        for b = 1:n_b
            if vertex_labels_a[a] == vertex_labels_b[b]
                add_vertex!(axb)
                ab_vertex_pair_to_axb_vertex[a, b] = nv(axb)
                set_props!(axb, nv(axb), Dict(:vertex_pair => (a, b)))
            end
        end
    end
    if verbose
        println("\t# nodes in dpg: ", nv(axb))
    end

    #= 
    assign edges between vertices of the direct product graph, axb.
        for speed, do this as a loop over pairs of edges in a and b, which
          are candidates for edges in axb.
        conditions for an edge between v1=(a1, b1) and v2=(a2, b2):
            * v1
            * 
    =#
    for (e_a, (a_1, a_2)) in enumerate(mol_a.edges)
        for (e_b, (b_1, b_2)) in enumerate(mol_b.edges)
            # only a candidate if edge labels are the same
            if edge_labels_a[e_a] != edge_labels_b[e_b]
                continue
            end

            # check candidate edge (v1, v2) in axb:
            #   v1 = (a1, b1)
            #   v2 = (a2, b2)
            v1 = ab_vertex_pair_to_axb_vertex[a_1, b_1]
            v2 = ab_vertex_pair_to_axb_vertex[a_2, b_2]
            if (v1 != 0) && (v2 != 0) # if v1 and v2 are vertices in axb...
                add_edge!(axb, v1, v2)
            end

            # check candidate edge (v1, v2) in axb:
            #   v1 = (a1, b2)
            #   v2 = (a2, b1)
            v1 = ab_vertex_pair_to_axb_vertex[a_1, b_2]
            v2 = ab_vertex_pair_to_axb_vertex[a_2, b_1]
            if (v1 != 0) && (v2 != 0) # if v1 and v2 are vertices in axb...
                add_edge!(axb, v1, v2)
            end
        end
    end
    if verbose
        println("# edges in dpg: ", ne(axb))
    end

    return axb
end

function grw_kernel(dpg::MetaGraph, γ::Float64)
    if γ >= 1 / Δ(dpg)
        error("γ is greater than 1 / Δ(dpg)")
    end
    A_x = Matrix(adjacency_matrix(dpg))
    B = I(size(A_x)[1]) - γ * A_x
    return sum(inv(B))
end

function grw_kernel(molecule_a::GraphMol, molecule_b::GraphMol, γ::Float64)
    dpg = direct_product_graph(molecule_a, molecule_b)
    return grw_kernel(dpg, γ)
end

function rw_kernel_fixed_length(dpg::MetaGraph, l::Int64)
    A_x = Matrix(adjacency_matrix(dpg))
    return sum(A_x ^ l)
end

function rw_kernel_fixed_length(molecule_a::GraphMol, molecule_b::GraphMol, l::Int64)
    dpg = direct_product_graph(molecule_a, molecule_b)
    return rw_kernel_fixed_length(dpg, l)
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

function fixed_point_grw_kernel(molecule_a::GraphMol, molecule_b::GraphMol, γ::Float64; ϵ::Float64=0.001)
    dpg = direct_product_graph(molecule_a, molecule_b)
    return fixed_point_grw_kernel(dpg, γ, ϵ = ϵ)
end

function compute_Gram_matrix(mols::Vector{GraphMol{SmilesAtom, SmilesBond}}, γ::Float64;
                            verbose::Bool=false)
    n_mol = length(mols) # number of molecules
    
    # for progress bar
    n_jobs = Int(n_mol * (n_mol - 1) / 2 + n_mol)
    pbar = Progress(n_jobs, 1)

    K = zeros(n_mol, n_mol) # Gram matrix
    for m = 1:n_mol
        for n = m:n_mol
            dpg = direct_product_graph(mols[m], mols[n], verbose=verbose)

            K[m, n] = grw_kernel(dpg, γ)
            K[n, m] = K[m, n]

            next!(pbar)
        end
    end

    return K
end

function centered_Gram_matrix(K::Matrix{Float64})
    n_mol = size(K)[1]
    C = I(n_mol) - 1 / n_mol * ones(n_mol, n_mol) # centering matrix
    return C * K * C
end

end
