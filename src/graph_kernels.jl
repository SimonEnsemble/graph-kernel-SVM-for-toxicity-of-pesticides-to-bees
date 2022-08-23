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

function fixed_length_rw_kernel(dpg::MetaGraph, l::Int64)
    A_x = Matrix(adjacency_matrix(dpg))
    return sum(A_x ^ l)
end

function fixed_length_rw_kernel(molecule_a::GraphMol, molecule_b::GraphMol, l::Int64)
    dpg = direct_product_graph(molecule_a, molecule_b)
    return fixed_length_rw_kernel(dpg, l)
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
