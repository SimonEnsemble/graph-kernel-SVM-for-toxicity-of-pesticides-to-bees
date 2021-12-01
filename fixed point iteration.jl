using LinearAlgebra

using RWK

function fixed_point_rwk(A_x::Matrix, γ::Float64; ϵ::Float64=0.001)
    # B = I - γ*A
    # to compute Inverse B 
    n_x = size(A_x)[1]
    y = rand(n_x)
    while true
        y_old = copy(y)
        y = ones(n_x) + γ*A_x*y
        if norm(y - y_old) < ϵ
            return sum(y)
        end
    end
end