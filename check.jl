using JLD2, LinearAlgebra, CairoMakie, CSV, DataFrames

K = load("BeeToxK.jld2")["K"]

N = size(K)[1]

C = I(N) - 1/N * ones(N, N)

Kcenter = C*K*C

all_ones = ones(size(K)) / size(K)[1]

tildeK = K - all_ones * K - K * all_ones + all_ones * K * all_ones

isapprox(Kcenter, tildeK)

λ, V = eigen(Kcenter)

λ = reverse(λ)

V = reverse(V, dims = 2)

L = 2

V_L = V[:,1:L]
Z = sqrt.(diagm(λ[1:L])) * transpose(V_L)

f = Figure()
Axis(f[1,1])
barplot!(1:N, λ)
f

df_contact = CSV.read("new_smiles.csv", DataFrame)

beetox = df_contact[:,"Outcome"]
beetox .== "Toxic"

f1 = Figure()
Axis(f1[1,1])
scatter!(Z[1, beetox .== "Toxic"], Z[2, beetox .== "Toxic"])
scatter!(Z[1, beetox .!= "Toxic"], Z[2, beetox .!= "Toxic"])
f1

