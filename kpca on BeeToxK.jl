using JLD2, LinearAlgebra, CairoMakie, CSV, DataFrames

Kcenter = load("BeeToxK.jld2")["Kcenter"]

λ, V = eigen(Kcenter)

λ = reverse(λ)

V = reverse(V, dims = 2)

L = 2

V_L = V[:,1:L]
Z = sqrt.(diagm(λ[1:L])) * transpose(V_L)

df_contact = CSV.read("new_smiles.csv", DataFrame)

beetox = df_contact[:,"Outcome"]

f = Figure()
Axis(f[1,1])
barplot!(1:N, λ)
f

f1 = Figure()
Axis(f1[1,1])
scatter!(Z[1, beetox .== "Toxic"], Z[2, beetox .== "Toxic"])
scatter!(Z[1, beetox .!= "Toxic"], Z[2, beetox .!= "Toxic"])
f1