using JLD2, LinearAlgebra, CairoMakie, CSV, DataFrames, UMAP

addhydrogens_flag = true

if addhydrogens_flag
    Kcenter = load("BeeToxKwh.jld2")["Kcenter"]
else
    Kcenter = load("BeeToxK.jld2")["Kcenter"]

λ, V = eigen(Kcenter)

λ = reverse(λ)

V = reverse(V, dims = 2)

# input L value here
L = 5

V_L = V[:,1:L]
Z = sqrt.(diagm(λ[1:L])) * transpose(V_L)

df_contact = CSV.read("new_smiles.csv", DataFrame)

beetox = df_contact[:,"Outcome"]

f = Figure()
Axis(f[1,1])
barplot!(1:length(λ), λ)
f

if L > 2
    embedding = umap(Z, 2, n_neighbors = 3)
elseif L == 2 
    embedding = Z
end

f1 = Figure()
Axis(f1[1,1])
scatter!(embedding[1, beetox .== "Toxic"], embedding[2, beetox .== "Toxic"])
scatter!(embedding[1, beetox .!= "Toxic"], embedding[2, beetox .!= "Toxic"])
f1