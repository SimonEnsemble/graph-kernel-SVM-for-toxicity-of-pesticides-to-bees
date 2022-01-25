using UMAP, LinearAlgebra, CairoMakie, CSV, DataFrames

Kcenter = load("BeeToxK.jld2")["Kcenter"]

Z = umap(Kcenter, 2)
# problem: why everytime running umap will give a different result?

df_contact = CSV.read("new_smiles.csv", DataFrame)

beetox = df_contact[:,"Outcome"]

f1 = Figure()
Axis(f1[1,1])
scatter!(Z[1, beetox .== "Toxic"], Z[2, beetox .== "Toxic"])
scatter!(Z[1, beetox .!= "Toxic"], Z[2, beetox .!= "Toxic"])
f1