using JLD2, LinearAlgebra, UnicodePlots

KDict = load("K.jld2")

K = KDict["K"]

all_ones = ones(size(K)) / size(K)[1]

# find the Kernel matrix K̃

tildeK = K - all_ones * K - K * all_ones + all_ones * K * all_ones

# Find the eigenvalues and eigenvectors of this matrix K̃

λₛ, Qₛ = eigen(tildeK)

λₛ

print(λₛ)

barplot(["λ<sub>$i</sub>" for i = 1:length(λₛ)], abs.(λₛ), title = "eigenvalues of ̃K")
