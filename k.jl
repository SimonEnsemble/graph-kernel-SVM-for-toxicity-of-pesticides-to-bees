using JLD2

KDict = load("K.jld2")

K = KDict["K"]

all_ones = ones(size(K)) / size(K)[1]

K̃ = K - all_ones * K - K * all_ones + all_ones * K * all_ones

