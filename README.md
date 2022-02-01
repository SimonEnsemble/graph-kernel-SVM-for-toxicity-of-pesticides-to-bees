# graph kernel applied to beetox data set

steps:
* `python RDkit_smiles_convert.py` (requires installation of RDKit)
* `julia compute_Gram_matrices.jl`
* open `kPCA_beetox.jl` in Pluto
* open `kSVM_beetox.jl` in Pluto
