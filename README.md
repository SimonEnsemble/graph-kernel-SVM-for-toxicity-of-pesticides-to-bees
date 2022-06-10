### classifying the toxicity of pesticides to honey bees via support vector machines and the random walk graph kernel

Python and Julia code to reproduce the results of our paper:
> P. Yang, A. Henle, X. Fern, C. Simon. Classifying the toxicity of pesticides to honey bees via support vector machines with random walk graph kernels. _The Journal of Chemical Physics_. 2022. DOI: 10.1063/5.0090573. [link](https://aip.scitation.org/doi/10.1063/5.0090573)

:bee: the raw data from the [BeeTox AI project](https://www.sciencedirect.com/science/article/pii/S2667318521000131) is in the folder, `BeeToxAI Data/`.
:bee: the folder `/src` contains our source code for the random walk graph kernel.
:bee: to reproduce our results:
* run `python RDkit_smiles_convert.py` to convert the SMILES strings used in the BeeToxAI project into Daylight SMILES or OpenSMILES that can be parsed by `MolecularGraph.jl`.
* run `julia compute_Gram_matrices.jl` to pre-compute the Gram matrices, giving the L-RWGK between each pesticide molecule for different L.
* open `kSVM_beetox.jl` in Pluto to train and evaluate the SVMs and produce the plots in the paper.
