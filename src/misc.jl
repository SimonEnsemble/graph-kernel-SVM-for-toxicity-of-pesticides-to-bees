"""
    cgm = centered_Gram_matrix(K)
Computes the centering matrix `C`` of Gram matrix `K` and returns the centered Gram matrix `CKC`
"""
function centered_Gram_matrix(K::Matrix{Float64})::Matrix{Float64}
    n_mol = size(K)[1]
    C = I(n_mol) - 1 / n_mol * ones(n_mol, n_mol)
    return C * K * C
end


"""
    RWGKSVM.banner()
Prints the stylized ASCII console banner for the package.
"""
function banner()
    FIGlet.render("RWGKSVM", FIGlet.availablefonts()[449])
end
