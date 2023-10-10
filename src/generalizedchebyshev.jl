"""
    generalizedchebyshev(sp::Space, nodes, n::Integer)

Return a tuple `(x, w)` of `n` nodes `x` sampled from `nodes` and quadrature weights `w`
that exactly integrate the functions in the space `sp`. For numerical stability of the
weights `w`, the nodes must be sampled densely from the domain of the space. In particular,
for spaces with singularities, the nodes should be clustered about the singularity. To check
the quality of the nodes, compare `sum(abs, w)` to the measure of the domain and see if they
are all positive.
"""
function generalizedchebyshev(sp::Space, nodes, n::Integer)
    @assert length(nodes) >= length(sp.funs) "There must be at least as many nodes as functions in the space"
    V = vandermonde(sp, nodes)
    @assert all(isfinite, V)  "The Vandermonde matrix has infinite values. Check that the nodes do not coincide with singularities."
    m = moments(sp)
    Q, R, p = qr(transpose(V), ColumnNorm())
    idx = sort(@view(p[1:n]))
    x = nodes[idx]
    w = transpose(view(V, idx, :)) \ m
    return x, w
end
