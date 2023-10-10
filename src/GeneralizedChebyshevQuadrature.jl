"""
A package to generate quadrature rules for functions of the form
``\\psi(x)*s(x-a) + \\phi(x)`` where ``\\psi, \\phi`` are smooth and well-represented by
polynomials, whereas ``s`` may be a singular function. Currently the library only implements
the following singular functions: logarithms and power laws.

## Quickstart

```julia
using GeneralizedChebyshevQuadrature
poly = legendre(20)                             # generate a set of Legendre polynomials of up to degree 20
sing = PowerLaw(-0.5, complex(1.0))             # consider a sqrt singularity at 1
space = Space(                                  # the function space
    vcat(poly, sing .* poly),                   # set of functions (just an iterator of functions)
    complex(0.0),                               # range element of the functions (to get right types)
    Interval(-1, 1))                            # domain of the functions
nodes = cospi.(range(-1, -1e-6, length=10^4))   # sample nodes clustered near endpoints
x, w = generalizedchebyshev(space, nodes, 40)   # compute a 40-point quadrature rule
```

Other options for polynomials include `chebyshev`, which appears to give less numerically
stable weights, and other options for singularities include `Logarithm`.
"""
module GeneralizedChebyshevQuadrature

using LinearAlgebra: qr, ColumnNorm

export Interval, Space, Polynomial, PowerLaw, Logarithm, legendre, chebyshev, generalizedchebyshev

include("spaces.jl")
include("functions.jl")
include("generalizedchebyshev.jl")

end # module GeneralizedChebyshevQuadrature
