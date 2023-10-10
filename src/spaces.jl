"""
    Interval(a, b)

Represent an interval between two numbers `a, b`. Typically an integration interval.
"""
struct Interval{T}
    a::T
    b::T
end
Interval(a, b) = Interval(promote(a, b)...)

"""
    Space(funs, range, domain)

Represent a function space in terms of:
- `funs`: an iterator of `ElementaryFunction`s defined by this package, e.g. `chebyshev(10)`
- `range`: a value with the same type as the range of the functions, e.g. `0.0 + 0.0im`
- `domain`: an `Interval` representing the integration domain
"""
struct Space{F,R,D}
    funs::F
    range::R
    domain::Interval{D}
end

function moments(sp::Space{F,R,D}) where {F,R,D}
    I = zero(R)*zero(D)
    out = Vector{typeof(I)}(undef, length(sp.funs))
    for (i,f) in enumerate(sp.funs)
        out[i] = integrate(f, sp.domain)
    end
    return out
end

function vandermonde(sp::Space{F,R}, nodes::AbstractVector) where {F,R}
    V = Matrix{R}(undef, length(nodes), length(sp.funs))
    for (j,f) in enumerate(sp.funs)
        for (i,x) in enumerate(nodes)
            V[i,j] = f(x)
        end
    end
    return V
end
