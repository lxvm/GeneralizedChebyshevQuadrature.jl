"""
    ElementaryFunction <: Function

Abstract type representing functions this library can integrate
"""
abstract type ElementaryFunction <: Function end

"""
    Polynomial <: ElementaryFunction
    Polynomial(coefs)

A polynomial defined by an iterator of coefficients. Evaluates to `x -> evalpoly(x, coefs)`.
"""
struct Polynomial{T} <: ElementaryFunction
    coefs::T
end
(p::Polynomial)(x) = evalpoly(x, p.coefs)
function integrate(p::Polynomial, dom::Interval{D}) where {D<:Real}
    a, b = dom.a, dom.b
    apow = one(a)
    bpow = one(b)
    I = zero(eltype(p.coefs))*(bpow - apow)/1
    for (i,c) in enumerate(p.coefs)
        I += c*((bpow *= b) - (apow *= a))/i
    end
    return I
end

"""
    PowerLaw <: ElementaryFunction
    PowerLaw(λ, a)

A shifted power law of the form ``(x-a)^\\lambda``. Evaluates to `x -> (x-a)^λ`
"""
struct PowerLaw{T,A} <: ElementaryFunction
    λ::T
    a::A
end
(p::PowerLaw)(x) = (x-p.a)^p.λ
function integrate(p::PowerLaw, dom::Interval{D}) where {D<:Real}
    a, b = dom.a, dom.b
    return ((a-p.a)^(p.λ+1) - (b-p.a)^(p.λ+1))/(p.λ+1)
end

"""
    Logarithm <: ElementaryFunction
    Logarithm(b, a)

A shifted logarithm of base `b`. Evaluates to `x -> log(b, x-a)`
"""
struct Logarithm{B,A} <: ElementaryFunction
    b::B
    a::A
end
Logarithm(a) = Logarithm(ℯ, a)
(p::Logarithm)(x) = log(p.b, x-p.a)
function integrate(p::Logarithm, dom::Interval{D}) where {D<:Real}
    a, b = dom.a, dom.b
    return ((a-p.a)^(p.λ+1) - (b-p.a)^(p.λ+1))/(p.λ+1)
end

# predefined polynomial spaces

"""
    legendre(n::Integer, [T=Float64])

Return a vector of monomials with the coefficients corresponding to the Legendre
monomials on ``[-1,1]`` up to degree `n`. Polynomial coefficients are of type `T`.
"""
function legendre(n::Integer, ::Type{T}=Float64) where {T}
    p = Vector{Polynomial{Vector{T}}}(undef, n+1)
    for i in 1:n+1
        c = zeros(T, i)
        sgn = 1
        for j in 1:div(i-1,2)+1
            c[i-2j+2] = binomial(i-1, j-1)*binomial(2*(i-1) - 2*(j-1), i-1)/2^(i-1)
            sgn *= -1
        end
        p[i] = Polynomial(c)
    end
    return p
end

"""
    chebyshev(n::Integer, [T=Float64])

Return a vector of monomials with the coefficients corresponding to the Chebyshev
monomials on ``[-1,1]`` up to degree `n`. Polynomial coefficients are of type `T`.
"""
function chebyshev(n::Integer, ::Type{T}=Float64) where {T}
    p = Vector{Polynomial{Vector{T}}}(undef, n+1)
    for i in 1:n+1
        c = zeros(T, i)
        if i == 1
            c[1] = 1
        else
            sgn = 1
            for j in 1:div(i-1,2)+1
                c[i-2j+2] = sgn*binomial(i-j, j-1)/(i-j)*(i-1)*2.0^(i-2j)
                sgn *= -1
            end
        end
        p[i] = Polynomial(c)
    end
    return p
end

# translations of functions

"""
    Translation(a)

Operator on `ElementaryFunction`s and numbers performing `f(x) -> f(x-a)` resp `x -> x-a`.
"""
struct Translation{T}
    a::T
end

(t::Translation)(x::Number) = x-t.a

function (t::Translation)(m::Polynomial)
    c = zeros(promote_type(eltype(m.coefs), typeof(t.a)), length(m.coefs))
    as = Vector{typeof(t.a)}(undef, length(m.coefs))
    for (n, Mn) in enumerate(m.coefs)
        as[n] = n == 1 ? one(t.a) : as[n-1]*(-t.a) # the sign determines the direction of translation
        for k in 1:n
            # TODO: cache the binomial coefficients, then reorder loops to avoid allocating as
            c[k] += Mn * binomial(n-1, k-1) * as[n-k+1]
        end
    end
    return Polynomial(c)
end

(t::Translation)(p::PowerLaw) = PowerLaw(p.λ, p.a+t.a)
(t::Translation)(l::Logarithm) = Logarithm(l.b, l.a+t.a)

Base.inv(t::Translation) = Translation(-t.a)

# sum functions

"""
    SumFun <: ElementaryFunction
    SumFun(f, g)

A container representing the sum of two `ElementaryFunction`s, `f+g`.
"""
struct SumFun{F<:ElementaryFunction,G<:ElementaryFunction} <: ElementaryFunction
    f::F
    g::G
end
(p::SumFun)(x) = p.f(x) + p.g(x)

integrate(p::SumFun, dom) = integrate(p.f, dom) + integrate(p.g, dom)

Base.:+(f::ElementaryFunction, g::ElementaryFunction) = SumFun(f, g)
function Base.:+(f::Polynomial, g::Polynomial)
    c = zeros(promote_type(eltype(f.coefs), eltype(g.coefs)), max(length(f.coefs), length(g.coefs)))
    for (i,fi) in enumerate(f.coefs)
        c[i] += fi
    end
    for (i,gi) in enumerate(g.coefs)
        c[i] += gi
    end
    return Polynomial(c)
end

# product functions


"""
    ProdFun <: ElementaryFunction
    ProdFun(f, g)

A container representing the product of two `ElementaryFunction`s, `f*g`.
"""
struct ProdFun{F<:ElementaryFunction,G<:ElementaryFunction} <: ElementaryFunction
    f::F
    g::G
end
(p::ProdFun)(x) = p.f(x) * p.g(x)

Base.:*(f::ElementaryFunction, g::ElementaryFunction) = ProdFun(f, g)
function Base.:*(f::Polynomial, g::Polynomial)
    c = zeros(promote_type(eltype(f.coefs), eltype(g.coefs)), length(f.coefs)+length(g.coefs)-1)
    for (j,fj) in enumerate(f.coefs), (i,gi) in enumerate(g.coefs)
        c[i+j-1] += fj*gi
    end
    return Polynomial(c)
end

# for any integration more complex than this, consider using SymbolicIntegration.jl

integrate(p::ProdFun, dom) = prodintegrate(p.f, p.g, dom)
function prodintegrate(f::ElementaryFunction, g::SumFun, dom)
    return prodintegrate(f, g.f, dom) + prodintegrate(f, g.g, dom)
end
prodintegrate(g::SumFun, f::ElementaryFunction, dom) = prodintegrate(f, g, dom)

function prodintegrate(p::PowerLaw, m::Polynomial, dom::Interval)
    a, b = dom.a, dom.b
    λ = p.λ
    t = Translation(-p.a)
    ta = inv(t)(float(a))
    tb = inv(t)(float(b))
    tm = t(m)
    apow = iszero(ta) ? ta^zero(λ) : ta^λ
    bpow = iszero(tb) ? tb^zero(λ) : tb^λ

    I = zero(eltype(tm.coefs))*(bpow - apow)/1
    for (i,c) in enumerate(tm.coefs)
        bpow *= tb
        apow *= ta
        # we integrate 1/x in the principal value sense
        I += c*(iszero(λ+i) ? log(tb)-log(ta) : (bpow - apow)/(λ+i))
    end
    return I
end
prodintegrate(m::Polynomial, p::PowerLaw, dom::Interval) = prodintegrate(p, m, dom)

xlogx(x, logx) = iszero(x) ? x : x*logx
function prodintegrate(l::Logarithm, m::Polynomial, dom::Interval)
    a, b = dom.a, dom.b
    t = Translation(-l.a)
    ta = inv(t)(float(a)) # we need float(a) since -oftype(-2.0, 0.0) -> -0.0 but
    tb = inv(t)(float(b)) # oftype(-2, 0.0) -> 0, which picks different logarithm branches
    tm = t(m)
    apow = one(ta)
    bpow = one(tb)
    logta = log(ta)
    logtb = log(tb)
    loga = zero(logta)
    logb = zero(logtb)
    den = log(l.b)

    I = zero(eltype(tm.coefs))*(xlogx(bpow,logb) - xlogx(apow,loga) - bpow + apow)/den
    for (i,c) in enumerate(tm.coefs)
        bpow *= tb
        apow *= ta
        logb += logtb
        loga += logta
        I += c * (xlogx(bpow,logb) - xlogx(apow, loga) - bpow + apow) / ((i^2)*den)
    end
    return I
end
prodintegrate(m::Polynomial, p::Logarithm, dom::Interval) = prodintegrate(p, m, dom)
