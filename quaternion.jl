module Quaternions

import Base: convert, zero, show, +, -, /, *, printstyled, promote_rule
import Base: abs, abs2, inv, conj, real, imag
import LinearAlgebra: normalize

export Quaternion, dot, cross, normalize, conj, inv
export abs, abs2, real, imag, proj_i, proj_j, proj_k
export qToMat, matToQ, qToR3, matToR3

# This constructs only Quaternions of the same type, e.g. only Int8 or Float64..
struct Quaternion{T<:Real} <: Number
    r::T
    i::T
    j::T
    k::T
end

# This promotes the entries to the same type and creates a quaternion
Quaternion(r::Real,i::Real,j::Real,k::Real) = Quaternion(promote(r,i,j,k)...)

# Some more constructors
Quaternion(x::Real) = Quaternion(x, zero(x), zero(x), zero(x))

Quaternion(a::Array{T}) where {T<:Real} =
    if length(a)==4
        Quaternion(a[1], a[2], a[3], a[4])
    elseif length(a)==3
        Quaternion(zero(T), a[1], a[2], a[3])
    else
        throw(DomainError)
    end

# Conversion and Promotion
convert(::Type{Quaternion{T}}, x::Real) where {T} =
    Quaternion(convert(T, x), convert(T, 0), convert(T, 0), convert(T, 0))
convert(::Type{Quaternion{T}}, q::Quaternion{T}) where {T <: Real} = q
convert(::Type{Quaternion{T}}, q::Quaternion) where {T} =
    Quaternion(convert(T, q.r), convert(T, q.i), convert(T, q.j), convert(T, q.k))

promote_rule(::Type{Quaternion{T}}, ::Type{T}) where {T <: Real} = Quaternion{T}
promote_rule(::Type{Quaternion}, ::Type{T}) where {T <: Real} = Quaternion
promote_rule(::Type{Quaternion{T}}, ::Type{S}) where {T <: Real, S <: Real} = Quaternion{promote_type(T,S)}
promote_rule(::Type{Quaternion{T}}, ::Type{Quaternion{S}}) where {T <: Real, S <: Real} = Quaternion{promote_type(T,S)}

function show(io::IO, q::Quaternion)
    pm(x) = x < 0 ? " - $(-x)" : " + $x"
    print(io, q.r, pm(q.i), "i", pm(q.j), "j", pm(q.k), "k")
end


### Here begin the functions
conj(z::Quaternion) = Quaternion(z.r, -z.i, -z.j, -z.k)
abs(z::Quaternion) = sqrt(z.r*z.r + z.i*z.i + z.j*z.j + z.k*z.k)
abs2(z::Quaternion) = z.r*z.r + z.i*z.i + z.j*z.j + z.k*z.k
inv(z::Quaternion) = conj(z)/abs2(z)


(+)(q::Quaternion, w::Real) =
    Quaternion(q.r + w, q.i, q.j, q.k)

(+)(w::Real, q::Quaternion) =
        Quaternion(q.r + w, q.i, q.j, q.k)

(+)(q::Quaternion, w::Quaternion) =
    Quaternion(q.r + w.r, q.i + w.i, q.j + w.j, q.k + w.k)

(-)(q::Quaternion) = Quaternion(-q.r, -q.i, -q.j, -q.k)

(-)(q::Quaternion, w::Real) =
    Quaternion(q.r - w, q.i, q.j, q.k)

(-)(w::Real, q::Quaternion) =
        Quaternion(w - q.r, -q.i, -q.j, -q.k)

(-)(q::Quaternion, w::Quaternion) =
    Quaternion(q.r - w.r, q.i - w.i, q.j - w.j, q.k - w.k)

(*)(q::Quaternion, x::Real) = Quaternion(q.r*x, q.i*x, q.j*x, q.k*x)

(*)(x::Real, q::Quaternion) = Quaternion(q.r*x, q.i*x, q.j*x, q.k*x)

(*)(q::Quaternion, w::Quaternion) = Quaternion(q.r * w.r - q.i * w.i - q.j * w.j - q.k * w.k,
                                               q.r * w.i + q.i * w.r + q.j * w.k - q.k * w.j,
                                               q.r * w.j - q.i * w.k + q.j * w.r + q.k * w.i,
                                               q.r * w.k + q.i * w.j - q.j * w.i + q.k * w.r)

(/)(q::Quaternion, x::Real) = Quaternion(q.r/x, q.i/x, q.j/x, q.k/x)

(/)(q::Quaternion, w::Quaternion) = q * inv(w)


real(q::Quaternion) = q.r
imag(q::Quaternion) = q-q.r
proj_i(q::Quaternion) = q.i
proj_j(q::Quaternion) = q.j
proj_k(q::Quaternion) = q.k
to_coords(q::Quaternion) = q.i, q.j, q.k


function qToMat(q::Quaternion{T}) where {T<:Real}
    return [q.r - im*q.k   -q.j - im*q.i;
            q.j - im*q.i    q.r + im*q.k]
end

function qToR3(q::Quaternion{T}) where {T<:Real}
    return [q.i, q.j, q.k]
end

function matToQ(m::Array{Complex{T}, 2}) where {T<:Real}
    return Quaternion(m[1,1].re, -m[1,2].im, m[2,1].re, m[2,2].im)
end

function matToR3(m::Array{Complex{T}, 2}) where {T<:Real}
    return qToR3(matToQ(m))
end

# This should only be used for purely imaginary quaternions
function dot(q1::Quaternion, q2::Quaternion)
    if real(q1)!=0 || real(q2)!=0
        printstyled(color=:yellow, "Warning in dot:
        The quaternions are not purely imaginary\n")
    end
    -0.5*real(q1*q2+q2*q1)::Real
end

# This is also true if there is some real part
function cross(q1::Quaternion, q2::Quaternion)
    0.5*imag(q1*q2-q2*q1)
end

normalize(q::Quaternion) = q/abs(q)

end
