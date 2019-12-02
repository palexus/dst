include("combinatoric.jl")
include("quaternion.jl")

using .Quaternions
using LinearAlgebra: det

function nextN(n::Quaternion, n1::Quaternion, n2::Quaternion)::Quaternion
    dot(n, n1 + n2) / (1 + dot(n1, n2)) * (n1 + n2) - n
end

function propagate(
    verlst::Array{Quaternion{T},1},
    horlst::Array{Quaternion{S},1},
    func,
) where {T<:Real,S<:Real}
    if verlst[1] != horlst[1]
        println("The initial value n0 must be the same in the the
        vertical and the horizontal list")
    else
        r, s = length(verlst), length(horlst)
        output = zeros(Quaternion, (r, s))
        output[:, 1] = verlst
        output[1, :] = horlst
        for i = 2:r, j = 2:s
            output[i, j] = func(
                output[i-1, j-1],
                output[i, j-1],
                output[i-1, j],
            )
        end
        return output
    end
end

function propagate_zigzag(
    layer1::Array{Quaternion{T}, 1},
    layer2::Array{Quaternion{S}, 1},
    func::Function,
    n::Int,
) where {T,S <:Real}
    m1 = length(layer1)
    m2 = length(layer2)
    if abs(m1-m2)>1
        println("Dimension mismatch for zigzag")
        throw(DomainError)
    end
    if m1!=m2
        println("Up to now propagate_zigzag is only
        implemented for periodic start conditions")
        return
    end
    output = zeros(Quaternion, (m1, n+2))
    #output = Array{Any, 2}(undef, m1, n+2)
    output[:, 1] = layer1
    output[:, 2] = layer2
    for j=1:n, i=1:m1
        output[i, j+2] = func(
            output[i, j],
            output[no(i, j, m1)...],
            output[nw(i, j, m1)...],
        )
    end
    output
end

function build_gauss(
    verlst::Array{Quaternion{T},1},
    horlst::Array{Quaternion{S},1},
) where {T<:Real,S<:Real}
    propagate(verlst, horlst, nextN)
end

function sample_small_circle(n::Quaternion, height::Real, d::Int; shift=false)
    n = normalize(n)
    radius = sqrt(1-height^2)
    if n.i!=0 || n.j!=0
        e1 = normalize(Quaternion([-n.j, n.i, 0]))
    else
        e1 = Quaternion([1, 0, 0])
    end
    e2 = cross(n, e1)
    if shift
        circle = [radius*(cos(2*pi/d*(i-1)+pi/d)*e1 + sin(2*pi/d*(i-1)+pi/d)*e2)+height*n
            for i = 1:d]
    else
        circle = [radius*(cos(2*pi/d*(i-1))*e1 + sin(2*pi/d*(i-1))*e2)+height*n
            for i = 1:d]
    end
    return circle
end



function nextF(
    n::Quaternion{T},
    n1::Quaternion{S},
    f::Quaternion{U}
) where {S,T,U<:Real}
    return f + cross(n, n1)
end


function knet!(g::MetaDiGraph; f0=zero(Quaternion{Float64}))
    set_prop!(g, 1, :surface, f0)
    trav = treeIter(trav_tree(g))
    for tup in trav
        src, dst = tup
        e = Edge(src, dst)
        f = get_prop(g, src, :surface)
        n = get_prop(g, src, :gauss)
        n1 = get_prop(g, dst, :gauss)
        if e in edges(g)
            if get_prop(g, e, :dir)=="right"
                set_prop!(g, dst, :surface, nextF(n1, n, f))
            else
                set_prop!(g, dst, :surface, nextF(n, n1, f))
            end
        else
            e = reverse(e)
            if get_prop(g, e, :dir)=="right"
                set_prop!(g, dst, :surface, nextF(n1, n, f))
            else
                set_prop!(g, dst, :surface, nextF(n, n1, f))
            end
        end
    end
end


function h_on_first_row!(g::MetaDiGraph)
    m, n = get_prop(g, :dim)
    S, T = 0, 0
    if get_prop(g, :periodic)==true
        S = 0.5*sum([-(-1)^k*get_prop(g, m+k, :iangle) for k=1:m])
        T = 0.5*sum([-(-1)^k*get_prop(g, k, :oangle) for k=1:m])
    end
    set_prop!(g, 1, :h, exp(im*S))
    set_prop!(g, m+1, :h, exp(im*T))
    for v=2:m
        out = sorted_outneighbors(g, v)
        Q = exp(im*get_prop(g, out[1], :iangle))
        H = get_prop(g, v-1, :h)
        H1 = Q/H  #Q/H
        set_prop!(g, v, :h, H1)
    end
end


function setup_h2!(g::MetaDiGraph)
    for v in vertices(g)
        v1 = right(g, v)
        v2 = left(g, v)
        v12 = opposite_vertex(g, v)
        if v1!=nothing && v2!=nothing
            Q = exp(im*get_prop(g, v, :oangle))
            if haskey(props(g,v1),:h) && !haskey(props(g,v2),:h)
                H1 = get_prop(g, v1, :h)
                H2 = Q/H1  #Q/H1
                set_prop!(g, v2, :h, H2)
            elseif haskey(props(g,v2),:h) && !haskey(props(g,v1),:h)
                H2 = get_prop(g, v2, :h)
                H1 = Q/H2
                set_prop!(g, v1, :h, H1)
            end
        end
        if v12!=nothing && v1!=nothing && v2!=nothing
            Q = exp(im*get_prop(g, v, :oangle))
            H = get_prop(g, v, :h)
            deltau = get_prop(g, v, v1, :delta)
            deltav = get_prop(g, v1, v12, :delta)
            k = tan(deltau/2)*tan(deltav/2)
            R = -(k*Q+1)/(Q+k)
            H12 = -1/H/R
            set_prop!(g, v12, :h, H12)
        end
    end
end


function setup_h!(g::MetaDiGraph)
    m, n = get_prop(g, :dim)
    set_eprops!(g, spherical_dist, :gauss, :delta)
    set_vprops_out!(g, spherical_angle, :gauss, :oangle)
    if get_prop(g, :type)=="zigzag"
        set_vprops_inn!(g, spherical_angle, :gauss, :iangle)
        h_on_first_row!(g)
    else
        set_prop!(g, 1, :h, 1)
        set_prop!(g, m+1, :h, exp(im*get_prop(g, 1, :oangle)/2))
    end
    setup_h2!(g)
end


function uMatrix0(delta, dH)
    [cot(delta/2)*dH im;im cot(delta/2)/dH]
end


function vMatrix0(delta, Hij)
    [1 im*tan(delta/2)*Hij;
    im*tan(delta/2)/Hij 1]
end


function setup_lax!(g::MetaDiGraph)
    if !haskey(props(g,1), :gauss)
        println("You need first to assign a Gauss map to the vertices")
        return false
    end
    for e in edges(g)
        s, t = src(e), dst(e)
        delta = get_prop(g, e, :delta)
        if get_prop(g, e, :dir)=="left"
            hij = get_prop(g, s, :h)*get_prop(g, t, :h)
            lax = vMatrix0(delta, hij)
            set_prop!(g, e, :lax, lax)
        else
            dh = get_prop(g, t, :h)/get_prop(g, s, :h)
            lax = uMatrix0(delta, dh)
            set_prop!(g, e, :lax, lax)
        end
    end
end


function setup_frame!(g::MetaDiGraph; t=0.0)
    phi0 = Complex{Float64}[1 0;0 1]
    set_prop!(g, 1, :frame, phi0)
    set_prop!(g, 1, :dframe, zeros(2,2))
    trav = treeIter(trav_tree(g))
    for tup in trav
        src, dst = tup
        e = Edge(src, dst)
        phi = get_prop(g, src, :frame)
        phit = get_prop(g, src, :dframe)
        if e in edges(g)
            lax = get_lax(g, e, t)
            dlax = get_dlax(g, e, t)
            set_prop!(g, dst, :frame, lax*phi)
            set_prop!(g, dst, :dframe, dlax*phi+lax*phit)
        else
            lax = get_lax(g, reverse(e), t)
            dlax = get_dlax(g, reverse(e), t)
            ilax = inv(lax)
            phi1 = ilax*phi
            set_prop!(g, dst, :frame, phi1)
            set_prop!(g, dst, :dframe, ilax*(phit-dlax*phi1))
        end
    end
end


function symBobenko(g::MetaDiGraph)
    for v in vertices(g)
        phi = get_prop(g, v, :frame)
        phit = get_prop(g, v, :dframe)
        surf = matToQ(2*inv(phi)*phit)
        set_prop!(g, v, :surface, surf)
    end
end

###### Initial conditions for Amsler
function generate_Amsler(
    n1::Quaternion{T},
    n2::Quaternion{S},
    d1::Int64,
    d2::Int64,
) where {S,T<:Real}
    n1 = normalize(n1)
    n2 = normalize(n2)
    e = cross(n1, n2)
    e1 = cross(n1, e)
    e2 = cross(n2, e)
    great1 = [cos(2 * pi / d1 * (i - 1)) * e + sin(2 * pi / d1 * (i - 1)) * e1 for i = 1:d1]
    great2 = [cos(2 * pi / d2 * (i - 1)) * e + sin(2 * pi / d2 * (i - 1)) * e2 for i = 1:d2]
    great1[end] = great1[1]
    great2[end] = great2[1]
    return great1, great2
end

spherical_dist(q1::Quaternion, q2::Quaternion) = acos(dot(q1, q2))

function spherical_angle(
    n::Quaternion, n1::Quaternion, n2::Quaternion)::Float64
    a = normalize(cross(n1,n))
    b = normalize(cross(n2,n))
    o = orientation(n, a, b)
    if n.k==1
        return o ? -acos(dot(a,b)) : acos(dot(a,b))
    else
        return o ? acos(dot(a,b)) : -acos(dot(a,b))
    end
end


stereo(x::Quaternion) = x.k!=1 ? 1/(1-x.k)*(x.i+im*x.j) : 1/(1+x.k)*(x.i+im*x.j)

function dstereoN(p::Quaternion, x::Quaternion)
    dF = [1 0 p.i/(1-p.k)^2;0 1 p.j/(1-p.k)^2]
    X = qToR3(x)
    return dF*X
end

function dstereoS(p::Quaternion, x::Quaternion)
    dF = [1 0 -p.i/(1+p.k)^2;0 1 -p.j/(1+p.k)^2]
    X = qToR3(x)
    return dF*X
end


function orientation(n::Quaternion, n1::Quaternion, n2::Quaternion)::Bool
    if n.k!=1
        x1 = dstereoN(n, n1)
        x2 = dstereoN(n, n2)
        x = [[x1, x2][i][j] for i=1:2, j=1:2]
        return sign(det(x))==1
    else
        x1 = dstereoS(n, n1)
        x2 = dstereoS(n, n2)
        x = [[x1, x2][i][j] for i=1:2, j=1:2]
        return sign(det(x))==1
    end
end
