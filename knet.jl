include("combinatoric.jl")
include("quaternion.jl")

using .Quaternions

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
    for v in vertices(g)
        try
            out = outneighbors(g, v)
            if haskey(props(g, v), :surface)
                for w in out
                    if get_prop(g, v, w, :dir)=="right" && !haskey(props(g, w), :surface)
                        f = get_prop(g, v, :surface)
                        n = get_prop(g, v, :gauss)
                        n1 = get_prop(g, w, :gauss)
                        f1 = nextF(n1, n, f)
                        set_prop!(g, w, :surface, f1)
                    elseif get_prop(g, v, w, :dir)=="left" && !haskey(props(g, w), :surface)
                        f = get_prop(g, v, :surface)
                        n = get_prop(g, v, :gauss)
                        n1 = get_prop(g, w, :gauss)
                        f1 = nextF(n, n1, f)
                        set_prop!(g, w, :surface, f1)
                    end
                end
            else
                w = out[findfirst(x->haskey(props(g, x), :surface), out)]
                f2 = get_prop(g, w, :surface)
                n2 = get_prop(g, w, :gauss)
                n = get_prop(g, v, :gauss)
                f = get_prop(g, v, w, :dir)=="left" ? nextF(n2, n, f2) : nextF(n,n2,f2)
                set_prop!(g, v, :surface, f)
                id = findfirst(x->!haskey(props(g, x), :surface), out)
                if id!=nothing
                    n1 = get_prop(g, out[id], :gauss)
                    f1 = nextF(n1, n, f)
                    set_prop!(g, out[id], :surface, f1)
                end
            end
        catch KeyError
            println("There is probably no gauss map")
        end
    end
end


function setup_h!(g::MetaDiGraph)
    m, n = get_prop(g, :dim)
    set_prop!(g, 1, :h, 0)
    set_prop!(g, 2*m, :h, 0)
    #fill lowest layer first with the inner angles
    if get_prop(g, :type)=="zigzag"
        S = 0.5*sum([(-1)^k*get_prop(g, m+k, :iangle) for k=1:m])
        T = 0.5*sum([(-1)^k*get_prop(g, k, :oangle) for k=1:m])
        set_prop!(g, 1, :h, S)
        set_prop!(g, 2m, :h, T)
        for v=2:m
            out = sorted_outneighbors(g, v)
            phi = mod(get_prop(g, out[1], :iangle)-pi, 2pi)
            h = get_prop(g, v-1, :h)
            h1 = -h-phi+pi
            set_prop!(g, v, :h, h1)
        end
    end
    for v in vertices(g)
        out = sorted_outneighbors(g,v)
        if length(out)>0
            deltau = get_prop(g, v, out[1], :delta)
            deltav = get_prop(g, v, out[2], :delta)
            k = tan(deltau/2)*tan(deltav/2)
            id = findfirst(x->haskey(props(g,x),:h), out)
            oid = findfirst(x->!haskey(props(g,x),:h), out)
            phi = mod(get_prop(g, v, :oangle)-pi, 2pi)
            if oid!=nothing
                h2 = get_prop(g, out[id], :h)
                set_prop!(g, out[oid], :h, -phi-h2+pi)
            end
            v12 = opposite_vertex(g, v)
            if v12!=nothing
                h = get_prop(g, v, :h)
                h1 = get_prop(g, out[1], :h)
                h2 = get_prop(g, out[2], :h)
                #phiR = angle((1-exp(im*phi)*k)/(k-exp(im*phi)))
                phiR = mod(pi+2*angle(1-k*exp(im*phi)),2pi)-phi
                h12 = phiR-h
                for i=-2:2
                    if sineGordon(h, h1, h12+i*pi, h2, k)
                        set_prop!(g, v12, :h, h12+i*pi)
                        break
                    end
                end
            end
        end
    end
end


function sineGordon(h, h1, h12, h2, k)
    diff = sin(0.5*(h1+h2-h-h12))-k*sin(0.5*(h+h1+h2+h12))
    return diff < 10^-7
end


function uMatrix0(delta, dh)
    [cot(delta/2)*exp(im*dh) im;im cot(delta/2)*exp(-im*dh)]
end


function vMatrix0(delta, hij)
    [1 im*tan(delta/2)*exp(im*hij);
    im*tan(delta/2)*exp(-im*hij) 1]
end


function setup_lax(g::MetaDiGraph)
    if !haskey(props(g,1), :gauss)
        println("You need first to assign a Gauss map to the vertices")
        return false
    end
    m, n = props(g)[:dim]
    set_eprops!(g, spherical_dist, :gauss, :delta)
    set_vprops_out!(g, get_angles, :gauss, :oangle)
    set_vprops_inn!(g, get_angles, :gauss, :iangle)
    setup_h!(g)
    for e in edges(g)
        s, t = src(e), dst(e)
        delta = get_prop(g, e, :delta)
        if get_prop(g, e, :dir)=="left"
            hij = get_prop(g, s, :h) + get_prop(g, t, :h)
            lax = vMatrix0(delta, hij)
            set_prop!(g, e, :lax, lax)
        else
            dh = get_prop(g, t, :h) - get_prop(g, s, :h)
            lax = uMatrix0(delta, dh)
            set_prop!(g, e, :lax, lax)
        end
    end
end


function setup_frame!(g::MetaDiGraph; t=0.0)
    if haskey(props(g, 1), :frame)
        rem_vprop!(g, :frame)
        rem_vprop!(g, :dframe)
    end
    phi0 = Complex{Float64}[1 0;0 1]
    n = nv(g)
    set_prop!(g, 1, :frame, phi0)
    set_prop!(g, 1, :dframe, zeros(2,2))
    all = collect(2:n)
    while !isempty(all)
        for v in vertices(g)
            out = outneighbors(g, v)
            for o in out
                if o in all
                    !haskey(props(g, o), :frame) || continue
                    haskey(props(g, v), :frame) || break
                    phi = get_prop(g, v, :frame)
                    phit = get_prop(g, v, :dframe)
                    lax = get_lax(g, v, o, t)
                    dlax = get_dlax(g, v, o, t)
                    set_prop!(g, o, :frame, lax*phi)
                    set_prop!(g, o, :dframe, dlax*phi+lax*phit)
                    filter!(el->el≠o, all)
                end
            end
            inn = inneighbors(g, v)
            for i in inn
                if i in all
                    !haskey(props(g, i), :frame) || continue
                    haskey(props(g, v), :frame) || break
                    phi = get_prop(g, v, :frame)
                    phit = get_prop(g, v, :dframe)
                    lax = get_lax(g, i, v, t)
                    ilax = inv(lax)
                    dlax = get_dlax(g, i, v, t)
                    set_prop!(g, i, :frame, ilax*phi)
                    set_prop!(g, i, :dframe, ilax*(phit-dlax*ilax*phi))
                    filter!(el->el≠i, all)
                end
            end
        end
    end
end


function symBobenko(g::MetaDiGraph)
    for v in vertices(g)
        phi = get_prop(g, v, :frame)
        phit = get_prop(g, v, :dframe)
        set_prop!(g, v, :surface, 2(inv(phi)*phit))
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

function get_angles(
    n::Quaternion, n1::Quaternion, n2::Quaternion)::Float64
    edge1 = normalize(n1-n)
    edge2 = normalize(n2-n)
    asin(abs(cross(edge1, edge2)))
end
