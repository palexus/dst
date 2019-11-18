include("combinatoric.jl")
include("knet_alg.jl")


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

function build_gauss(
    verlst::Array{Quaternion{T},1},
    horlst::Array{Quaternion{S},1},
) where {T<:Real,S<:Real}
    propagate(verlst, horlst, nextN)
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
                    elseif get_prop(g, v, w, :dir)=="right" && haskey(props(g, w), :surface)
                        f = get_prop(g, v, :surface)
                        n = get_prop(g, v, :gauss)
                        n1 = get_prop(g, w, :gauss)
                        newf1 = nextF(n1, n, f)
                        oldf1 = get_prop(g, w, :surface)
                        @show abs(newf1-oldf1)
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


function testknet(g::MetaDiGraph)
    for e in edges(g)
        println("---------It starts_-------")
        v,w = src(e),dst(e)
        f = get_prop(g, v, :surface)
        f1 = get_prop(g, w, :surface)
        n = get_prop(g, v, :gauss)
        n1 = get_prop(g, w, :gauss)
        if get_prop(g, e, :dir)=="left"
            println("__________left___________")
            @show e
            println("Der Abstand ist ",abs(f1-f-cross(n, n1))<10^-7)
            if abs(f1-f-cross(n,n1))>10^-7
                @show abs(f1-f-cross(n,n1))
                @show f1-f
                @show cross(n,n1)
            end
        else
            println("____________right_________")
            @show e
            println("Der Abstand ist ",abs(f1-f-cross(n1, n))<10^-7)
            if abs(f1-f-cross(n1,n))>10^-7
                @show abs(f1-f-cross(n1,n))
                @show f1-f
                @show cross(n1,n)
            end
        end
    end
end

function myplot!(g::MetaDiGraph, color)
    surf = get_vprops(g, :surface)
    verts = qToR3.(surf)
    verts = [verts[i][j] for i=1:length(verts), j=1:3]
    conn = get_triangles(g)
    conn = get_quads(g)
    conn = [conn[i][j] for i=1:length(conn), j=1:4]

    scene = mesh(verts, conn, color=color, shading=false)
    wireframe!(scene[end][1], color = (:black, 0.6), linewidth = 3)
    display(scene)
end


####### Initial conditions for zigzag
m = 10
n = 20
q = Quaternion([0,0,1])
c1 = sample_small_circle(q, 0.3, m)
c2 = sample_small_circle(q, 0.4, m; shift=true)
gauss = propagate_zigzag(c1, c2, nextN, n)


####### plot knet
g = zigzag(m,n+2)
set_vprops!(g, gauss, :gauss)
knet!(g)
color = zeros(m*(n+2))
color[1]=-0.2
myplot!(g, color)


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

m, n = 30, 30
great1, great2 = generate_Amsler(Quaternion([1, 0, 0]),
                                 Quaternion([0, -1, 0]), m, n)
gauss = build_gauss(great1, great2)


###### Plot Amsler
g = di_grid(m,n)
set_vprops!(g, gauss, :gauss)
knet!(g)

color = zeros(m*n)
color[1]=-0.2
myplot!(g, color)
