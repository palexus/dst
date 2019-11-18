include("combinatoric.jl")
include("knet_alg.jl")


m = 5
n = 4
q = Quaternion([0,0,1])
c1 = sample_small_circle(q, 0.3, m)
c2 = sample_small_circle(q, 0.4, m; shift=false)
gauss = propagate_zigzag(c1, c2, nextN, n)
verts = verts_zigzag(gauss)
conn = faces_zigzag(m,n+2)


g = zigzag(m,n)

set_vprops!(g, gauss, :gauss)
print_vprop(g, :gauss)

delta = dataFromGauss(gauss, get_lengths)
set_eprops!(g, delta, :delta)
dic = props(g)

print_eprop(g, :delta)

set_vprops!(g, get_angles, :gauss, :angle)

print_vprop(g, :angle)
print_vprops(g)
props(g)

function nextF(
    n::Quaternion{T},
    n1::Quaternion{S},
    f::Quaternion{U}
) where {S,T,U<:Real}
    return f + cross(n1, n)
end

function knet(g::MetaDiGraph; f0=zero(Quaternion))
    try
        set_prop!(g, 1, :surface, f0)
        for v in vertices(g)
            out = outneighbors(g, v)
            if haskey(props(g, v), :surface)
                if length(out)==2
                    w1, w2 = out
                    if get_prop(g, v, w1, :dir)=="left"
                        w1, w2 = w2, w1
                    end
                    f = get_prop(g, v, :surface)
                    n = get_prop(g, v, :gauss)
                    if !haskey(props(g, w1), :surface)
                        n1 = get_prop(g, w1, :gauss)
                        f1 = nextF(n, n1, f)
                        set_prop!(g, w1, :surface, f1)
                    end
                    if !haskey(props(g, w2), :surface)
                        n2 = get_prop(g, w2, :gauss)
                        f2 = nextF(n2, n, f)
                        set_prop!(g, w2, :surface, f2)
                    end
                end
            end
        end
    catch KeyError
        println("There is probably no gauss map")
    end
end

gr = di_grid([21,21])
great1, great2 = generate_Amsler(
    Quaternion([1, 0, 0]),
    Quaternion([0, 1, 0]),
    20, 20)

N = build_gauss(great1, great2)



set_vprops!(gr, N, :gauss)

knet(gr)

print_vprop(gr, :surface)
