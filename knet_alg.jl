using Makie: mesh, wireframe, wireframe!

include("ddg.jl")
using .Ddg, .Quaternions

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

function nextN(n::Quaternion, n1::Quaternion, n2::Quaternion)::Quaternion
    dot(n, n1 + n2) / (1 + dot(n1, n2)) * (n1 + n2) - n
end

function nextF(
    n::Quaternion{T},
    n1::Quaternion{S},
    f::Quaternion{U}
) where {S,T,U<:Real}
    return f + cross(n1, n)
end

spherical_dist(q1::Quaternion, q2::Quaternion) = acos(dot(q1, q2))

function get_lengths(
    n::Quaternion, n1::Quaternion, n2::Quaternion)::Tuple{Float64, Float64}
    return spherical_dist(n,n1), spherical_dist(n,n2)
end

function get_angles(
    n::Quaternion, n1::Quaternion, n2::Quaternion)::Float64
    edge1 = normalize(n1-n)
    edge2 = normalize(n2-n)
    asin(abs(cross(edge1, edge2)))
end

function get_k(
    n::Quaternion, n1::Quaternion, n2::Quaternion)::Float64
    s1, s2 = spherical_dist(n,n1), spherical_dist(n,n2)
    return tan(s1/2)*tan(s2/2)
end

function dataFromGauss(gauss::Array{Quaternion, 2}, func::Function)
    m, n = size(gauss)
    output = Array{Any, 2}(undef, m, n-1)
    for i=1:m
        for j=1:n-1
            output[i, j] = func(gauss[i,j], gauss[nw(i,j,m)...],gauss[no(i,j,m)...])
        end
    end
    output
end

function geoKnet(gauss::Array{Quaternion, 2}; f0=zero(Quaternion)::Quaternion,
    zigzag=true)
    m, n = size(gauss)
    surf = zeros(Quaternion, (m,n))
    surf[1, 1] = f0
    if zigzag
        for j=1:2:n
            k, l = 1, j
            for i=1:2*m-1
                if i%2==1
                    g = gauss[k, l]
                    g1 = gauss[no(k,l,m)...]
                    surf[no(k,l,m)...] = nextF(g, g1, surf[k,l])
                    k, l  = no(k,l,m)
                elseif i%2==0
                    g = gauss[k, l]
                    g1 = gauss[so(k,l,m)...]
                    surf[so(k,l,m)...] = nextF(g1, g, surf[k,l])
                    k, l  = so(k,l,m)
                end
            end
            if j+1<n
                g = gauss[1, j+1]
                g2 = gauss[1, j+2]
                surf[1, j+2] = nextF(g2, g, surf[1, j+1])
            end
        end
        if n%2==1
            for i=1:m
                k, l= se(i,n,m)
                g = gauss[k,l]
                g1 = gauss[i,n]
                surf[i,end] = nextF(g1,g,surf[k,l])
            end
        end
        return surf
    else
        for i = 1:m-1
            surf[i+1, 1] = nextF(gauss[i, 1], gauss[i+1, 1], surf[i, 1])
        end
        for i = 1:m, j = 1:n-1
            surf[i, j+1] = nextF(
                gauss[i, j],
                gauss[i, j+1],
                surf[i, j];
                first = false,
            )
        end
        return surf
    end
end

function uMatrix(delta, dh)
    t -> [cot(delta/2)*exp(im*dh) im*exp(t);im*exp(t) cot(delta/2)*exp(-im*dh)]
end

function vMatrix(delta, hij)
    t -> [1 im*exp(-t)*tan(delta/2)*exp(im*hij);
          im*exp(-t)*tan(delta/2)*exp(-im*hij) 1]
end

diff_u(t) = [0 im*exp(t);im*exp(t) 0]

function diff_v(delta, hij)
    t -> [0 -im*exp(-t)*tan(delta/2)*exp(im*hij);
          -im*exp(-t)*tan(delta/2)*exp(-im*hij) 0]
end



#=

m = 24
n = 30
q = Quaternion([0,0,1])
c1 = sample_small_circle(q, 0.3, m)
c2 = sample_small_circle(q, 0.4, m; shift=true)
gauss = propagate_zigzag(c1, c2, nextN, n)
verts = verts_zigzag(gauss)
conn = faces_zigzag(m,n+2)

#delta = dataFromGauss(gauss, get_lengths)
#psi = dataFromGauss(gauss, get_angles)
#k = dataFromGauss(gauss, get_k)

surf = geoKnet(gauss)
sverts = verts_zigzag(surf)

color = zeros(m*(n+2))
color[10]=-0.2
scene = mesh(sverts, conn, color=color, shading=false)
wireframe!(scene[end][1], color = (:black, 0.6), linewidth = 3)

show(surf);
sx, h = textslider(-1:0.01:1, "height", start=0.0)
q = Quaternion([0,1,1])
cw = lift(x -> to_coords.(sample_small_circle(q, x, 12)), h)
sc = lines(cw, showaxis = true, limits = HyperRectangle(Vec3f0(-3), Vec3f0(3)))
final = hbox(sx, sc)

=#
