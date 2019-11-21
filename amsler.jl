using Makie
using GeometryTypes
using Colors

#include("quaternion.jl")
#using .Quaternions

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
    f::Quaternion{U};
    first = true,
) where {S,T,U<:Real}
    if first
        return f + cross(n1, n)
    else
        return f + cross(n, n1)
    end
end

function build_Knet(
    gauss::Array{Quaternion,2};
    f0 = zero(Quaternion)::Quaternion,
)
    m, n = size(gauss)
    surf = zeros(Quaternion, m, n)
    surf[1, 1] = f0
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


"""
    generate_Amsler(n1, n2, d1, d2)

Compute the Cauchy data of the Gauss map of the Amsler surface. n1 and n2 define
two great circles and d1 and d2 prescribe the number of equally distributed
points on the great circle.

# Examples
```julia-repl
julia> great1, great2 = generate_Amsler(Quaternion([1,0,0]), Quaternion([0,1,0]), 12, 12)

```
"""
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
    great1 = [cos(2 * pi / d1 * (i - 1)) * e + sin(2 * pi / d1 * (i - 1)) * e1 for i = 1:d1+1]
    great2 = [cos(2 * pi / d2 * (i - 1)) * e + sin(2 * pi / d2 * (i - 1)) * e2 for i = 1:d2+1]
    great1[end] = great1[1]
    great2[end] = great2[1]
    return great1, great2
end


# Evolve in 4 quadrants and glue together
function build_Amsler()
    great1, great2 = generate_Amsler(
        Quaternion([1, 0, 0]),
        Quaternion([0, 1, 0]),
        20, 20)
    N = build_gauss(great1, great2)
    F = build_Knet(N)
    N1 = build_gauss(reverse(great1), great2)
    F1 = build_Knet(N1)
    N2 = build_gauss(reverse(great1), reverse(great2))
    F2 = build_Knet(N2)
    N3 = build_gauss(great1, reverse(great2))
    F3 = build_Knet(N3)
    bigFN = vcat([F1[end:-1:2,:], F]...)
    bigFS = vcat([F2[end:-1:2,:], F3]...)
    surf = hcat([bigFS[:, end:-1:2], bigFN]...)
    x = proj_i.(surf)
    y = proj_j.(surf)
    z = proj_k.(surf)
    return x,y,z
end

function my_plot(x,y,z)
    n = length(x)
    scene =  Makie.surface(x,y,z,show_axis = false,
        color = fill(RGBA(1.,1.,1.,0.4), n, n),
        limits = HyperRectangle(Vec3f0(-5), Vec3f0(5)))
    cam3d_cad!(scene)
    wireframe!(scene, x, y, z, show_axis = false,
               linewidth = 2, color = RGBA(0.0,0.0,0.2,0.4))
    scene.center = false

    #s1 = slider(LinRange(0.01, 1, 100), raw = true, camera = campixel!, start = 0.3)
    #s2 = slider(LinRange(-2pi, 2pi, 100), raw = true, camera = campixel!)
end

x, y, z = build_Amsler()

my_plot(x,y,z)

#Makie.save("plot2.png", scene)



#=
### This is how to plot the great circles
great1, great2 = generate_Amsler(Quaternion([-1,5.9,-4]), Quaternion([1,1,0]), 6, 24)
x = proj_i.(great1)
y = proj_j.(great1)
z = proj_k.(great1)

x1 = proj_i.(great2)
y1 = proj_j.(great2)
z1 = proj_k.(great2)

function my_plot()
    scene = lines(x, y, z, color = :blue)
    lines!(scene, x1, y1, z1, color = :black)
end

my_plot()
=#
