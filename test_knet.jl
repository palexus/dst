include("combinatoric.jl")
include("knet.jl")
using Makie: textslider, lift, mesh, wireframe!, hbox
using Test

function test_knet!(g::MetaDiGraph)
    for e in edges(g)
        v,w = src(e),dst(e)
        f = get_prop(g, v, :surface)
        f1 = get_prop(g, w, :surface)
        n = get_prop(g, v, :gauss)
        n1 = get_prop(g, w, :gauss)
        if get_prop(g, e, :dir)=="left"
            diff = f1-f-cross(n,n1)
            @test zero(Quaternion) ≈ diff atol=10^-7
        else
            diff = f1-f-cross(n1,n)
            @test zero(Quaternion) ≈ diff atol=10^-7
        end
    end
end


function test_setup_h!(g::MetaDiGraph)
    faces = get_quads(g)
    test = []
    for f in faces
        i, i1, i2, i12 = sort(f)
        h = get_prop(g, i, :h)
        h1 = get_prop(g, i1, :h)
        h2 = get_prop(g, i2, :h)
        h12 = get_prop(g, i12, :h)
        deltau = get_prop(g, i, i1, :delta)
        deltav = get_prop(g, i, i2, :delta)
        k = tan(deltau/2)*tan(deltav/2)
        push!(test, sineGordon(h, h1, h12, h2, k))
    end
    return test
end

function myplot!(g::MetaDiGraph)
    conn = get_triangles(g)
    conn = [conn[i][j] for i=1:length(conn), j=1:3]
    color = zeros(nv(g))
    color[1]=-0.2
    surf = get_vprops(g, :surface)
    verts = []
    if typeof(get_prop(g, 1, :surface))==Matrix{Complex{Float64}}
        verts = matToR3.(surf)
    else
        verts = qToR3.(surf)
    end
    verts = [verts[i][j] for i=1:length(verts), j=1:3]
    scene = mesh(verts, conn, color=color, shading=false)
    wireframe!(scene[end][1], color = (:black, 0.6), linewidth = 3)
    display(scene)
end

function plot_gauss(g::MetaDiGraph)
    conn = get_triangles(g)
    conn = [conn[i][j] for i=1:length(conn), j=1:3]
    color = zeros(nv(g))
    color[1]=-0.2
    gauss = get_vprops(g, :gauss)
    verts = qToR3.(gauss)
    verts = [verts[i][j] for i=1:length(verts), j=1:3]
    scene = mesh(verts, conn, color=color, shading=false)
    wireframe!(scene[end][1], color = (:black, 0.6), linewidth = 3)
    display(scene)
end


function initial_condition_zigzag(m::Int, n::Int)
    q = Quaternion([0,0,1])
    c1 = sample_small_circle(q, 0.3, m)
    c2 = sample_small_circle(q, 0.4, m; shift=true)
    gauss = propagate_zigzag(c1, c2, nextN, n-2)
    gauss
end

####### plot knet
m = 7
n = 3
g = zigzag(m,n; periodic=true)
gauss = initial_condition_zigzag(m, n)
set_vprops!(g, gauss, :gauss)
knet!(g)
test_knet!(g)



setup_lax(g)
@show test_setup_h!(g)

setup_frame!(g)
symBobenko(g)

myplot!(g)

rem_vprop!(g, :surface)
knet!(g)
myplot!(g)

plot_gauss(g)


###### Plot Amsler
m, n = 30, 30
great1, great2 = generate_Amsler(Quaternion([1, 0, 0]),
                                 Quaternion([0, -1, 0]), m, n)
gauss = build_gauss(great1, great2)

gr = di_grid(m,n)
set_vprops!(g, gauss, :gauss)

knet!(gr)
myplot!(gr)
