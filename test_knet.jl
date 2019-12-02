include("combinatoric.jl")
include("knet.jl")
using Makie: textslider, lift, mesh!, wireframe!, hbox, vbox, campixel!, slider
using Makie: Scene, scatter!
using Test
using LinearAlgebra: det


function test_knet!(g::MetaDiGraph)
    passed = true
    for e in edges(g)
        v,w = src(e),dst(e)
        f = get_prop(g, v, :surface)
        f1 = get_prop(g, w, :surface)
        n = get_prop(g, v, :gauss)
        n1 = get_prop(g, w, :gauss)
        phi = get_prop(g, v, :frame)
        phi1 = get_prop(g, w, :frame)
        n = matToQ(-im*inv(phi)*[1 0;0 -1]*phi)
        n1 = matToQ(-im*inv(phi1)*[1 0;0 -1]*phi1)
        if get_prop(g, e, :dir)=="left"
            try
                diff = imag(f1-f)-cross(n,n1)
                @test zero(Quaternion) ≈ diff atol=10^-7
            catch TestSetException
                passed = false
                println("Test failed for the edge $e")
                @show f1, f
                @show n1, n
            end
        else
            try
                diff = imag(f1-f)-cross(n1,n)
                @test zero(Quaternion) ≈ diff atol=10^-7
            catch TestSetException
                passed = false
                println("Test failed for the edge $e")
                @show f1, f
                @show n1, n
            end
        end
    end
    return passed
end


function test_setup_lax!(g::MetaDiGraph)
    passed = true
    for v in vertices(g)
        out = outneighbors(g, v)
        v12 = opposite_vertex(g, v)
        if v12!=nothing && left(g, v)!=nothing && right(g,v)!=nothing
            U = get_prop(g, v, out[1], :lax)
            V1 = get_prop(g, out[1], v12, :lax)
            V = get_prop(g, v, out[2], :lax)
            U2 = get_prop(g, out[2], v12, :lax)
            try
                @test V1*U-U2*V ≈ zeros(2,2) atol=10^-5
            catch TestSetException
                passed = false
                @show v, out[1], out[2], v12
            end
        end
    end
    return passed
end


function test_setup_h!(g::MetaDiGraph)
    passed = true
    for v in vertices(g)
        v2=left(g,v)
        v1=right(g,v)
        v12=opposite_vertex(g,v)
        if v1!=nothing && v2!=nothing && v12!=nothing
            deltau=get_prop(g, v, v1, :delta)
            deltav=get_prop(g, v, v2, :delta)
            k = tan(deltau/2)*tan(deltav/2)
            H = get_prop(g, v, :h)
            H1 = get_prop(g, v1, :h)
            H2 = get_prop(g, v2, :h)
            H12 = get_prop(g, v12, :h)
            try
                @test H12*H ≈ (k+H1*H2)/(1+k*H1*H2) atol=10^-7
            catch e
                if typeof(e)<:TestSetException
                    passed = false
                    println("Test failed at the face")
                    @show v, v1, v12, v2
                    @show h, h1, h12, h2, k
                elseif typeof(e)<:UndefVarError
                    @show v, v1, v12, v2
                else
                    rethrow(e)
                end
            end
        end
    end
    return passed
end


function test_setup_frame!(g::MetaDiGraph)
    passed = true
    m, n = get_prop(g, :dim)
    for v=1:m
        next = mod(v, m)+1
        phi = get_prop(g, v, :frame)
        phi1 = get_prop(g, next, :frame)
        laxr = get_prop(g, v, right(g,v), :lax)
        laxl = get_prop(g, next, left(g,next), :lax)
        actphi = get_prop(g, right(g,v), :frame)
        # Normalize for the test (overlap does not work of course without norm)
        phi = 1/sqrt(det(phi))*phi
        phi1 = 1/sqrt(det(phi1))*phi1
        laxl = 1/sqrt(det(laxl))*laxl
        laxr = 1/sqrt(det(laxr))*laxr
        actphi = 1/sqrt(det(actphi))*actphi
        try
            @test abs.(laxr*phi) ≈ abs.(laxl*phi1) atol=10^-7
            @test abs.(laxr*phi) ≈ abs.(actphi) atol=10^-7
            @test abs.(laxl*phi1) ≈ abs.(actphi) atol=10^-7
        catch TestSetException
            passed = false
            println("Test failed at the pair")
            @show v, next
            @show laxr, phi
            @show laxl, phi1
            @show get_prop(g, right(g,v), :frame)
        end
    end
    return passed
end

function plotverts(g::MetaDiGraph, t::Real)
    setup_frame!(g, t=t)
    symBobenko(g)
    surf = get_vprops(g, :surface)
    verts = qToR3.(surf)
    verts = [verts[i][j] for i=1:length(verts), j=1:3]
    return verts
end

function myplot!(g::MetaDiGraph)
    conn = get_triangles(g)
    conn = [conn[i][j] for i=1:length(conn), j=1:3]
    color = zeros(nv(g))
    color[1]=-0.2

    s1 = slider(LinRange(-1.0, 1.0, 1001),
          raw = true, camera = campixel!, start = 0.0)
    kx = s1[end][:value]

    scene = Scene()
    verts =  lift(x->plotverts(g, x), kx)
    mesh!(scene, verts,
          conn, color=color, shading=false, transparency=true)
    wireframe!(scene[end][1], color = (:black, 0.6), linewidth = 3)
    scatter!(scene, verts, markersize = .03)

    hbox(scene, s1, parent = Scene(resolution = (800, 600)))
end

myplot!(g)

function my_plot!(g::MetaDiGraph)
    conn = get_triangles(g)
    conn = [conn[i][j] for i=1:length(conn), j=1:3]
    color = zeros(nv(g))
    color[1]=-0.2
    surf = get_vprops(g, :surface)
    verts = qToR3.(surf)
    verts = [verts[i][j] for i=1:length(verts), j=1:3]
    initialgauss = get_vprops(g, :gauss)
    nverts = qToR3.(initialgauss)
    nverts = [nverts[i][j] for i=1:length(nverts), j=1:3]
    actualgauss = []
    for v in vertices(g)
        phi = get_prop(g, v, :frame)
        push!(actualgauss, -im*inv(phi)*[1 0;0 -1]*phi)
    end
    sverts = matToR3.(actualgauss)
    sverts = [sverts[i][j] for i=1:length(sverts), j=1:3]
    geog = copy(g)
    knet!(geog)
    geoverts = qToR3.(get_vprops(geog, :surface))
    geoverts = [geoverts[i][j] for i=1:length(geoverts), j=1:3]

    scene1 = mesh(verts, conn, color=color, shading=false)
    wireframe!(scene1[end][1], color = (:black, 0.6), linewidth = 3)
    scene4 = mesh(geoverts, conn, color=color, shading=false)
    wireframe!(scene4[end][1], color = (:black, 0.6), linewidth = 3)
    scene2 = mesh(nverts, conn, color=color, shading=false)
    wireframe!(scene2[end][1], color = (:black, 0.6), linewidth = 3)
    scene3 = mesh(sverts, conn, color=color, shading=false)
    wireframe!(scene3[end][1], color = (:black, 0.6), linewidth = 3)
    display(hbox(vbox(scene4, scene1), vbox(scene2, scene3)))
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


function initial_condition_zigzag(m::Int, n::Int, periodic::Bool)
    q = Quaternion([0,0,1])
    if periodic
        c1 = sample_small_circle(q, 0.3, m)
        c2 = sample_small_circle(q, 0.1, m; shift=true)
        gauss = propagate_zigzag(c1, c2, nextN, n-2)
    else
        c1 = sample_small_circle(q, 0.3, m-1)
        c2 = sample_small_circle(q, 0.1, m-1; shift=true)
        gauss = propagate_zigzag(c1, c2, nextN, n-2)
        gauss = vcat(gauss, reshape(gauss[1, :], (1, n)))
    end
    gauss
end


####### plot knet
m = 20
n = 40
periodic = false
g = zigzag(m,n, periodic=periodic)

# Initial Data
gauss = initial_condition_zigzag(m, n, periodic)
set_vprops!(g, gauss, :gauss)

# setup everything
setup_h!(g)
setup_lax!(g)
setup_frame!(g, t=0.0)
symBobenko(g)

myplot!(g)


# Here are the tests
println("---------------")
@test test_setup_h!(g)
println("---------------")
@test test_setup_lax!(g)
println("---------------")
@test test_knet!(g)
println("---------------")

# Alg surface plot and plot geo/alg Gauss/surface plot
myplot!(g)
my_plot!(g)





###### Plot Amsler
m, n = 20, 20
great1, great2 = generate_Amsler(Quaternion([1, 0, 0]),
                                 Quaternion([0, 1, 0]), m, n)
gauss = build_gauss(great1, great2)[1:8,1:8]

gr = di_grid(8,8)

set_vprops!(gr, gauss, :gauss)

setup_h!(gr)

setup_lax!(gr)
setup_frame!(gr)
symBobenko(gr)
myplot!(gr)


my_plot!(gr)

gauss = build_gauss(great1, great2)
geogr=di_grid(m,n)
set_vprops!(geogr, gauss, :gauss)
knet!(geogr)
myplot!(geogr)
