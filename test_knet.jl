include("combinatoric.jl")
include("knet.jl")
using Makie: textslider, lift, mesh, wireframe!, hbox, vbox
using Test


function test_knet!(g::MetaDiGraph)
    for e in edges(g)
        v,w = src(e),dst(e)
        f = get_prop(g, v, :surface)
        f1 = get_prop(g, w, :surface)
        n = get_prop(g, v, :gauss)
        n1 = get_prop(g, w, :gauss)
        #phi = get_prop(g, v, :frame)
        #phi1 = get_prop(g, w, :frame)
        #n = matToQ(-im*inv(phi)*[1 0;0 -1]*phi)
        #n1 = matToQ(-im*inv(phi1)*[1 0;0 -1]*phi1)
        if get_prop(g, e, :dir)=="left"
            diff = imag(f1-f)-cross(n,n1)
            @show e
            @test zero(Quaternion) ≈ diff atol=10^-7
        else
            @show e
            diff = imag(f1-f)-cross(n1,n)
            @test zero(Quaternion) ≈ diff atol=10^-7
        end
    end
end

function test_setup_lax!(g::MetaDiGraph)
    for v in vertices(g)
        out = outneighbors(g, v)
        v12 = opposite_vertex(g, v)
        U = get_prop(g, v, out[1], :lax)
        V1 = get_prop(g, out[1], v12, :lax)
        V = get_prop(g, v, out[2], :lax)
        U2 = get_prop(g, out[2], v12, :lax)
        @show v
        @test V1*U-U2*V ≈ zeros(2,2) atol=10^-7
    end
end


function test_setup_h!(g::MetaDiGraph)
    faces = get_quads(g)
    for f in faces
        i, i1, i2, i12 = sort(f)
        h = get_prop(g, i, :h)
        h1 = get_prop(g, i1, :h)
        h2 = get_prop(g, i2, :h)
        h12 = get_prop(g, i12, :h)
        deltau = get_prop(g, i, i1, :delta)
        deltav = get_prop(g, i, i2, :delta)
        k = tan(deltau/2)*tan(deltav/2)
        try
            @test sineGordon(h, h1, h12, h2, k)
        catch TestSetException
            println("Test failed at the face")
            @show f
            @show h, h1, h2, h12, k
            println(sin(0.5*(h1+h2-h-h12))-k*sin(0.5*(h+h1+h2+h12)))
        end
    end
end


function myplot!(g::MetaDiGraph)
    conn = get_triangles(g)
    conn = [conn[i][j] for i=1:length(conn), j=1:3]
    color = zeros(nv(g))
    color[1]=-0.2
    surf = get_vprops(g, :surface)
    verts = qToR3.(surf)
    verts = [verts[i][j] for i=1:length(verts), j=1:3]
    scene = mesh(verts, conn, color=color, shading=false)
    wireframe!(scene[end][1], color = (:black, 0.6), linewidth = 3)
    display(scene)
end

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

    scene1 = mesh(verts, conn, color=color, shading=false)
    wireframe!(scene1[end][1], color = (:black, 0.6), linewidth = 3)
    scene2 = mesh(nverts, conn, color=color, shading=false)
    wireframe!(scene2[end][1], color = (:black, 0.6), linewidth = 3)
    scene3 = mesh(sverts, conn, color=color, shading=false)
    wireframe!(scene3[end][1], color = (:black, 0.6), linewidth = 3)
    display(hbox(scene1, vbox(scene2, scene3)))
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
    c1 = sample_small_circle(q, 0.15, m)
    c2 = sample_small_circle(q, 0.1, m; shift=true)
    gauss = propagate_zigzag(c1, c2, nextN, n-2)
    gauss
end


####### plot knet
m = 11
n = 10
g = zigzag(m,n, periodic=true)
gauss = initial_condition_zigzag(m, n)
set_vprops!(g, gauss, :gauss)
#knet!(g)
#test_knet!(g)
#myplot!(g)
#plot_gauss(g)

setup_lax!(g)
println("---------------------------")
test_setup_h!(g)

setup_frame!(g)

symBobenko(g)
#test_knet!(g)

myplot!(g)
my_plot!(g)


psi1 = get_prop(g, 11, :oangle)
psi2 = get_prop(g, 11, :langle)
psi3 = get_prop(g, 11, :iangle)
psi4 = get_prop(g, 11, :rangle)

psi1+psi2+psi3+psi4











###### Plot Amsler
m, n = 20, 20
great1, great2 = generate_Amsler(Quaternion([1, 0, 0]),
                                 Quaternion([0, 1, 0]), m, n)
gauss = build_gauss(great1, great2)

gr = di_grid(m,n)
set_vprops!(gr, gauss, :gauss)

#knet2!(gr)
setup_lax!(gr)
setup_frame!(gr)
symBobenko(gr)
myplot!(gr)


out = outneighbors(gr, 24)
phiU1 = get_prop(gr, 34, :oangle)
phiU2 = get_prop(gr, 34, :langle)
phiU3 = get_prop(gr, 34, :iangle)
phiU4 = get_prop(gr, 34, :rangle)

phiU1+phiU2+phiU3+phiU4
2pi

phiU = get_prop(gr, 14, :oangle)
phiL = get_prop(gr, 15, :langle)
phiR = get_prop(gr, 34, :rangle)
phiD = get_prop(gr, 35, :iangle)

deltau = get_prop(gr, 14, 15, :delta)
deltav = get_prop(gr, 14, 34, :delta)
k = tan(deltau/2)*tan(deltav/2)

phiR2 = angle(-(k*exp(im*(phiU))+1)/(exp(im*(phiU))+k))
phi = angle(-(k*exp(im*(phiR2))+1)/(exp(im*(phiR2))+k))



geogr=di_grid(m,n)
set_vprops!(geogr, gauss, :gauss)
knet!(geogr)
myplot!(geogr)
