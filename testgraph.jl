include("combinatoric.jl")
include("knet_alg.jl")
using Makie: textslider, lift


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


function setup_h!(g::MetaDiGraph)
    m, n = get_prop(g, :dim)
    set_prop!(g, 1, :h, 0)
    set_prop!(g, 2*m, :h, 0)
    #fill lowest layer first with the inner angles
    if get_prop(g, :type)=="zigzag"
        S = mod(0.5*sum([(-1)^k*get_prop(g, m+k, :iangle) for k=1:m]), 2pi)
        T = mod(0.5*sum([(-1)^k*get_prop(g, k, :oangle) for k=1:m]), 2pi)
        set_prop!(g, 1, :h, S)
        set_prop!(g, 2m, :h, T)
        for v=2:m
            out = sorted_outneighbors(g, v)
            psi = get_prop(g, out[1], :iangle)
            h = get_prop(g, v-1, :h)
            h1 = mod(-h-psi, 2pi)
            set_prop!(g, v, :h, h1)
        end
    end
    for v in vertices(g)
        out = sorted_outneighbors(g,v)
        if length(out)>0
            k = 1
            if length(out)==2
                deltau = get_prop(g, v, out[1], :delta)
                deltav = get_prop(g, v, out[2], :delta)
                k = tan(deltau/2)*tan(deltav/2)
            end
            id = findfirst(x->haskey(props(g,x),:h), out)
            oid = findfirst(x->!haskey(props(g,x),:h), out)
            psi = get_prop(g, v, :oangle)
            if oid!=nothing
                h2 = get_prop(g, out[id], :h)
                set_prop!(g, out[oid], :h, -psi-h2)
            end
            v12 = opposite_vertex(g, v)
            if v12!=nothing
                h = get_prop(g, v, :h)
                psiR = angle(-(exp(im*psi)*k+1)/(exp(im*psi)+k))
                set_prop!(g, v12, :h, psiR-h-pi)
            end
        end
    end
end


function sineGordon(h, h1, h12, h2, k)
    diff = sin(0.5*(h1+h2-h-h12))-k*sin(0.5*(h+h1+h2+h12))
    return diff < 10^-7
end


function test_setup_h!(g::MetaDiGraph)
    faces = get_quads(g)
    test = []
    for f in faces
        i, i1, i2, i12 = f
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
        println(e)
        @show delta
        if get_prop(g, e, :dir)=="left"
            hij = get_prop(g, s, :h) + get_prop(g, t, :h)
            lax = vMatrix0(delta, hij)
            println("left")
            @show hij, lax
            set_prop!(g, e, :lax, lax)
        else
            dh = get_prop(g, t, :h) - get_prop(g, s, :h)
            lax = uMatrix0(delta, dh)
            println("right")
            @show dh, lax
            set_prop!(g, e, :lax, lax)
        end
    end
end


function setup_frame!(g::MetaDiGraph; t=0.0)
    phi0 = Complex{Float64}[1 0;0 1]
    n = nv(g)
    set_prop!(g, 1, :frame, phi0)
    set_prop!(g, 1, :dframe, zeros(2,2))
    all = collect(2:n)
    while !isempty(all)
        @show all
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
                    @show o, all
                    filter!(el->el≠o, all)
                    @show o, all
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

function initial_condition_zigzag(m::Int, n::Int)
    q = Quaternion([0,0,1])
    c1 = sample_small_circle(q, 0.3, m)
    c2 = sample_small_circle(q, 0.4, m; shift=true)
    gauss = propagate_zigzag(c1, c2, nextN, n-2)
    gauss
end

####### plot knet
m = 8
n = 5
g = zigzag(m,n)
gauss = initial_condition_zigzag(m, n)
set_vprops!(g, gauss, :gauss)

setup_lax(g)
@show test_setup_h!(g)

setup_frame!(g)
symBobenko(g)

myplot!(g)

print_eprop(g, :delta)
print_vprop(g, :h)
print_vprop(g, :frame)
print_eprop(g, :lax)


function myplot!(g::MetaDiGraph)
    conn = get_triangles(g)
    conn = [conn[i][j] for i=1:length(conn), j=1:3]
    color = zeros(nv(g))
    color[1]=-0.2

    #sx, h = textslider(-1:0.01:1, "height", start=0.0)
    #cw = lift(x -> setup_frame!(g, t=x), h)

    #symBobenko(g)
    surf = get_vprops(g, :surface)
    verts = matToR3.(surf)
    verts = [verts[i][j] for i=1:length(verts), j=1:3]

    scene = mesh(verts, conn, color=color, shading=false)
    wireframe!(scene[end][1], color = (:black, 0.6), linewidth = 3)

    #sc = lines(cw, showaxis = true, limits = HyperRectangle(Vec3f0(-3), Vec3f0(3)))
    #final = hbox(sx, sc)
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
setup_lax(g)


color = zeros(m*n)
color[1]=-0.2
myplot!(g, color)
