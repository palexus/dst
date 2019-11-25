using LightGraphs, MetaGraphs
import Base: reverse

function nw(i::Int, j::Int, m::Int)
    if j%2 == 0
        return i, j+1
    else
        if i==1
            return m, j+1
        else
            return i-1, j+1
        end
    end
end

function so(i::Int, j::Int, m::Int)
    if j%2 == 1
        return i, j-1
    else
        if i==m
            return 1, j-1
        else
            return i+1, j-1
        end
    end
end

function no(i::Int, j::Int, m::Int)
    if j%2 == 1
        return i, j+1
    else
        if i==m
            return 1, j+1
        else
            return i+1, j+1
        end
    end
end

function sw(i::Int, j::Int, m::Int)
    if j%2 == 0
        return i, j-1
    else
        if i==1
            return m, j-1
        else
            return i-1, j-1
        end
    end
end

function toIndex(i::Int, j::Int, row::Int)::Int
    i+(j-1)*row
end

function indexTo(index::Int, row::Int)
    i = mod(index, row)
    j = convert(Int, floor(index/row)+1)
    return i, j
end


"""
    di_grid(dims::AbstractVector{T}) where {T<:Integer}

Creates a directed simple graph with the combinatorics
of a grid with dimension [d1, d2, ...].

# Examples
```julia-repl
julia> g = di_grid(5,4)
{20, 31} directed Int64 metagraph with Float64 weights defined by :weight (default weight 1.0)
```
"""
function di_grid(dims::AbstractVector{T}) where {T<:Integer}
    g = path_digraph(dims[1])
    for d in dims[2:end]
        g = cartesian_product(path_digraph(d), g)
    end
    g = MetaDiGraph(g)
    set_prop!(g, :type, "grid")
    set_prop!(g, :periodic, false)
    set_prop!(g, :dim, dims)
    if length(dims)==2
        for v in vertices(g)
            out = outneighbors(g, v)
            if length(out)==2
                w1, w2 = sort(out)
                set_prop!(g, v, w2, :dir, "left")
                set_prop!(g, v, w1, :dir, "right")
            elseif length(out)==1
                if v%dims[1]==0
                    set_prop!(g, v, out[1], :dir, "left")
                else
                    set_prop!(g, v, out[1], :dir, "right")
                end
            end
        end
    end
    return g
end

di_grid(m::Int, n::Int) = di_grid([m,n])

"""
    cylinder(dims::AbstractVector{T}) where {T<:Integer}

Creates a directed graph with the combinatorics of cylinder S¹xR^n
of dimension [d1, d2, ...]. The first direction is periodic.

# Examples
```julia-repl
julia> g = cylinder(5,4)
{20, 35} directed Int64 metagraph with Float64 weights defined by :weight (default weight 1.0)
```
"""
function cylinder(dims::AbstractVector{T}) where {T<:Integer}
    g = cycle_digraph(dims[1])
    for d in dims[2:end]
        g = cartesian_product(path_digraph(d), g)
    end
    g = MetaDiGraph(g)
    set_prop!(g, :type, "cylinder")
    set_prop!(g, :periodic, true)
    set_prop!(g, :dim, dims)
    if length(dims)==2
        for v in vertices(g)
            out = outneighbors(g, v)
            if length(out)==2
                w1, w2 = sort(out)
                set_prop!(g, v, w2, :dir, "left")
                set_prop!(g, v, w1, :dir, "right")
            elseif length(out)==1
                if v%dims[1]==0
                    set_prop!(g, v, out[1], :dir, "left")
                else
                    set_prop!(g, v, out[1], :dir, "right")
                end
            end
        end
        set_prop!(g, dims[1]*dims[2], dims[1]*(dims[2]-1)+1, :dir, "right")
    end
    return g
end

cylinder(m::Int, n::Int) = cylinder([m,n])


"""
    zigzag(dims::AbstractVector{T}; periodic=false) where {T<:Integer}

Creates a directed Metagraph with the combinatorics of a mesh over a zigzag-path.
One can see it as a grid where the edges are interchanged with the diagonal
edges. Note that this is implemented only for 2d and can be periodic in the first
direction. One zigzag with 2m points is called zigzag([m, 2]). The second entry
discribes the number of "layers" and the first the number of "ground points".

# Examples
```julia-repl
julia> g = zigzag(5,4;periodic=true)
{20, 30} directed Int64 metagraph with Float64 weights defined by :weight (default weight 1.0)
```
"""
function zigzag(dims::AbstractVector{T}; periodic=false) where {T<:Integer}
    if length(dims)!=2
        println("This is only implemented for 2d. Feel free to generalize this method...")
        return
    end
    m, n = dims
    g = MetaDiGraph(m*n)
    if periodic
        for i=1:n-1
            for v=1:m
                index = (i-1)*m+v
                add_edge!(g, index, toIndex(nw(v, i, m)..., m))
                add_edge!(g, index, toIndex(no(v, i, m)..., m))
                set_prop!(g, index, toIndex(nw(v, i, m)..., m), :dir, "left")
                set_prop!(g, index, toIndex(no(v, i, m)..., m), :dir, "right")
            end
        end
        set_prop!(g, :periodic, true)
    else
        for i=1:n-1
            for v=1:m
                index = (i-1)*m+v
                if i%2==1
                    if v==1
                        add_edge!(g, index, toIndex(no(v, i, m)..., m))
                        set_prop!(g, index, toIndex(no(v, i, m)..., m), :dir, "right")
                    else
                        add_edge!(g, index, toIndex(nw(v, i, m)..., m))
                        add_edge!(g, index, toIndex(no(v, i, m)..., m))
                        set_prop!(g, index, toIndex(nw(v, i, m)..., m), :dir, "left")
                        set_prop!(g, index, toIndex(no(v, i, m)..., m), :dir, "right")
                    end
                else
                    if v!=m
                        add_edge!(g, index, toIndex(nw(v, i, m)..., m))
                        add_edge!(g, index, toIndex(no(v, i, m)..., m))
                        set_prop!(g, index, toIndex(nw(v, i, m)..., m), :dir, "left")
                        set_prop!(g, index, toIndex(no(v, i, m)..., m), :dir, "right")
                    else
                        add_edge!(g, index, toIndex(nw(v, i, m)..., m))
                        set_prop!(g, index, toIndex(nw(v, i, m)..., m), :dir, "left")
                    end
                end
            end
        end
        set_prop!(g, :periodic, false)
    end
    set_prop!(g, :type,"zigzag")
    set_prop!(g, :dim, dims)
    return g
end


function zigzag(m::Int, n::Int; periodic=false)
    return zigzag([m,n]; periodic=periodic)
end


function SimpleGraph(g::MetaDiGraph)
    res = SimpleGraph(nv(g))
    for e in edges(g)
        add_edge!(res, e)
    end
    return res
end


function reverse(g::SimpleDiGraph, e::Edge)
    rem_edge!(g, src(e), dst(e))
    add_edge!(g, dst(e), src(e))
end


function trav_tree(g::MetaDiGraph)
    get_prop(g, :type)=="zigzag" || return bfs_tree(g, 1)
    tree = bfs_tree(SimpleGraph(g), 1)
    for e in edges(tree)
        has_edge(g, e) || reverse(tree, e)
    end
    return tree
end


function left(g::MetaDiGraph, v::Int)
    out = sorted_outneighbors(g, v)
    return out[1]
end


function right(g::MetaDiGraph, v::Int)
    out = sorted_outneighbors(g, v)
    return out[2]
end


function opposite_vertex(g::MetaDiGraph, v::Int)
    try
        left(g, right(g, v))
    catch KeyError
        try
            right(g, left(g, v))
        catch KeyError
        end
    end
end


"""
    set_vprops!(g::MetaDiGraph{T,S}, list::AbstractArray, propname) where {T<:Int, S<:Any}

Set properties on the vertices. Each vertex gets attached an list element with
the same index. The property goes under the name propname.

# Examples
```julia-repl
julia> g = zigzag([5,4])
{20, 30} directed Int64 metagraph with Float64 weights defined by :weight (default weight 1.0)

julia> myProp = collect(1:20);

julia> set_vprops!(mg, myProp, :index)

julia> print_vprop(mg, :index)
The vertex 1 has index == 1
The vertex 2 has index == 2
The vertex 3 has index == 3
The vertex 4 has index == 4
...

```
"""
function set_vprops!(g::MetaDiGraph{T,S}, list::AbstractArray, propname) where {T<:Int, S<:Any}
    set_prop!(g, :properties, propname)
    for v in vertices(g)
        set_prop!(g, v, propname, list[v])
    end
end


"""
    set_vprops!(g::MetaDiGraph{T,S}, func::Function, argname, propname)  where {T<:Int, S<:Any}

Set properties on the vertices. The function func takes as values the properties
under the name argname of the vertex and its two neighbors (if two exist) and
evaluates something which is stored under the name propname.

# Examples
```julia-repl
julia> g = zigzag([3,2])
{6, 6} directed Int64 metagraph with Float64 weights defined by :weight (default weight 1.0)

julia> myfunc(x, x1, x2)=x1+x2-x
myfunc (generic function with 1 method)

julia> index = collect(1:6);

julia> set_vprops!(mg, index, :index)

julia> set_vprops!(mg, myfunc, :index, :myprop)

julia> print_vprop(mg, :myprop)
The vertex 1 has myprop == 9
The vertex 2 has myprop == 7
The vertex 3 has myprop == 8
The vertex 4 does not have the key myprop
The vertex 5 does not have the key myprop
The vertex 6 does not have the key myprop

```
"""
function set_vprops!(g::MetaDiGraph{T,S}, func::Function, argname, propname)  where {T<:Int, S<:Any}
    for v in vertices(g)
        out = outneighbors(g, v)
        if length(out)==2
            w2, w1 = sorted_outneighbors(g, v)
            arg = get_prop(g, v, argname)
            arg1 = get_prop(g, w1, argname)
            arg2 = get_prop(g, w2, argname)
            set_prop!(g, v, propname, func(arg, arg1, arg2))
        end
    end
end


function set_vprops_left!(g::MetaDiGraph{T,S}, func::Function, argname, propname) where {T<:Int, S<:Any}
    for v in vertices(g)
        out = outneighbors(g, v)
        inn = inneighbors(g, v)
        if length(out)>0 && length(inn)>0
            inn = inneighbors(g, v)
            ll = inn[findfirst(x->right(g, x)==v, inn)]
            ur = left(g, v)
            arg = get_prop(g, v, argname)
            arg1 = get_prop(g, ur, argname)
            arg2 = get_prop(g, ll, argname)
            set_prop!(g, v, propname, func(arg, arg1, arg2))
        end
    end
end


function set_vprops_out!(g::MetaDiGraph{T,S}, func::Function, argname, propname; argnum=3)  where {T<:Int, S<:Any}
    if argnum == 3
        for v in vertices(g)
            out = outneighbors(g, v)
            if length(out)==2
                w2, w1 = sorted_outneighbors(g, v)
                arg = get_prop(g, v, argname)
                arg1 = get_prop(g, w1, argname)
                arg2 = get_prop(g, w2, argname)
                set_prop!(g, v, propname, func(arg, arg1, arg2))
            end
        end
    elseif argnum == 2
        for e in edges(g)
            s, t = src(e), dst(e)
            arg1 = get_prop(g, s, argname)
            arg2 = get_prop(g, t, argname)
            set_prop!(g, s, :propname, func(arg1, arg2))
        end
    end
end


function set_vprops_inn!(g::MetaDiGraph{T,S}, func::Function, argname, propname)  where {T<:Int, S<:Any}
    for v in vertices(g)
        inn = inneighbors(g, v)
        if length(inn)==2
            w1, w2 = sorted_innneighbors(g, v)
            arg = get_prop(g, v, argname)
            arg1 = get_prop(g, w1, argname)
            arg2 = get_prop(g, w2, argname)
            set_prop!(g, v, propname, func(arg, arg2, arg1))
        end
    end
end


function rem_vprop!(g::MetaDiGraph, name::Symbol)
    for v in vertices(g)
        rem_prop!(g, v, name)
    end
end


function rem_eprop!(g::MetaDiGraph, name::Symbol)
    for e in edges(g)
        rem_prop!(g, e, name)
    end
end


"""
    print_vprop(mg::MetaDiGraph{T,S}, name) where {T<:Int, S<:Any}

Print all vertices and its value of the property name.

# Examples
```julia-repl
julia> g = zigzag([3,2])
{6, 6} directed Int64 metagraph with Float64 weights defined by :weight (default weight 1.0)

julia> myfunc(x, x1, x2)=x1+x2-x
myfunc (generic function with 1 method)

julia> index = collect(1:6);

julia> set_vprops!(mg, index, :index)

julia> set_vprops!(mg, myfunc, :index, :myprop)

julia> print_vprop(mg, :myprop)
The vertex 1 has myprop == 9
The vertex 2 has myprop == 7
The vertex 3 has myprop == 8
The vertex 4 does not have the key myprop
The vertex 5 does not have the key myprop
The vertex 6 does not have the key myprop

```
"""
function print_vprop(mg::MetaDiGraph{T,S}, name) where {T<:Int, S<:Any}
    for v in vertices(mg)
        try
            println("The vertex $v has $name == $(get_prop(mg, v, name))")
        catch KeyError
            println("The vertex $v does not have the key $name")
        end
    end
end


"""
    print_vprops(mg::MetaDiGraph{T,S}) where {T<:Int, S<:Any}

Print all vertices and all its values. They are stored in Dictionaries.

# Examples
```julia-repl
julia> g = zigzag([3,2])
{6, 6} directed Int64 metagraph with Float64 weights defined by :weight (default weight 1.0)

julia> myfunc(x, x1, x2)=x1+x2-x
myfunc (generic function with 1 method)

julia> index = collect(1:6);

julia> set_vprops!(mg, index, :index)

julia> set_vprops!(mg, myfunc, :index, :myprop)

julia> print_vprops(mg)
props(mg, v) = Dict{Symbol,Any}(:index => 1,:myprop => 9)
props(mg, v) = Dict{Symbol,Any}(:index => 2,:myprop => 7)
props(mg, v) = Dict{Symbol,Any}(:index => 3,:myprop => 8)
props(mg, v) = Dict{Symbol,Any}(:index => 4)
props(mg, v) = Dict{Symbol,Any}(:index => 5)
props(mg, v) = Dict{Symbol,Any}(:index => 6)
```
"""
function print_vprops(mg::MetaDiGraph{T,S}) where {T<:Int, S<:Any}
    for v in vertices(mg)
        @show props(mg, v)
    end
end


function get_vprops(g::MetaDiGraph, propname)
    out = []
    for v in vertices(g)
        try
            push!(out, get_prop(g, v, propname))
        catch KeyError
            println("There are vertices without the property $propname")
        end
    end
    return out
end


function sorted_outneighbors(g::MetaDiGraph, v::Int)
    out = outneighbors(g,v)
    res = []
    for o in out
        if get_prop(g, v, o, :dir)=="left"
            pushfirst!(res, o)
        else
            push!(res, o)
        end
    end
    return res
end


function sorted_innneighbors(g::MetaDiGraph, v::Int)
    inn = inneighbors(g,v)
    res = []
    for o in inn
        if get_prop(g, o, v, :dir)=="left"
            pushfirst!(res, o)
        else
            push!(res, o)
        end
    end
    return res
end


function set_eprops!(g::MetaDiGraph, func::Function, argname::Symbol, propname::Symbol)
    for e in edges(g)
        s, t = src(e), dst(e)
        arg1 = get_prop(g, s, argname)
        arg2 = get_prop(g, t, argname)
        set_prop!(g, e, propname, func(arg1, arg2))
    end
end


function print_eprop(g::MetaDiGraph{T,S}, name) where {T<:Int, S<:Any}
    for e in edges(g)
        println("$e has $name: $(get_prop(g, e, name))")
    end
end


function print_eprops(mg::MetaDiGraph{T,S}) where {T<:Int, S<:Any}
    for e in edges(mg)
        @show props(mg, e)
    end
end


function get_eprop(g::MetaDiGraph, name)
    out = []
    for e in edges(g)
        push!(out, get_prop(g, e, name))
    end
end


function get_lax(g::MetaDiGraph, e::Edge, t::Real)
    lax = copy(get_prop(g, e, :lax))
    if get_prop(g, e, :dir)=="right"
        lax[1, 2]*=exp(t)
        lax[2, 1]*=exp(t)
    else
        lax[1, 2]*=exp(-t)
        lax[2, 1]*=exp(-t)
    end
    lax
end


get_lax(g::MetaDiGraph, v::Int, w::Int, t::Real) = get_lax(g, Edge(v,w), t)


function get_dlax(g::MetaDiGraph, e::Edge, t::Real)
    lax = copy(get_prop(g, e, :lax))
    lax[1, 1] = 0.0
    lax[2, 2] = 0.0
    if get_prop(g, e, :dir)=="right"
        lax[1, 2]*=exp(t)
        lax[2, 1]*=exp(t)
    else
        lax[1, 2]*=(-1)*exp(-t)
        lax[2, 1]*=(-1)*exp(-t)
    end
    lax
end


get_dlax(g::MetaDiGraph, v::Int, w::Int, t::Real)=get_dlax(g, Edge(v,w), t)


function label_boundary!(g::MetaDiGraph)
    set_prop!(g, :boundary, "")
    for v in vertices(g)
        inn = length(inneighbors(g,v))
        out = length(outneighbors(g,v))
        if out<=2
            if inn+out==2
                if out==2
                    set_prop!(g, v, :boundary, "outcorner")
                elseif inn==2
                    set_prop!(g, v, :boundary, "inncorner")
                else
                    set_prop!(g, v, :boundary, "in and out corner")
                end
            elseif inn+out==3
                set_prop!(g, v, :boundary, "on boundary line")
            elseif inn+out==4
                set_prop!(g, v, :boundary, "inner point")
            else
                set_prop!(g, v, :boundary, "umbilic")
            end
        else
            set_prop!(g, v, :boundary, "umbilic")
        end
    end
end


"""
    get_triangles(g::MetaDiGraph)

return a list of triangles. For zigzag path there are extra triangles on the top
and the bottom but this is fine.

# Examples
```julia-repl
julia> g = zigzag(5,4)
{20, 30} directed Int64 metagraph with Float64 weights defined by :weight (default weight 1.0)

julia> get_triangles(g)
30-element Array{Any,1}:
 [1, 6, 10]
 [2, 6, 7]
 [3, 7, 8]
 [4, 8, 9]
 [5, 9, 10]
 [6, 11, 12]
 ⋮
 [15, 9, 10]
 [16, 11, 12]
 [17, 12, 13]
 [18, 13, 14]
 [19, 14, 15]
 [20, 11, 15]
```
"""
function get_triangles(g::MetaDiGraph)::Array{Array{Int, 1}, 1}
    triangles = []
    for v in vertices(g)
        out = outneighbors(g,v)
        v12 = opposite_vertex(g, v)
        v12 != nothing || break
        if length(out)==2
            push!(triangles, [v, out[1], out[2]])
            push!(triangles, [out[1], out[2], v12])
        end
    end
    return triangles
end


"""
    get_quads(g::MetaDiGraph)

return a list of quads.

# Examples
```julia-repl
julia> g = zigzag(5,4)
{20, 30} directed Int64 metagraph with Float64 weights defined by :weight (default weight 1.0)

julia> get_quads(g)
10-element Array{Any,1}:
 [1, 6, 10, 11]
 [2, 6, 7, 12]
 [3, 7, 8, 13]
 [4, 8, 9, 14]
 [5, 9, 10, 15]
 [6, 11, 12, 16]
 [7, 12, 13, 17]
 [8, 13, 14, 18]
 [9, 14, 15, 19]
 [10, 11, 15, 20]
```
"""
function get_quads(g::MetaDiGraph)::Array{Array{Int, 1}, 1}
    quads = []
    for v in vertices(g)
        out = outneighbors(g,v)
        if length(out)==2
            out1 = outneighbors(g, out[1])
            out2 = outneighbors(g, out[2])
            isface = any(x -> x in out1, out2)
            if isface
                indx = findfirst(in(out1), out2)
                push!(quads, [v, out[1], out[2], out2[indx]])
            end
        end
    end
    return quads
end
