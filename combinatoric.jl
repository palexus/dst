using LightGraphs, MetaGraphs

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
julia> di_grid([5,4])
{20, 31} directed simple Int64 graph
```
"""
function di_grid(dims::AbstractVector{T}) where {T<:Integer}
    g = path_digraph(dims[1])
    for d in dims[2:end]
        g = cartesian_product(path_digraph(d), g)
    end
    mg = MetaDiGraph(g)
    set_prop!(mg, :type, "grid")
    set_prop!(mg, :dim, dims)
    if length(dims)==2
        lrlabel = [("left", "right") for v=1:(dims[1]*dims[2])]
        set_eprops!(mg, lrlabel, :dir)
    end
    return mg
end

di_grid(m::Int, n::Int) = di_grid([m,n])

"""
    cylinder(dims::AbstractVector{T}) where {T<:Integer}

Creates a directed graph with the combinatorics of cylinder S¹xR^n
of dimension [d1, d2, ...]. The first direction is periodic.

# Examples
```julia-repl
julia> cylinder([5,4])
{20, 35} directed simple Int64 graph
```
"""
function cylinder(dims::AbstractVector{T}) where {T<:Integer}
    g = cycle_digraph(dims[1])
    for d in dims[2:end]
        g = cartesian_product(path_digraph(d), g)
    end
    mg = MetaDiGraph(g)
    set_prop!(mg, :type, "cylinder")
    set_prop!(mg, :dim, dims)
    if length(dims)==2
        lrlabel = [("left", "right") for v=1:(dims[1]*dims[2])]
        set_eprops!(mg, lrlabel, :dir)
    end
    return mg
end

cylinder(m::Int, n::Int) = cylinder([m,n])


"""
    zigzag(dims::AbstractVector{T}; periodic=true) where {T<:Integer}

Creates a directed Metagraph with the combinatorics of a mesh over a zigzag-path.
One can see it as a grid where the edges are interchanged with the diagonal
edges. Note that this is implemented only for 2d and is periodic in the first
direction. One zigzag with 2m points is called zigzag([m, 2]). The second entry
discribes the number of "layers" and the first the number of "ground points".

# Examples
```julia-repl
julia> zigzag([7,4])
{28, 42} directed simple Int64 graph
```
"""
function zigzag(dims::AbstractVector{T}; periodic=true) where {T<:Integer}
    if !periodic
        println("only implemented for periodic zigzag")
        return
    end
    if length(dims)!=2
        println("This is only implemented for 2d. Feel free to generalize this method...")
        return
    end
    m, n = dims
    graphs = [SimpleDiGraph{T}(m) for j=1:n]
    g = graphs[1]
    for (i,gi) in enumerate(graphs[2:end])
        g = blockdiag(g, gi)
        for v in vertices(gi)
            index = (i-1)*m+v
            add_edge!(g, index, toIndex(nw(v, i, m)..., m))
            add_edge!(g, index, toIndex(no(v, i, m)..., m))
        end
    end
    mg = MetaDiGraph(g)
    set_prop!(mg, :type, "zigzag")
    set_prop!(mg, :dim, dims)
    if length(dims)==2
        lrlabel = [("left", "right") for v=1:(dims[1]*dims[2])]
        set_eprops!(mg, lrlabel, :dir)
    end
    return mg
end

function zigzag(m::Int, n::Int)
    return zigzag([m,n])
end


"""
    traversion_tree(graph::SimpleDiGraph{T}) where {T<:Int}

Creates a spanning tree of the given graph.

# Examples
```julia-repl
julia> g = di_grid([3,2])
{6, 7} directed simple Int64 graph

julia> t = traversion_tree(g)
{6, 5} directed simple Int64 graph

julia> collect(edges(t))
5-element Array{LightGraphs.SimpleGraphs.SimpleEdge{Int64},1}:
 Edge 1 => 2
 Edge 1 => 4
 Edge 2 => 3
 Edge 2 => 5
 Edge 3 => 6

```
"""
function traversion_tree(graph::SimpleDiGraph{T}) where {T<:Int}
    g = SimpleDiGraph(length(vertices(graph)))
    visited = []
    for e in edges(graph)
        if !in(dst(e), visited)
            push!(visited, dst(e))
            add_edge!(g, e)
        end
    end
    g
end


"""
    set_vprops!(g::MetaDiGraph{T,S}, list::AbstractArray, propname) where {T<:Int, S<:Any}

Set properties on the vertices. Each vertex gets attached an list element with
the same index. The property goes under the name propname.

# Examples
```julia-repl
julia> g = zigzag([5,4])
{20, 30} directed simple Int64 graph

julia> mg = MetaDiGraph(g)
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
{6, 6} directed simple Int64 graph

julia> mg = MetaDiGraph(g)
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
            w2, w1 = sort(outneighbors(g, v))
            arg = get_prop(g, v, argname)
            arg1 = get_prop(g, w1, argname)
            arg2 = get_prop(g, w2, argname)
            set_prop!(g, v, propname, func(arg, arg1, arg2))
        end
    end
end


"""
    print_vprop(mg::MetaDiGraph{T,S}, name) where {T<:Int, S<:Any}

Print all vertices and its value of the property name.

# Examples
```julia-repl
julia> g = zigzag([3,2])
{6, 6} directed simple Int64 graph

julia> mg = MetaDiGraph(g)
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
{6, 6} directed simple Int64 graph

julia> mg = MetaDiGraph(g)
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

function set_eprops!(g::MetaDiGraph{T,S}, data::AbstractArray, propname;
                     type="none", dim=zeros(2)) where {T<:Int, S<:Any}
    m, n = dim
    try
        type = get_prop(g, :type)
        m, n = get_prop(g, :dim)
    catch KeyError
        if dim==zeros(2)
            _, m = sort(outneighbors(g, 1))
            m-=1
            n = convert(Int, length(vertices(g))÷m)
        end
    end
    if type=="grid" || type=="none"
        for v in vertices(g)
            out = outneighbors(g, v)
            if length(out)==2
                w1, w2 = sort(out)
                set_prop!(g, v, w2, propname, data[v][1])
                set_prop!(g, v, w1, propname, data[v][2])
            elseif length(out)==1
                if v%m==0
                    set_prop!(g, v, out[1], propname, data[v][1])
                else
                    set_prop!(g, v, out[1], propname, data[v][2])
                end
            end
        end
    elseif type=="cylinder"
    elseif type=="zigzag"
        for v in vertices(g)
            out = outneighbors(g, v)
            if length(out)==2
                w1, w2 = sort(out)
                mode = (mod(v,m)==1 && (v÷m)%2==0) || (mod(v,m)==0 && (v÷m)%2==0)
                if mode
                    set_prop!(g, v, w2, propname, data[v][1])
                    set_prop!(g, v, w1, propname, data[v][2])
                else
                    set_prop!(g, v, w1, propname, data[v][1])
                    set_prop!(g, v, w2, propname, data[v][2])
                end
            end
        end
    end
end

function set_eprops_from_v(g::MetaDiGraph{T,S}, func::Function, propname) where {T<:Int, S<:Any}
    for e in edges(g)
        set_prop!(g, e, propname, func(src(e), dst(e)))
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

function get_eprops(g::MetaDiGraph, name)
    out = []
    for e in edges(g)
        out.append(get_prop(g, e, name))
    end
end

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
        inn = inneighbors(g,v)
        out = outneighbors(g,v)
        if length(out)==2
            push!(triangles, [v, out[1], out[2]])
        end
        if length(inn)==2
            push!(triangles, [v, inn[1], inn[2]])
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


g = zigzag([5,4])

tri = get_triangles(g)

ed = collect(edges(g))

gr = di_grid([5,4])

tri2 = get_triangles(gr)
cyl = cylinder([5,4])

quads = get_quads(cyl)
