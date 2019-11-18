include("quaternion.jl")

module Ddg

using ..Quaternions: Quaternion, qToR3

export propagate_zigzag, periodic_zigzag, ne, nw, sw, se, verts_zigzag, faces_zigzag

function dual_edge(F_1::Quaternion, F::Quaternion)
    l = abs2(F_1-F)
    if l==0
        print("Empty edge")
        F_1-F
    else
        (F_1-F)/l
    end
end

function dualize_grid(grid::Array{Quaternion, 2}, start=zero(Quaternion))
    m, n = size(grid)
    dual_grid = zeros(Quaternion, m, n)
    dual_grid[1,1] = start
    for i=1:m-1
        dual_grid[i+1, 1] = dual_grid[i, 1] + dual_edge(grid[i+1,1], grid[i,1])
    end
    for i=1:m, j=1:n-1
        dual_grid[i,j+1] = dual_grid[i,j] - dual_edge(grid[i,j+1], grid[i,j])
    end
    return dual_grid
end

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

function propagate_zigzag(
    layer1::Array{Quaternion{T}, 1},
    layer2::Array{Quaternion{S}, 1},
    func::Function,
    n::Int,
) where {T,S <:Real}
    m1 = length(layer1)
    m2 = length(layer2)
    if abs(m1-m2)>1
        println("Dimension mismatch for zigzag")
        throw(DomainError)
    end
    if m1!=m2
        println("Up to now propagate_zigzag is only
        implemented for periodic start conditions")
        return
    end
    output = zeros(Quaternion, (m1, n+2))
    #output = Array{Any, 2}(undef, m1, n+2)
    output[:, 1] = layer1
    output[:, 2] = layer2
    for j=1:n, i=1:m1
        output[i, j+2] = func(
            output[i, j],
            output[no(i, j, m1)...],
            output[nw(i, j, m1)...],
        )
    end
    output
end

function periodic_zigzag(points::Array{Quaternion, 2})
    m, n = size(points)
    verts = reshape(qToR3.(points), (1,m*n))
    verts = [verts[i][j] for i=1:m*n, j=1:3]
    faces = zeros(Int64, (2*m*(n-2),3))
    counter = 1
    for j=1:n-2, i=1:m
        if j%2==1
            if i==1
                faces[counter,:] = [(j-1)*m+i,j*m+1,(j+1)*m]
                faces[counter+1,:] = [(j+1)*m+i,j*m+1,(j+1)*m]
            else
                faces[counter,:] = [(j-1)*m+i, i+j*m-1, i+1+j*m-1]
                faces[counter+1,:] = [(j+1)*m+i, i+j*m-1, i+1+j*m-1]
            end
        else
            if i==m
                faces[counter,:] = [j*m, j*m+1, (j+1)*m]
                faces[counter+1,:] = [(j+2)*m, j*m+1, (j+1)*m]
            else
                faces[counter,:] = [(j-1)*m+i, j*m+i, j*m+i+1]
                faces[counter+1,:] = [(j+1)*m+i,j*m+i, j*m+i+1]
            end
        end
        counter+=2
    end
    return verts, faces
end

function verts_zigzag(points::Array{Quaternion, 2})
    m, n = size(points)
    verts = reshape(qToR3.(points), (1,m*n))
    verts = [verts[i][j] for i=1:m*n, j=1:3]
    return verts
end

function faces_zigzag(m::Int, n::Int)::Array{Int64, 2}
    faces = zeros(Int64, (2*m*(n-2),3))
    counter = 1
    for i=1:m
        for j=1:n
            if j<n-1
                i1, uj = nw(i, j, m)
                i2, _ = no(i, j, m)
                faces[counter,:] = [(j-1)*m+i, (uj-1)*m+i1, (uj-1)*m+i2]
                counter+=1
            end
            if j>2
                i1, dj = sw(i,j,m)
                i2, _ = so(i,j,m)
                faces[counter,:] = [(j-1)*m+i, (dj-1)*m+i1, (dj-1)*m+i2]
                counter+=1
            end
        end
    end
    faces
end
end #end module
