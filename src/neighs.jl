
# exports

export find_neighbors, find_neighbors!, distance


"""
    neighbor_cells(i, j, k, N)

Calculate all neighboring cells for i, j, k cell.
"""
function neighbor_cells(i, j, k, N)
    a = Vector{Int}()
    for kk in k-1:k+1
        for jj in j-1:j+1
            for ii in i-1:i+1
                if (0<ii*jj*kk) && (ii<=N[1]) && (jj<=N[2]) && (kk<=N[3])
                    push!(a, cell_number(ii, jj, kk, N))
                end
            end
        end
    end
    return a
end

"""
    cell_number(i, j, k, N)

Calculates cell number for given i, j, k cell indices (when in a list).
"""
function cell_number(i, j, k, N)
    return i+(j-1)*N[1]+(k-1)*N[1]*N[2]
end

"""
    get_cells(x::Array{Float64, 2}, horizon::Float64)

Fill cells with particles.
"""
function get_cells(x::Array{Float64, 2}, horizon::Float64)
    _min = minimum(x, dims=2)
    _max = maximum(x, dims=2)
    N = Int.(1 .+ floor.((_max-_min)/horizon))
    cells = [Vector{Int}() for i in 1:prod(N)]
    cell_neighs = Vector{Vector{Int}}(undef, prod(N))
    for k in 1:N[3]
        for j in 1:N[2]
            for i in 1:N[1]
                cell_neighs[cell_number(i, j, k, N)] = neighbor_cells(i, j, k, N)
            end
        end
    end
    for i in 1:size(x, 2)
        ii, jj, kk = Int.(1 .+ floor.((x[:, i].-_min)/horizon))
        push!(cells[cell_number(ii, jj, kk, N)], i)
    end
    return cells, cell_neighs
end


"""
    find_neighbors!(neighbors::Array{Int64, 2}, x::Array{Float64, 2}, horizon::Float64)

Calculate neighbors for each particle (inplace).
"""
function find_neighbors!(neighbors::Array{Int64, 2}, x::Array{Float64, 2}, horizon::Float64)

    cells, cell_neighs = get_cells(x, horizon)
    for cell_i in 1:length(cells)
        for ca_id in 1:length(cells[cell_i])
            ind = 1
            ca = cells[cell_i][ca_id]
            a1, b1, c1 = x[1, ca], x[2, ca], x[3, ca]
            for neigh_id in 1:length(cell_neighs[cell_i])
                neighs = cell_neighs[cell_i][neigh_id]
                for fa_id in 1:length(cells[neighs])
                    fa = cells[neighs][fa_id]
                    a2, b2, c2 = x[1, fa], x[2, fa], x[3, fa]
                    if 1e-4<((a1-a2)*(a1-a2)+(b1-b2)*(b1-b2)+(c1-c2)*(c1-c2))<horizon^2
                        if ind<=size(neighbors)[1]
                            neighbors[ind, ca] = fa
                        else
                            error("Neighbors are more than capacity allocated. Please increase neighbors maximum capacity.")
                        end
                        ind += 1
                    end
                end
            end
        end
    end
    neighbors
end


"""
    find_neighbors(x::Array{Float64, 2}, horizon::Float64)

Calculate neighbors for each particle.
"""
function find_neighbors(x::Array{Float64, 2}, horizon::Float64, max_neigh::Int64; hard_max::Int64=0)::Array{Int64, 2}
    max_neigh = max(hard_max, min(max_neigh, size(x, 2)))
    neighbors = zeros(Int64, max_neigh, size(x, 2))
    find_neighbors!(neighbors, x, horizon)
    max_row = max(hard_max, maximum(sum(neighbors.!=0,dims=1)))
    new_neigh = neighbors[1:max_row,:]
    replace!(new_neigh,0=>999999999)
    sort!(new_neigh, dims=1)
    replace!(new_neigh,999999999=>0)
end
