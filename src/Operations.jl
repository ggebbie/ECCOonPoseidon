function sum(ma::MeshArrays.gcmarray{T, N, Matrix{T}}) where T<:AbstractFloat where N
    tmp = zeros(T, 1)
    for a in eachindex(ma)
        tmp .+= sum(T, ma.f[a])
    end

    return tmp[1]
end

function vertical_sum(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:AbstractFloat where N
    γ = ma.grid #grid info
    sum_temp = MeshArray(γ,T) #make a flat mesharray
    fill!(sum_temp, zero(T))
    for a in eachindex(ma)
        sum_temp.f[a[1]] .+= ma.f[a] #ijk[1] is the face
    end
    return sum_temp 
end

#must specify dims
flat_sum(x::T; dims = Inf) where T<:Matrix = vec(sum(x, dims = dims))

function lateral_sum(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:AbstractFloat where N
    nz = size(ma, 2)
    faces_sum = zeros(T, nz)
    for a in eachindex(ma)
        faces_sum[a[2]] += sum(ma.f[a])
    end
    return faces_sum 
end



function zonal_sum(ma::MeshArrays.gcmarray{T, 1, Matrix{T}}) where T<:Real
    temp = zeros(T, 270)

    temp[1:270] .+= flat_sum(ma.f[1], dims = 1)
    temp[1:270] .+= flat_sum(ma.f[2], dims = 1) 

    temp[1:270] .+= flat_sum(reverse(ma.f[4], dims = 1), dims = 2)
    temp[1:270] .+= flat_sum(reverse(ma.f[5], dims = 1), dims = 2) 

    return temp 
end

function zonal_sum(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real
    nz = size(ma, 2)
    temp = zeros(T, nz, 270)
    for lvl in 1:nz
        temp[lvl, :] .= ma_zonal_sum(ma[:, lvl])
    end

    return temp 
end

function zonal_average(ma::MeshArrays.gcmarray{T, 1, Matrix{T}},
    weights ::MeshArrays.gcmarray{S, 1, Matrix{S}}) where {T<:Real, S<:Real}

    return zonal_sum(ma .* weights) ./ zonal_sum(weights)
end


function zonal_average(ma::MeshArrays.gcmarray{T, 2, Matrix{T}},
    weights ::MeshArrays.gcmarray{S, 2, Matrix{S}}) where {T<:Real, S<:Real}
    nz = size(ma, 2)

    numerator = zeros(T, nz, 270)
    denom = zeros(T, nz, 270)

    for lvl in 1:nz
        numerator[lvl, :] .= zonal_sum(ma[:, lvl] .* weights[:, lvl])
        denom[lvl, :] .= zonal_sum(weights[:, lvl])
    end
    return temp_num ./ temp_denom
end

function zonal_average(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real
    nz = size(ma, 2)

    numerator = zeros(T, nz, 270)
    denom = zeros(T, nz, 270)

    for lvl in 1:nz
        numerator[lvl, :] .= zonal_sum(ma[:, lvl])
        denom[lvl, :] .= zonal_sum((!iszero).(ma[:, lvl]))
    end
    return temp_num ./ temp_denom
end
