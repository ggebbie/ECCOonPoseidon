#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings, GibbsSeaWater, Plots


function extract_sX(ETAN::MeshArrays.gcmarray{T, 1, Matrix{T}}, 
    X::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    inv_depths::MeshArrays.gcmarray{T, 1, Matrix{T}}) where T 
    sX = similar(X)
    s1 = ETAN .* inv_depths
    s1 .+= 1
    sX = X .* s1
    return sX
end

function pgrid(z,lat=30)

    # if p doesn't exist as argument but z does,
    # transfer z to p.
    if !isempty(z)
        # shift depth to pressure
        #use 30deg lat as the default location to convert depth to pressure, from GibbsSeaWater.jl
        p = gsw_p_from_z.(-z, lat) #z must be postive
    else
        # if p and z don't exist, error. 
        error("Input vertical coordinate, pressure or depth, not found")
    end
end

function volume_average_by_depth(ds::MeshArray, ΔV::MeshArray)
    nz = size(ds, 2)
    vol_avg = zeros(Float32, nz)
    V = zeros(Float32, nz)

    for k=1:nz
        V[k] = Float32(sum(ΔV[:, k]))
    end

    for ff=1:5, k=1:nz
        vol_avg[k] += sum(ds[ff, k] .* ΔV[ff, k]) / V[k]
    end
    return vol_avg
end

#define area averaged depth (does this make sense?)
function area_average_by_depth(ds::MeshArray, ΔA::MeshArray)
    nz = size(ds, 2)
    vol_avg = zeros(Float32, nz)
    A = sum(ΔA)

    for ff=1:5, k=1:nz
        ds_mask = ds[ff, k]
        ds_mask[isnan.(ds_mask)] .= 0.0
        vol_avg[k] += sum(ds_mask .* ΔA[ff]) / A
    end

    return vol_avg
end

function area_average_by_σ(ds::MeshArray, ΔA::MeshArray)
    nz = size(ds, 2)
    vol_avg = zeros(Float32, nz)
    
    for k=1:nz 
        ΔA_mask = deepcopy(ΔA) #this should shrink with changing isopycnal surfaces
        for ff=1:5
            ds_mask = deepcopy(ds[ff, k])
            drypts = isnan.(ds_mask)
            ds_mask[drypts] .= 0.0 #better way to do this? 
            ΔA_mask[ff][drypts] .= 0.0 
            vol_avg[k] += sum(ds_mask .* ΔA[ff]) 
        end
        vol_avg[k] = vol_avg[k] / sum(ΔA_mask)
    end
    vol_avg[vol_avg .== 0.0] .= NaN
    return vol_avg
end

function depth_average_σ(ds::MeshArray)
    depth_avg = MeshArray(γ,Float32); 
    fill!(depth_avg, 0.0)
    nz = size(ds, 2)
    for ff=1:5, k=1:nz
        depth_avg[ff] .+= (ds[ff, k]) ./ nz
    end
    return depth_avg
end

function zonal_avg_3D(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}) where {T<:Real}
    nlev = size(ma, 2)
    temp_num = zeros(nlev, 270)
    temp_denom = zeros(nlev, 270)

    weights = similar(ma[:, 1]); fill!(weights, 0.0)

    for lvl in 1:nlev
        ma_copy = ma[:, lvl]
        for ff = 1:5
            where_nan = (isnan).(ma_copy[ff])
            ma_copy[ff][where_nan] .= 0.0
            weights[ff] .= (!).(where_nan)
        end

        temp_num[lvl, :] .= ma_zonal_sum(ma_copy .* weights)
        temp_denom[lvl, :] .= ma_zonal_sum(weights)
    end
    temp_denom[temp_denom .== 0] .= NaN
    return temp_num ./ temp_denom
end

function ma_zonal_sum(ma::MeshArrays.gcmarray{T, 1, Matrix{T}}) where T<:Real
    temp = zeros(270)
    for ff in eachindex(ma) 
	    if (ff==1) || (ff == 2)
            temp[1:270] .+= @views sum( ma[ff], dims = 1)[:] 
        elseif ff == 4 || ff == 5 
            temp[1:270] .+= @views sum( reverse(ma[ff], dims = 1), dims = 2)[:] 
        # elseif ff == 3
        #     temp[271:end] .+= @views sum( reverse(ma[ff], dims = 1), dims = 2)[:]
        end
    end
    return temp 
end

function zonal_sum_3D(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real
    nz = size(ma, 2)
    temp = zeros(nz, 270)
    for lvl in 1:nz
        ma_copy = ma[:, lvl]
        for ff = 1:5
            where_nan = (isnan).(ma_copy[ff])
            ma_copy[ff][where_nan] .= 0.0
        end

        temp[lvl, :] .= ma_zonal_sum(ma[:, lvl])
    end

    return temp 
end

σ1datadir(f="") = datadir() * "/sigma1/" * f
σ2datadir(f="") = datadir() * "/sigma2/" * f
