"""
    PAC_mask(Γ, basins, basin_list, ϕ, λ; region = "", extent = "model-defined", include_bering = true)

Create a mask for the Pacific Ocean region in a global ocean model.

# Arguments
- `Γ`: The global ocean model grid.
- `basins`: A matrix representing the basins in the model.
- `basin_list`: A list of basin names.
- `ϕ`: Latitude grid.
- `λ`: Longitude grid.
- `region`: A string specifying the region of interest. Default is an empty string, which corresponds to the entire Pacific Ocean.
- `extent`: A string specifying the extent of the region. Default is "model-defined". If "full", includes marginal seas.
- `include_bering`: A boolean specifying whether to include the Bering Sea. Default is true.

# Returns
- `PAC_msk`: A mask for the specified region of the Pacific Ocean.

# Example
```julia
"""
function PAC_mask(Γ, basins, basin_list, ϕ, λ; region = "", 
    extent = "model-defined", include_bering = true)
    ocean_mask = wet_pts(Γ)
    basin_name="Pacific"
    # IndID = findall(basin_list.=="Indian")[1]
    basinID=findall(basin_list.==basin_name)[1]
    Okhotsk = findall(basin_list.=="Okhotsk Sea")[1]      
    basinID2 = findall(basin_list.== "Japan Sea")[1]   
    basinID3 = findall(basin_list.== "East China Sea")[1]   
    basinID4 = findall(basin_list.== "South China Sea")[1]   
    basinID5 =  findall(basin_list.== "Bering Sea")[1]   
    PAC_msk=similar(basins)
    bounds = [0.0, 0.0]
    if region == "NPAC"
        # bounds[1] = 51.0
        bounds[1] = 70.0
        bounds[2] = 23.0
    elseif region == "ESPAC"
        bounds[1] = 17.0
        bounds[2] = -35.0
    elseif region == "NPAC30"
        bounds[1] = 70.0
        bounds[2] = 30.0  
    elseif region == "SPAC"
        bounds[1] = 25.0
        bounds[2] = -40.0
    elseif region == "EPAC"
        bounds[1] = 25.0
        bounds[2] = -25.0
    elseif region == "PAC56"
        bounds[1] = 65.0
        bounds[2] = -56.0
    elseif region == "PAC"
        bounds[1] = 70.0
        bounds[2] = -40.0
    else 
        bounds[1] = 70.0
        bounds[2] = -40.0
    end
    full_extent = (extent == "full")

    for ff in 1:5
        above_SO = (bounds[2] .< ϕ[ff] .< bounds[1]) #make this go to bering strait
        marginal_seas = (basins[ff] .== Okhotsk) .+ (basins[ff] .== basinID2) .+ 
                        (basins[ff] .== basinID3) .+  (basins[ff] .== basinID4)
        Okhotsk_sea = (basins[ff] .== Okhotsk)
        SouthChinaSea = 0 
        BeringSea     = (basins[ff] .== basinID5) .* include_bering
        model_PAC     = (basins[ff] .== basinID) 
        temp = model_PAC .+ BeringSea .+ Okhotsk_sea .+ (full_extent .* marginal_seas)
        temp = ocean_mask[ff] .* temp .* above_SO 
        PAC_msk[ff] .= temp 
    end

    PAC_msk[findall(PAC_msk .== 0)] = 0f0
    PAC_msk[findall(PAC_msk .> 0)] = 1f0 

    return PAC_msk
end

"""
    get_min_lat(ϕ, mask)

Get the minimum latitude from a masked latitude grid.

# Arguments
- `ϕ`: Latitude grid.
- `mask`: Mask to apply to the grid.

# Returns
- Minimum latitude from the masked grid.
"""

function get_min_lat(ϕ, mask)
    ϕ_mask_Inf = ϕ .* mask; 
    for ff = 1:5
        ϕ_mask_Inf.f[ff][ϕ_mask_Inf[ff] .== 0] .= Inf
    end
    return minimum(ϕ_mask_Inf)
end

"""
    get_max_lat(ϕ, mask)

Get the maximum latitude from a masked latitude grid.

# Arguments
- `ϕ`: Latitude grid.
- `mask`: Mask to apply to the grid.

# Returns
- Maximum latitude from the masked grid.
"""

function get_max_lat(ϕ, mask)
    ϕ_mask_Inf = ϕ .* mask; 
    [ϕ_mask_Inf.f[ff][ϕ_mask_Inf[ff] .== 0] .= -Inf for ff = 1:5]
    ϕ_mask_max = maximum(ϕ_mask_Inf)
    return ϕ_mask_max
end

"""
    within_lon(λ, λ1, λ2)

Check if longitudes in a grid are within a specified range.

# Arguments
- `λ`: Longitude grid.
- `λ1`: Lower bound of the range.
- `λ2`: Upper bound of the range.

# Returns
- A mask where true values indicate longitudes within the range.
"""
function within_lon(λ, λ1, λ2)
    λ_mask = λ .> Inf
    for ff in eachindex(λ)
        λ_mask.f[ff] .= λ1 .< λ.f[ff] .< λ2
    end
    return λ_mask
end
"""
    get_cs_and_sn(γ)

Get cosine and sine arrays for a given grid.

# Arguments
- `γ`: The grid.

# Returns
- `cs`: Cosine array.
- `sn`: Sine array.
"""
function get_cs_and_sn(γ)
    cs = MeshArray(γ,Float32); sn = MeshArray(γ,Float32)
    fill!(cs, 0.0); fill!(sn, 0.0)
    cs.f[1] .= 1; cs.f[2] .= 1
    sn.f[4] .= -1; sn.f[5] .= -1
    return cs, sn
end
"""
    rotate_UV_native(uvel, vvel, cs, sn)

Rotate u and v velocity components to eastward and northward components.

# Arguments
- `uvel`: U velocity component.
- `vvel`: V velocity component.
- `cs`: Cosine array.
- `sn`: Sine array.

# Returns
- `evel`: Eastward velocity component.
- `nvel`: Northward velocity component.
"""
function rotate_UV_native(uvel::MeshArrays.gcmarray{T,N,Matrix{T}},
    vvel::MeshArrays.gcmarray{T,N,Matrix{T}},
    cs::MeshArrays.gcmarray{T,1,Matrix{T}}, sn::MeshArrays.gcmarray{T,1,Matrix{T}}) where T<:AbstractFloat where N
    zeros(T, )
    evel = T.(similar(uvel))
    nvel = T.(similar(vvel))
    for ff in eachindex(evel)
        evel.f[ff] .= uvel.f[ff] .* cs.f[ff[1]] .- vvel.f[ff] .*sn.f[ff[1]]
        nvel.f[ff] .= uvel.f[ff] .* sn.f[ff[1]] .+ vvel.f[ff] .*cs.f[ff[1]]
    end
    return evel,nvel
end


"""
    get_ϕ_max_min_mask(region, Γ, λ, ϕ, basins, basin_list)

Get masks for maximum and minimum latitudes for a specified region.

# Arguments
- `region`: The region of interest.
- `Γ`: The global ocean model grid.
- `λ`: Longitude grid.
- `ϕ`: Latitude grid.
- `basins`: A matrix representing the basins in the model.
- `basin_list`: A list of basin names.

# Returns
- `ϕ_min_mask`: Mask for minimum latitudes.
- `ϕ_max_mask`: Mask for maximum latitudes.
"""

function get_ϕ_max_min_mask(region, Γ, λ, ϕ, basins, basin_list)
    abs_dist(x, r) = abs(x) < r

    if region == "NPAC"
        PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
        region2 = "PAC56"
        PAC56_mask = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region2)
    
        ϕ_mask_min = Float32(get_min_lat(ϕ, PAC_msk)); ϕ_min_mask = ϕ .> -Inf
        ϕ_min_mask.f[3] .= 0.0
        [ϕ_min_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- 23.8, 0.1) .* PAC56_mask.f[ff] for ff in 1:2]
        [ϕ_min_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- 22.9, 0.1) .* PAC56_mask.f[ff] for ff in 4:5]

        ϕ_mask_max = Float32(get_max_lat(ϕ, PAC_msk)); ϕ_max_mask = ϕ .> -Inf; 
        # [ϕ_max_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- 51.0, 0.1) .* PAC56_mask.f[ff] for ff in 1:2]
        # [ϕ_max_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- 50.3, 0.1) .* PAC56_mask.f[ff] for ff in 4:5]
        return ϕ_min_mask, 0.0 .* ϕ_max_mask
    elseif region == "ESPAC"
        PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
        region2 = "PAC56"
        ϕ_min_mask.f[3] .= 0.0
        PAC56_mask = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region2)
        ϕ_mask_min = Float32(get_min_lat(ϕ, PAC_msk)); ϕ_min_mask = ϕ .> -Inf

        [ϕ_min_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- (-34.3), 0.1) .* PAC56_mask.f[ff] for ff in 1:2]
        [ϕ_min_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- (-35.1), 0.1) .* PAC56_mask.f[ff] for ff in 4:5]
        ϕ_mask_max = Float32(get_max_lat(ϕ, PAC_msk)); ϕ_max_mask = ϕ .> -Inf; 

        [ϕ_max_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- 17.2, 0.1) .* PAC56_mask.f[ff] for ff in 1:2]
        [ϕ_max_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- 16.3, 0.1) .* PAC56_mask.f[ff] for ff in 4:5]
        return ϕ_min_mask, ϕ_max_mask
    end
    
end



