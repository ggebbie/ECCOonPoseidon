region2 = "PAC56"
PAC56_mask = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region2, extent = "full")


function latitude_section_mask(ϕ, ϕ_mask_min, PAC_mask)
    ϕ_min_mask = ϕ .> -Inf
    ϕ_min_mask.f[3] .= 0.0
    [ϕ_min_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- ϕ_mask_min, 0.5) .* PAC_mask.f[ff] for ff in 1:2]
    [ϕ_min_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- ϕ_mask_min, 0.5) .* PAC_mask.f[ff] for ff in 4:5]
    return ϕ_min_mask
end

function longitude_section_mask(λ, λ_mask_min, PAC_mask)
    λ_min_mask = λ .> -Inf
    λ_min_mask.f[3] .= 0.0
    [λ_min_mask.f[ff] .= abs_dist.(λ.f[ff] .- λ_mask_min, 0.5) .* PAC_mask.f[ff] for ff in 1:2]
    [λ_min_mask.f[ff] .= abs_dist.(λ.f[ff] .- λ_mask_min, 0.5) .* PAC_mask.f[ff] for ff in 4:5]
    return λ_min_mask
end

function plot_basin_mask!(λ, ϕ, mask, proj)    
    for ff in 1:5
        ax.pcolormesh(λ[ff], ϕ[ff],  mask[ff],
        vmin = 0, vmax = 2, shading="nearest", transform=proj, 
        rasterized = true, cmap = "binary")               
    end

end