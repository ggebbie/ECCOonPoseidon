"""
    region_mask(wet_pts, λ, ϕ, latlim, lonlim)

Create a regional mask based on latitude and longitude limits.

# Arguments
- `wet_pts`: Ocean mask obtained from `wet_pts(Γ)`.
- `λ`: Latitude mesh array for each face of the grid.
- `ϕ`: Longitude mesh array for each face of the grid.
- `latlim`: Latitude limits (tuple) defining the region of interest.
- `lonlim`: Longitude limits (tuple) defining the region of interest.

# Output
- Returns a regional mask where entries are 1.0 
inside the specified latitude and longitude limits and 0.0 outside.
"""

function region_mask(wet_pts, λ, ϕ, latlim, lonlim)
    reg_msk = similar(wet_pts)
    for ff in 1:5
        in_lon = lonlim[1] .< λ[ff] .< lonlim[2]
        in_lat = latlim[1] .< ϕ[ff] .< latlim[2]
        reg_msk.f[ff] .= wet_pts.f[ff] .* in_lon .* in_lat
    end
    return reg_msk
end
