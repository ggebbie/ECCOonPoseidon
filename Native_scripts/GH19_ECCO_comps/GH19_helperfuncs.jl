
function zonalgriddist(lon, lat)
    dx = similar(lat)
    for j in eachindex(lat)
        dx[j] = haversine((lon[1],lat[j])
                         ,(lon[2],lat[j]))
    end
    return dx
end

function cellarea(lon, lat)
    dx = zonalgriddist(lon, lat)
    dy = haversine((lon[1],lat[1])
                  ,(lon[1],lat[2]))
    dy = repeat([dy], length(lon) )
    area = dx .* dy'
    return area
end

function GH19_cell_volumes(depth, lons, lats)
#construct the masked volume 
    nz = length(depth)
    areas = cellarea(lons, lats)
    volumes = zeros(Float32, length(depth), size(areas)[1], size(areas)[2])
    Δd = depth[2:end] .- depth[1:end-1]
    Δd = vcat(Δd, 500)
    for k in 1:nz
        volumes[k, :, :] .= Δd[k] .* areas
    end
    return volumes
end
    