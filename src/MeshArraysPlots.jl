#cartopy projections
crs = ECCOonPoseidon.cartopy.crs
function replace_zeros!(ma::MeshArray, fill_val::Real)
    for a in eachindex(ma)
        ma.f[a][iszero.(ma.f[a])] .= fill_val
    end
end

function wrap_λ(λ::MeshArray)
    λ_wrap = deepcopy(λ)
    [λ_wrap.f[ff][λ.f[ff] .< 0.0] .+= 360 for ff = 1:5]
    return λ_wrap
end

function LLCcropC360(x::MeshArray, γ; modify_λ = false)
    xcrop = LLCcropC(x, γ)
    xwrap = vcat(xcrop[181:end, :], xcrop[1:180, :])
    if modify_λ
        xwrap[xwrap .< 0.0] .+= 360
    end
    return xwrap
end


function pcolormesh_ma(ax, λ, ϕ, var, cmap; add_gridlines = false, bounds = nothing)
    projPC = crs.PlateCarree()
    
    #obtain maximal value 
    maxval = maximum(abs.(var))
    isnothing(bounds) && (vmin = -maxval; vmax = maxval)
    (!isnothing)(bounds) && (vmin = bounds[1]; vmax = bounds[2])

    CF = Any[]
    for ff = 1:5
        cf = ax.pcolormesh(λ[ff], ϕ[ff],  var[ff], vmin = vmin, vmax = vmax, 
                            transform=projPC, cmap = cmap, 
                            shading = "nearest")
        push!(CF, cf)
    end
    if add_gridlines
        gl = ax.gridlines(crs=projPC, draw_labels=true, linewidth=2, 
        color="gray", alpha=0, linestyle="--")
        gl.top_labels = false; gl.right_labels = false
    end
    return CF[1]
end
