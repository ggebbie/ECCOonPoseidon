"""
Extracts surface wind stress components τx and τy from diagnostic files.

Parameters:
- `diagpath`: Dictionary containing paths to different experiments.
- `expname`: Name of the experiment.
- `fname`: File name of the diagnostic data.
- `γ`: Grid information.

Returns:
- `τx`: Array of τx component.
- `τy`: Array of τy component.
"""
function extract_ocnTAU(diagpath, expname, fnameτ, γ)
    flux_forcing_exps = ["seasonalclimatology", "climatological_tau", "clim_tau_iter0"]
    #could use 
    #  read_bin(diagpath[expname]*fnameτ,Float32,γ) and pull the last two elements
    if expname ∈ flux_forcing_exps
        @time EXF = γ.read(diagpath[expname]*fnameτ, MeshArray(γ, Float32, 10))
        τx = EXF[:, 9]
        τy = EXF[:, 10]
        return τx, τy
    else
        @time EXF = γ.read(diagpath[expname]*fnameτ, MeshArray(γ, Float32, 15))
        τx = EXF[:, 14]
        τy = EXF[:, 15]
        return τx, τy
    end
end

"""
Extracts Eulerian velocities u, v, and w from diagnostic files.

Parameters:
- `diagpath`: Dictionary containing paths to different experiments.
- `expname`: Name of the experiment.
- `fnameuvw`: File name of the diagnostic data.
- `γ`: Grid information.

Returns:
- `u`: Array of u-component of velocity.
- `v`: Array of v-component of velocity.
- `w`: Array of w-component of velocity.
"""
function extract_eulerian_velocities(diagpath, expname, fnameuvw, γ)
    UV = γ.read(diagpath[expname]*fnameuvw, MeshArray(γ, Float32, 150))
    u = UV[:, 1:50]
    v = UV[:, 51:100]
    w = UV[:, 101:150]
    return u, v, w
end

"""
Extracts Eulerian and bolus velocities u, v, w, Ub, Vb, and Wb from diagnostic files.

Parameters:
- `diagpath`: Dictionary containing paths to different experiments.
- `expname`: Name of the experiment.
- `fnameuvw`: File name of the diagnostic data.
- `γ`: Grid information.
- `Γ`: Grid information.
- `mskC`: Mask data for cell centers.
- `mskW`: Mask data for west cell interfaces.
- `mskS`: Mask data for south cell interfaces.

Returns:
- `u`: Array of u-component of Eulerian velocity.
- `v`: Array of v-component of Eulerian velocity.
- `w`: Array of w-component of Eulerian velocity.
- `Ub`: Array of u-component of bolus velocity.
- `Vb`: Array of v-component of bolus velocity.
- `Wb`: Array of w-component of bolus velocity.
"""
function extract_eulerian_and_bolus_velocities(
    diagpath, expname, fnameuvw, γ, Γ,
    mskC::MeshArrays.gcmarray{T, 2, Matrix{T}},
    mskW::MeshArrays.gcmarray{T, 2, Matrix{T}},
    mskS::MeshArrays.gcmarray{T, 2, Matrix{T}}
) where T<:Real
    UV = γ.read(diagpath[expname]*fnameuvw, MeshArray(γ, Float32, 250))
    u = UV[:, 1:50]
    v = UV[:, 51:100]
    w = UV[:, 101:150]
    GM_PsiX = UV[:, 151:200]
    GM_PsiY = UV[:, 201:250]
    Ub, Vb, Wb = calc_bolus(GM_PsiX, GM_PsiY, Γ, mskC, mskW, mskS)
    return u, v, w, Ub, Vb, Wb
end

"""
Extracts κUθ, κVθ, Uθ, and Vθ from diagnostic files for lateral heat budget.

Parameters:
- `diagpath`: Dictionary containing paths to different experiments.
- `expname`: Name of the experiment.
- `fnameHθ`: File name of the diagnostic data.
- `γ`: Grid information.

Returns:
- `κUθ`: U-component of diffusive heat flux
- `κVθ`: V-component of diffusive heat flux
- `Uθ`:  U-component of advective heat flux
        Uses Eulerian and Bolus velocities. 
- `Vθ`:  U-component of advective heat flux
        Uses Eulerian and Bolus velocities.
"""
function extract_lateral_heatbudget(diagpath, expname, fnameHθ, γ)
    dθλ = γ.read(diagpath[expname]*fnameHθ, MeshArray(γ, Float32, 200))
    κUθ = dθλ[:, 1:50]
    κVθ = dθλ[:, 51:100]
    Uθ = dθλ[:, 101:150]
    Vθ = dθλ[:, 151:200]
    return κUθ, κVθ, Uθ, Vθ
end

"""
Extracts κzθ and wθ from diagnostic files for vertical heat budget.

Parameters:
- `diagpath`: Dictionary containing paths to different experiments.
- `expname`: Name of the experiment.
- `fnameRθ`: File name of the diagnostic data.
- `γ`: Grid information.

Returns:
- `κzθ`: Array of diffusion in the "R" direction.
- `wθ`: Array of wθ component.
        Uses Eulerian and Bolus velocities

"""
function extract_vertical_heatbudget(diagpath, expname, fnameRθ, γ)
    # Vertical convergences
    dθr = γ.read(diagpath[expname]*fnameRθ, MeshArray(γ, Float32, 150))
    wθ = dθr[:, 1:50]
    κzθ = dθr[:, 51:100] .+ dθr[:, 101:150]
    return κzθ, wθ
end

"""
Calculates the bolus velocities Ub, Vb, and Wb from GM_PsiX, GM_PsiY, and grid information.

Parameters:
- `GM_PsiX`: X-component of GM streamfunction.
- `GM_PsiY`: Y-component of GM streamfunction.
- `Γ`: Grid information.
- `mskC`: Mask data for cell centers.
- `mskW`: Mask data for west cell interfaces.
- `mskS`: Mask data for south cell interfaces.

Returns:
- `Ub`: Array of u-component of bolus velocity.
- `Vb`: Array of v-component of bolus velocity.
- `Wb`: Array of w-component of bolus velocity.
"""
function calc_bolus(
    GM_PsiX::MeshArrays.gcmarray{T, 2, Matrix{T}},
    GM_PsiY::MeshArrays.gcmarray{T, 2, Matrix{T}},
    Γ,
    mskC::MeshArrays.gcmarray{T, 2, Matrix{T}},
    mskW::MeshArrays.gcmarray{T, 2, Matrix{T}},
    mskS::MeshArrays.gcmarray{T, 2, Matrix{T}}
) where {T<:Real}
    nr = length(Γ.RC)
    bolusU = T.(similar(Γ.hFacW))
    bolusV = T.(similar(Γ.hFacS))

    @inbounds for k = 1:nr-1
        bolusU[:, k] = (GM_PsiX[:, k+1] .- GM_PsiX[:, k]) ./ Γ.DRF[k]
        bolusV[:, k] = (GM_PsiY[:, k+1] .- GM_PsiY[:, k]) ./ Γ.DRF[k]
    end
    @inbounds bolusU.f[:, nr] = -GM_PsiX.f[:, nr] ./ Γ.DRF[nr]
    @inbounds bolusV.f[:, nr] = -GM_PsiY.f[:, nr] ./ Γ.DRF[nr]

    # And its vertical part
    # (seems correct, leading to 0 divergence)
    # tmpx and tmpy are the BOLUS transports
    tmp_x = T.(similar(GM_PsiX))
    tmp_y = T.(similar(GM_PsiX))
    bolusW = T.(similar(GM_PsiX))

    @inbounds for a in eachindex(tmp_x)
        tmp_x.f[a] .= GM_PsiX.f[a] .* Γ.DYG.f[a[1]]
        tmp_y.f[a] .= GM_PsiY.f[a] .* Γ.DXG.f[a[1]]
    end

    calc_UV_conv3D!(tmp_x, tmp_y, bolusW)
    @inbounds for a in eachindex(tmp_x)
        # Negative for divergence instead of convergence
        # Need to rescale bolus W
        bolusW.f[a] .= -1 .* bolusW.f[a] ./ Γ.RAC.f[a[1]]
        bolusU.f[a] .= bolusU.f[a] .* mskW.f[a]
        bolusV.f[a] .= bolusV.f[a] .* mskS.f[a]
        bolusW.f[a] .= bolusW.f[a] .* mskC.f[a]
    end

    return bolusU, bolusV, bolusW
end
