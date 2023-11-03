include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
area32 = Float32.(area)

ocean_mask = wet_pts(Γ)
region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
H = vertical_sum(cell_depths); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H
cell_volumes = get_cell_volumes(area, cell_depths);

runpath,diagpath = listexperiments(exprootdir());

include(srcdir("plot_and_dir_config.jl"))

ϕ_min_mask, ϕ_max_mask = get_ϕ_max_min_mask(region, Γ, λ, ϕ, basins, basin_list)
mskC, mskW, mskS = get_msk(Γ)

#define control volume 
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= -z[:].<= uplvl)
ctrl_vol = sum(cell_volumes[:, lvls])

function interpolate_to_vertical_faces!(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    ma_interp::MeshArrays.gcmarray{T, 2, Matrix{T}})  where T<:Real
    for k = 1:49 
        ma_interp[:, k+1] = (ma[:, k] .+ ma[:, k+1]) ./2 #linear interpolate to faces
    end
end


function filter_heat_budget(diagpath::Dict{String, String}, expname::String, γ::gcmgrid)

    filelist = searchdir(diagpath[expname], "state_2d_set1") # first filter for state_3d_set1
    datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"trsp_3d_set2")

    sθinterp= MeshArray(γ,Float32,50); 
    sθinterp0= MeshArray(γ,Float32,50); 

    nt = 312;
    var_names = ["wθTop", "wθTop_bol", "wθBot", "wθBot_bol", 
                 "ΔW_θTop", "ΔW_θBot", "W_ΔθTop", "W_ΔθBot", "ΔW_ΔθTop", "ΔW_ΔθBot"]
    vars = Dict(varname => zeros(Float32, nt) for varname in var_names)
    println(expname)
    @time for tt = 1:nt
        println(tt)
        fnameS = datafilelist_S[tt]
        fnameθ = datafilelist_θ[tt]
        sθ = extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_H)
        sθ_0 = extract_sθ("iter0_bulkformula",diagpath, γ, fnameS, fnameθ, inv_H)

        θ_volavg = sum(sθ[:, lvls] .* cell_volumes[:, lvls]) / ctrl_vol

        θ_volavg0 = sum(sθ_0[:, lvls] .* cell_volumes[:, lvls]) / ctrl_vol

        fill!(sθinterp, 0.0)
        fill!(sθinterp0, 0.0)

        interpolate_to_vertical_faces!(sθ, sθinterp) #interpolate temperatures
        interpolate_to_vertical_faces!(sθ_0, sθinterp0) #interpolate temperatures

        sθinterp = sθinterp .- θ_volavg
        sθinterp0 = sθinterp0 .- θ_volavg0

        _, _, w, _, _, Wb = extract_eulerian_and_bolus_velocities(diagpath, expname, datafilelist_uvw[tt], 
                                                                      γ, Γ, mskC, mskW, mskS)

        _, _, w0, _, _, Wb0 = extract_eulerian_and_bolus_velocities(diagpath, "iter0_bulkformula", datafilelist_uvw[tt], 
                                                                      γ, Γ, mskC, mskW, mskS)
        for i in eachindex(w)
          w.f[i] .= w.f[i] .+ Wb.f[i]
          w0.f[i] .= w0.f[i] .+ Wb0.f[i]

        end  

        W = w .* area32; W_bol = Wb .* area32
        W0 = w0 .* area32; W_bol0 = Wb0 .* area32
        
        ΔW = W .- W0; ΔWbol = W_bol .- W_bol0; Δθ = sθinterp .- sθinterp0

        ΔW_θ = ΔW .* sθinterp
        W_Δθ = W .* Δθ
        ΔW_Δθ = ΔW .* Δθ

        Wθ = W .* sθinterp
        Wθ_bol = W_bol .* sθinterp

        vars["wθTop"][tt]     = sum(Wθ[:, lvls[1]] .* PAC_msk )/ ctrl_vol
        vars["wθBot"][tt]     = sum(Wθ[:, lvls[end] + 1] .* PAC_msk )/ ctrl_vol

        vars["ΔW_θTop"][tt]     = sum(ΔW_θ[:, lvls[1]] .* PAC_msk )/ ctrl_vol
        vars["ΔW_θBot"][tt]     = sum(ΔW_θ[:, lvls[end] + 1] .* PAC_msk )/ ctrl_vol

        vars["W_ΔθTop"][tt]     = sum(W_Δθ[:, lvls[1]] .* PAC_msk )/ ctrl_vol
        vars["W_ΔθBot"][tt]     = sum(W_Δθ[:, lvls[end] + 1] .* PAC_msk )/ ctrl_vol

        vars["ΔW_ΔθTop"][tt]     = sum(ΔW_Δθ[:, lvls[1]] .* PAC_msk )/ ctrl_vol
        vars["ΔW_ΔθBot"][tt]     = sum(ΔW_Δθ[:, lvls[end] + 1] .* PAC_msk )/ ctrl_vol

    
        vars["wθTop_bol"][tt] = sum(Wθ_bol[:, lvls[1]] .* PAC_msk )/ ctrl_vol
        vars["wθBot_bol"][tt] = sum(Wθ_bol[:, lvls[end] + 1] .* PAC_msk )/ ctrl_vol
        
    end
    return vars
end

vars =  ["only_init", "only_kappa", "iter129_bulkformula",  "iter0_bulkformula", "only_buoyancy", "only_wind"]
for expt in vars
    println(expt)
    vars = filter_heat_budget(diagpath, expt, γ)
    savename = datadir("native/" * expt * region * "_THETA_budget_Wθ_eul_bol_ΔDecomp" * suffix * ".jld2")
    jldsave(savename, dθ = vars)
end
