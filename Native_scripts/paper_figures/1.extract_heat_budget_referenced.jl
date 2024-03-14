include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
import NaNMath as nm
import PyPlot as plt
import NumericalIntegration

cumul_integrate = NumericalIntegration.cumul_integrate
include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");

(ϕ,λ) = latlonC(γ)
area = readarea(γ)


ocean_mask = wet_pts(Γ)
region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)

H = vertical_sum(cell_depths); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H

include(srcdir("plot_and_dir_config.jl"))

ϕ_min_mask, ϕ_max_mask = get_ϕ_max_min_mask(region, Γ, λ, ϕ, basins, basin_list)

#define control volume 
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= -z[:].<= uplvl)

function filter_heat_budget(diagpath::Dict{String, String}, expname::String, γ::gcmgrid, 
    cell_volumes, H, lvls, mask)
    filelist = searchdir(diagpath[expname], "state_2d_set1") # first filter for state_3d_set1
    datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"trsp_3d_set2")
    datafilelist_H  = filter(x -> occursin("data",x),filelist) 
    filelist = searchdir(diagpath[expname],"trsp_3d_set3")
    datafilelist_R  = filter(x -> occursin("data",x),filelist)

    cs, sn = get_cs_and_sn(γ)
    GTF3d= get_geothermalheating(Γ, γ)

    nt = 312;
    var_names = ["θ", "κxyθ", "κzθ", "GTH", "VθSouth", "VθNorth", "wθTop", "wθBot", "wθ", "uvθ"]
    vars = Dict(varname => zeros(Float32, nt) for varname in var_names)

    κxyθ_conv3D = MeshArray(γ,Float32,50);
    κzθ_conv3D = MeshArray(γ,Float32,50);
    ctrl_vol = sum(cell_volumes[:, lvls])

    uvθ_conv3D = MeshArray(γ,Float32,50);
    wθ_conv3D = MeshArray(γ,Float32,50);

    mskC, mskW, mskS = get_msk(Γ)

    area32 = Float32.(area)

    vars["GTH"] .= sum(GTF3d[:, lvls] .* cell_volumes[:, lvls]) / ctrl_vol
    θ_V = MeshArray(γ,Float32,50); fill!(θ_V, 0.0)
    θ_U = MeshArray(γ,Float32,50); fill!(θ_U, 0.0)
    θ_W = MeshArray(γ,Float32,50); fill!(θ_W, 0.0)

    for tt = 1:nt
        println(tt)
        fnameS = datafilelist_S[tt]
        fnameθ = datafilelist_θ[tt]
        sθ = extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_H)
        vars["θ"][tt] = sum(sθ[:, lvls] .* cell_volumes[:, lvls]) / ctrl_vol
        sθ_ref = sθ .- vars["θ"][tt]
        κUθ, κVθ, Uθ, Vθ = extract_lateral_heatbudget(diagpath, expname , datafilelist_H[tt], γ)
        κzθ, wθ = extract_vertical_heatbudget(diagpath, expname , datafilelist_R[tt], γ)
        u, v, w, Ub, Vb, Wb = extract_eulerian_and_bolus_velocities(diagpath, expname, 
        datafilelist_uvw[tt], γ, Γ, mskC, mskW, mskS)

        u .+= Ub
        v .+= Vb
        w .+= Wb
        
        Utrsp, Vtrsp = UVtoTrsp(u, v, Γ, γ)
        Wtrsp = w .* area32


        Uθ_ref = Utrsp .* vars["θ"][tt]
        Vθ_ref = Vtrsp .* vars["θ"][tt]
        Wθ_ref = Wtrsp .* vars["θ"][tt]

        # Gael's function produces numerical error! Can't use for sensitive budgets!
        # Eθ, Nθ = rotate_uv(Uθ .- Uθ_ref, Vθ .- Vθ_ref, Γ) #should rotate things
        Eθ, Nθ = rotate_UV_native(Uθ .- Uθ_ref, Vθ .- Vθ_ref, cs, sn) 
        #doesnt seem to work, need to try interpolating to cells 

        wθref_top = wθ[:, lvls[1]] .- Wθ_ref[:, lvls[1]] 
        wθref_bot = wθ[:, lvls[end]+1] .- Wθ_ref[:, lvls[end]+1] 

        calc_UV_conv3D!(κUθ, κVθ, κxyθ_conv3D);
        calc_W_conv3D!(κzθ, κzθ_conv3D)

        calc_UV_conv3D!(Uθ, Vθ, uvθ_conv3D);
        calc_W_conv3D!(wθ, wθ_conv3D)

        vars["VθSouth"][tt] = sum(Nθ[:, lvls] .* ϕ_min_mask ) / ctrl_vol
        vars["VθNorth"][tt] = sum(Nθ[:, lvls] .* ϕ_max_mask ) / ctrl_vol
        vars["wθTop"][tt] = sum(wθref_top  .* mask )/ ctrl_vol
        vars["wθBot"][tt] = sum(wθref_bot .* mask )/ ctrl_vol
        
        vars["wθ"][tt] = sum(uvθ_conv3D[:, lvls] .* mask ) / ctrl_vol
        vars["uvθ"][tt] = sum(wθ_conv3D[:, lvls] .* mask ) / ctrl_vol

        vars["κxyθ"][tt] = sum(κxyθ_conv3D[:, lvls] .* mask ) / ctrl_vol
        vars["κzθ"][tt] = sum(κzθ_conv3D[:, lvls] .* mask ) / ctrl_vol

    end
    return vars
end

# vars = ["iter129_bulkformula",  "iter0_bulkformula", "seasonalclimatology"]
exps =  ["only_init", "only_kappa", "only_sfc", "iter129_bulkformula",  "iter0_bulkformula"]
exps =  ["only_buoyancy", "only_wind"]

exps =  [ "only_init", "only_wind", "only_buoyancy", "only_kappa"]

for expname in exps
    vars = filter_heat_budget(diagpath, expname, γ, cell_volumes, H, lvls, PAC_msk)
    savename = datadir("native/" * expname * region * "_THETA_budget_ref_with_Bolus_wextra" * suffix * ".jld2")
    jldsave(savename, dθ = vars)
end
