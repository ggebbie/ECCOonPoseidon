#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
using .OHC_helper
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
cell_depths = OHC_helper.get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
H = OHC_helper.smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area .* PAC_msk, cell_depths));

uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= z[:].<= uplvl)

runpath,diagpath = listexperiments(exprootdir());

ϕ_min_mask, ϕ_max_mask = OHC_helper.get_ϕ_max_min_mask(region, Γ, λ, ϕ, basins, basin_list)
sum(ϕ_min_mask)
#define the control volume of interest 
mskC, mskW, mskS = OHC_helper.get_msk(Γ)
cs, sn = OHC_helper.get_cs_and_sn(γ)
GTF3d= OHC_helper.get_geothermalheating(Γ, γ)

function filter_heat_budget(diagpath::Dict{String, String}, expname::String, γ::gcmgrid)
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

    nt = 100; ctrl_vol = sum(cell_volumes[:, lvls])
    var_names = ["θ", "κxyθ", "GTH", "VθSouth", "VθNorth", "wθTop", "wθBot", "κzθ", "VθConv", "WθConv"]
    vars = Dict(varname => zeros(Float32, nt) for varname in var_names)
    vars["GTH"] .= sum(GTF3d[:, lvls] .* cell_volumes[:, lvls]) / ctrl_vol
    uxyθ_conv3D = MeshArray(γ,Float32,50);
    κxyθ_conv3D = MeshArray(γ,Float32,50);
    κzθ_conv3D = MeshArray(γ,Float32,50);

    area32 = Float32.(area)
    @time for tt = 1:nt
        println(tt)
        fnameS = datafilelist_S[tt]
        fnameθ = datafilelist_θ[tt]
        sθ = OHC_helper.extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_H)
        vars["θ"][tt] = sum(sθ[:, lvls] .* cell_volumes[:, lvls]) / ctrl_vol

        κUθ, κVθ, Uθ, Vθ = OHC_helper.extract_heatbudgetH(diagpath, expname , datafilelist_H[tt], γ)
        κzθ, wθ = OHC_helper.extract_heatbudgetR(diagpath, expname , datafilelist_R[tt], γ)
        u, v, w, Ub, Vb, Wb = OHC_helper.extract_velocities_and_bolus(diagpath, expname , 
        datafilelist_uvw[tt], γ, Γ, mskC, mskW, mskS)
        for i in eachindex(u)
          u.f[i] .= u.f[i] .+ Ub.f[i]
          v.f[i] .= v.f[i] .+ Vb.f[i]
          w.f[i] .= w.f[i] .+ Wb.f[i]
        end

        Utrsp, Vtrsp = OHC_helper.UVtoTrsp(u, v, Γ)
        Wtrsp = w .* area32
        Uθ_ref = Utrsp .* vars["θ"][tt];
        Vθ_ref = Vtrsp .* vars["θ"][tt];
        Wθ_ref = Wtrsp .* vars["θ"][tt];
        #this does change what happens at the faces, but not the 3D convergence! 
        # Eθ, Nθ = rotate_uv(Uθ .- Uθ_ref, Vθ .- Vθ_ref, Γ) #should rotate things
        Eθ, Nθ = OHC_helper.rotate_UV_native(Uθ .- Uθ_ref, Vθ .- Vθ_ref, cs, sn) 

        OHC_helper.calc_UV_conv3D!(Uθ .- Uθ_ref, Vθ .- Vθ_ref, uxyθ_conv3D);
        OHC_helper.calc_UV_conv3D!(κUθ, κVθ, κxyθ_conv3D);
        OHC_helper.calc_Wconv3D!(κzθ, κzθ_conv3D)
        wθref_top = wθ[:, lvls[1]] .- Wθ_ref[:, lvls[1]] 
        wθref_bot = wθ[:, lvls[end]+1] .- Wθ_ref[:, lvls[end]+1] 
        wθ_conv = (wθ[:, lvls[end]+1] .- wθ[:, lvls[1]]) .- (Wθ_ref[:, lvls[end]+1] .- Wθ_ref[:, lvls[1]])
        vars["VθSouth"][tt] = sum(Nθ[:, lvls] .* ϕ_min_mask ) / ctrl_vol
        vars["VθNorth"][tt] = sum(Nθ[:, lvls] .* ϕ_max_mask ) / ctrl_vol
        vars["VθConv"][tt] = sum(uxyθ_conv3D[:, lvls] .* PAC_msk ) / ctrl_vol

        vars["wθTop"][tt] = sum(wθref_top  .* PAC_msk )/ ctrl_vol
        vars["wθBot"][tt] = sum(wθref_bot .* PAC_msk )/ ctrl_vol
        vars["WθConv"][tt] = sum(wθ_conv .* PAC_msk ) / ctrl_vol

        vars["κxyθ"][tt] = sum(κxyθ_conv3D[:, lvls] .* PAC_msk ) / ctrl_vol
        vars["κzθ"][tt] = sum(κzθ_conv3D[:, lvls] .* PAC_msk ) / ctrl_vol

    end
    return vars
end

for expname in ["iter129_bulkformula"]
    vars = filter_heat_budget(diagpath, expname, γ)
    savename = datadir("native/" * expname * region * "_THETA_budget_ref_TEST" * suffix * ".jld2")
    jldsave(savename, dθ = vars)
end
