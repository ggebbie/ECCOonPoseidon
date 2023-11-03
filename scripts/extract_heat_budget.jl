
include("../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings

include(srcdir("config_exp.jl"))

area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());

ocean_mask = wet_pts(Γ)

cell_depths = get_cell_thickness(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths);

H = OHC_helper.sum_vertical(cell_depths, γ); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H

#get the geothermal heating term
GTF = get_geothermalheating(Γ, γ)

#define experimennt
expname = "iter129_bulkformula"
get_data_files(xx, expname) = filter(x -> occursin("data",x),searchdir(diagpath[expname],xx))  
#define the file names

datafilelist_H  = get_data_files("trsp_3d_set2", expname)
datafilelist_R  = get_data_files("trsp_3d_set3", expname)
datafilelist_θ  = get_data_files("state_3d_set1", expname)
datafilelist_S  = get_data_files("state_2d_set1", expname)

#specify time-step
tt = 1

fnameθ = datafilelist_θ[tt]
fnameH = datafilelist_H[tt]
fnameR = datafilelist_R[tt]
fnameS = datafilelist_S[tt]

#Load in temperature and heat budget terms
sθ = OHC_helper.extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_H)

κUθ, κVθ, Uθ, Vθ = extract_lateral_heatbudget(diagpath, expname , datafilelist_H[tt], γ)
κzθ, wθ = extract_vertical_heatbudget(diagpath, expname , datafilelist_R[tt], γ)

#compute convergences for all ocean gridcells
κθ_conv3D = MeshArray(γ,Float32,50);
calc_UV_conv3D!(κUθ, κVθ, κθ_conv3D); #horizontal diffusive converegence

uθ_conv3D = MeshArray(γ,Float32,50);
calc_UV_conv3D!(Uθ, Vθ, uθ_conv3D); #horizontal advective converegence

wθ_conv = MeshArray(γ,Float32,50);
calc_W_conv3D!(wθ, wθ_conv) #vertical diffusive converegence

κzθ_conv = MeshArray(γ,Float32,50);
calc_W_conv3D!(κzθ, κzθ_conv) #vertical advective converegence

#dθ/dt = κθ_conv3D + uθ_conv3D + wθ_conv + κzθ_conv + GTF