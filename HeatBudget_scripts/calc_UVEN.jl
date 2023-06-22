#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings,
Plots
gr()
import NaNMath as nm
using .OHC_helper
import PyPlot as plt
using PyCall

@pyimport seaborn as sns;
@pyimport pandas as pd;
colors =  sns.color_palette("deep")[1:4]
sns.set_theme(context = "talk", font_scale = 1.0,
              palette = sns.color_palette("deep"));#sns.set_context("talk")
pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;
include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
area = readarea(γ)
runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
labels_L = [L"\theta^{129}", L"\theta^{\Delta F}",L"\theta^{\Delta T}", L"\theta^{0}"]
tecco = 1992+1/24:1/12:2018

#load in velocities 
expname = "iter129_bulkformula"
filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"
tt = 20 * 13
fname_uvw = datafilelist_uvw[tt]
UVW = γ.read(diagpath[expname]*fname_uvw,MeshArray(γ,Float32,150))
ut = UVW[:, 1:50];
vt = UVW[:, 51:100];
u, v = UVtoTransport(ut, vt, Γ)
#load in latitude mask 
ocean_mask = wet_pts(Γ)
region = "NPAC30"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "full")
suffix = "2to3"
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)
ϕ_mask_Inf = ϕ .* PAC_msk; ϕ_mask_Inf[findall(ϕ_mask_Inf .== 0 )] .= Inf
ϕ_mask_min = minimum(ϕ_mask_Inf) 
#this should actually correspond to the latitude just before 
ϕ_S = (ϕ .≈ ϕ_mask_min)
ϕ_S = ϕ_S .* PAC_msk


u1 = deepcopy(u); v1 = deepcopy(v); 
E1 = MeshArray(γ,Float32,50); N1 = MeshArray(γ,Float32,50)
@time for k in 1:50
    Etr_tmp, Ntr_tmp = UVtoUEVN(u1[:, k],v1[:, k],Γ)
    E1.f[:, k] .= Etr_tmp.f; N1.f[:, k] .= Ntr_tmp.f
end

u2 = deepcopy(u); v2 = deepcopy(v); 
@time E2, N2 = UVtoUEVN3D(u2,v2,Γ);
u3 = deepcopy(u); v3 = deepcopy(v); 
@time E3, N3  = rotate_uv(u3,v3,Γ); #is wrong!!!


for ff in eachindex(u)
    N1.f[ff][isnan.(N1.f[ff])].=0
    N2.f[ff][isnan.(N2.f[ff])].=0
    N3.f[ff][isnan.(N3.f[ff])].=0
end

function at_lat(X, ϕ, ϕi)
    nz = size(X, 2)
    temp = [0.0]
    for ff ∈ 1:5, k ∈ 1:nz
        temp .+= sum(X.f[ff, k][ϕ[ff] .≈ ϕi])
    end
    return temp 
end
# println("N1[2, 20] ", nm.sum.(N1.f[:, 37]), "\n", 
# " N2[2, 20]  ",nm.sum.(N2.f[:, 37]), "\n", " N3[1, 20]", "\n", nm.sum.(N3.f[:, 37]))

println(at_lat(N1[:, lvls].*PAC_msk, ϕ, ϕ_mask_min))
println(at_lat(N2[:, lvls].*PAC_msk, ϕ, ϕ_mask_min))
println(at_lat(N3[:, lvls].*PAC_msk, ϕ, ϕ_mask_min))