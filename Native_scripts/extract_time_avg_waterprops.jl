#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LaTeXStrings, PyCall, GibbsSeaWater 
import PyPlot as plt, NaNMath as nm
@pyimport seaborn as sns
@pyimport pandas as pd
# colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
colors =  sns.color_palette("deep")[1:4]
include(srcdir("config_exp.jl"))
cm = pyimport("cmocean.cm");colorway = cm.balance;
sns.set_theme(context = "poster", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = OHC_helper.wet_pts(Γ)
cell_depths = OHC_helper.get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);

runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

tecco= 1992+1/24:1/12:2018 # ecco years
nz = length(z)
nt = length(tecco)
P = MeshArray(γ,Float32, nz)
for ijk in eachindex(P)
    P[ijk] .= -pstdz[ijk[2]]
end

p₀ = 2000
ρθSavg = Dict()
for expname in ["iter129_bulkformula", "iter0_bulkformula", "seasonalclimatology"]
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    avg = Dict()

    start = 1
    ntt = length(tecco) - start + 1

    θz_mean = MeshArray(γ,Float32, nz); fill!(θz_mean, 0.0)
    σz_mean = MeshArray(γ,Float32, nz); fill!(σz_mean, 0.0)
    Sz_mean = MeshArray(γ,Float32, nz); fill!(Sz_mean, 0.0)

    for tt in start:nt
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        θname = datafilelist_θ[tt]

        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        @time θSz = γ.read(diagpath[expname]*θname,MeshArray(γ,Float32,2*nz))
        θz = θSz[:, 1:nz]; Sz = θSz[:, nz+1:end]
        for ijk in eachindex(P)
            σ = OHC_helper.densityJMD95.(θz.f[ijk],Sz.f[ijk], P[ijk], p₀) #EOS from MITGCM 
            σ = σ .- 1000
            θz_mean[ijk] .+= θz.f[ijk] ./ ntt
            Sz_mean[ijk] .+= Sz.f[ijk] ./ ntt
            σz_mean[ijk] .+= σ ./ ntt
        end 

    end

    avg["θ"] = θz_mean
    avg["σ2"] = σz_mean
    avg["S"] = Sz_mean

    ρθSavg[expname] = avg
end

svename = datadir("native/native_sigma2_timeavg_1992_2017.jld2")
jldsave(svename, ρθSavg = ρθSavg)