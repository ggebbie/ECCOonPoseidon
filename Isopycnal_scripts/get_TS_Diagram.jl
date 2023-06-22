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

include(srcdir("config_exp.jl"))
include("./IsopycnalHelpers.jl")

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "NPAC"
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);

#computed on standard pressures, 
#may have a different value if using pressures diagnosed from ECCO
#https://ecco-v4-python-tutorial.readthedocs.io/VectorCalculus_ECCO_barotropicVorticity.html
runpath,diagpath = listexperiments(exprootdir());

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)

#read in the first time step of S and θ
expname = "iter129_bulkformula"

sig1grid = sigma1grid()
nσ = length(sig1grid)
# σlvls = findall( σtop .<= sig1grid .<= σbot)

#get root for salinity on sigma1
#THETA for iter129
Troot = "THETA_on_sigma1" 
Sroot = "SALT_on_sigma1" 

filelist = searchdir(runpath[expname]*"sigma1/",Troot) # first filter for state_3d_set1
datafilelistT  = filter(x -> occursin("data",x),filelist) # second filter for "data"

filelist = searchdir(runpath[expname]*"sigma1/",Sroot) # first filter for state_3d_set1
datafilelistS  = filter(x -> occursin("data",x),filelist) # second filter for "data"

nt = length(datafilelistT)

fname = expname * region * "_AVG_P_sigma1.jld2" 
!ispath(σ1datadir(fname)) && (include("get_avgP_on_sigma1.jl"))
Pσ = load(σ1datadir(fname), "P")
sigma2p = mean(Pσ, dims = 2)
p2z = gsw_z_from_p.(abs.(sigma2p), 30., 0., 0.)[:]
plvls = findall( -1500 .<= p2z .<= -1000)
iσ = 66

tstep1 = Dict(v => Float32[] for v in ["θ", "S"])
tstepend = Dict(v => Float32[] for v in ["θ", "S"])
sts = [1, 301]; ds = [tstep1, tstepend]

for (st, d) in zip(sts, ds)
    θiσ1 = MeshArray(γ,Float32); fill!(θiσ1, 0.0)
    Siσ1 = MeshArray(γ,Float32); fill!(Siσ1, 0.0)
    tts =  st .+ collect(0:11)
    for tt in tts
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        Tname = datafilelistT[tt]
        Sname = datafilelistS[tt]

        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        @time θσ1 = γ.read(runpath[expname]*"sigma1/"*Tname,MeshArray(γ,Float32,nσ))
        @time Sσ1 = γ.read(runpath[expname]*"sigma1/"*Sname,MeshArray(γ,Float32,nσ))
        msk_θ = θσ1[:, iσ]
        msk_S = Sσ1[:, iσ]
        [msk_θ[ff][ PAC_msk[ff] .== 0] .= NaN for ff = 1:5]
        [msk_S[ff][ PAC_msk[ff] .== 0] .= NaN for ff = 1:5]

        θiσ1 .+= msk_θ ./ 12
        Siσ1 .+= msk_S ./ 12
    end
        #filter for nans
    for ff=1:5
        push!(d["θ"], filter(!isnan, θiσ1[ff])... )
        push!(d["S"], filter(!isnan, Siσ1[ff])...)
    end

end
using Unitful
σ_sel = string(sig1grid[iσ])

import PyPlot as plt
fig, axes = plt.subplots(ncols=2, sharex = true, sharey = true)

axes[1].scatter(tstep1["θ"] , tstep1["S"], label = "1992", c= "black", s = 2)
axes[2].scatter(tstepend["θ"] , tstepend["S"], label = "2017", c="black", s=2, marker = "*")

[ax.set_xlabel("Temperature [°C]") for ax in axes]
[ax.legend() for ax in axes]

axes[1].set_ylabel("Salinity [psu]")
fig.suptitle("TS Distribution for σ = " * σ_sel * "\n for " * expname)
fig.tight_layout()
fig
