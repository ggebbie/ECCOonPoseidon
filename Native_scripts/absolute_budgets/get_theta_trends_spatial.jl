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
colors =  sns.color_palette("colorblind")[1:4]
labels = ["ECCO", "ECCO Forcing", "ECCO Init. and κ", "CTRL"]
sns.set_theme(context = "poster", style = "ticks", font_scale = 1.2,
              palette = sns.color_palette("colorblind"));
pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;
include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
tecco = 1992+1/24:1/12:2018
runpath,diagpath = listexperiments(exprootdir());

ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

#load in latitude mask 
ocean_mask = wet_pts(Γ)
suffix = "2to3"
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)

#create volume mask
area = readarea(γ)
cell_depths = get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
H = smush(cell_depths[:, lvls]); 
H[findall(H .==0 )] = Inf
inv_depths = 1 ./ H

#load trend trend_matricesß, F is LS estimator
E,F = trend_matrices(tecco)

cf = Vector{Any}(undef ,1)

β = MeshArray(γ,Float32); fill!(β, 0.0)

expname = "nosfcadjust"
filelist = searchdir(diagpath[expname], "state_2d_set1") # first filter for state_3d_set1
datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
nt = length(datafilelist_S)

@time for tt = 1:nt
    fnameS = datafilelist_S[tt]
    fnameθ = datafilelist_θ[tt]
    sθ = extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_depths)
    sθH = depth_average(sθ[:, lvls], cell_depths[:, lvls], H, γ)
    for ff in 1:5
        β[ff] .+= F[2,tt] .* sθH[ff] 
    end
end

β_reg = LLCcropC(β,γ) 
regular_λ = LLCcropC(λ,γ)
regular_ϕ = LLCcropC(ϕ,γ)

fname = expname * "_θ_trends_2to3km.jld2"
jldsave(datadir(fname), β = β_reg, λ = regular_λ, ϕ = regular_ϕ)

