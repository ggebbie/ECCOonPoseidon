#this script uses T-S data to create Hovmoller diagrams at a single point 

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings,
    PyCall
using .OHC_helper
import PyPlot as plt
@pyimport seaborn as sns;
@pyimport pandas as pd;
sns.set_theme(context = "talk", style = "ticks",
              palette = sns.color_palette("colorblind"));
              
cm = pyimport("cmocean.cm");


include(srcdir("config_exp.jl"))

#using . is faster 
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());

# abbreviations for each experiment for labels, etc.
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
nexp = length(shortnames) # number of experiments
 
ocean_mask = OHC_helper.wet_pts(Γ)
region = "Sargasso"
msk = region_mask(ocean_mask, λ, ϕ, (22, 37), (-73, -60))
nz = 50

cell_depths = OHC_helper.get_cell_depths(msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));

H = OHC_helper.smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H
tecco = 1992+1/24:1/12:2018

monthsperyear = 12
fcycle = 1 # units: yr^{-1}
# for removing seasonal cycle from monthly averaged timeseries
overtones= 4; # semi-annual, quad-annual, etc.: helps capture an asymmetric seasonal cycle
Ecycle,Fcycle = seasonal_matrices(fcycle,tecco,overtones)

#pre-allocate

# vm = maximum(abs.(θ["iter129_bulkformula"] - mean(θ["iter129_bulkformula"], dims = 2)))
savename = datadir("native/iter129_bulkformula" * region * "_WaterProp_budget_z.jld2")
θS = load(savename)["θS"]
kmax = 22; z[k]
θ = θS["θ"][1:kmax, :]; S = θS["S"][1:kmax, :]
σ1 = OHC_helper.densityJMD95.(θ,S, NaN, 1000) .- 1000. #EOS from MITGCM 

[θ[i, :] .= remove_seasonal(θ[i, :][:],Ecycle,Fcycle) for i=1:kmax]
[S[i, :] .= remove_seasonal(S[i, :][:],Ecycle,Fcycle) for i=1:kmax]
[σ1[i, :] .= remove_seasonal(σ1[i, :][:],Ecycle,Fcycle) for i=1:kmax]

# θ = θ .- mean(θ, dims = 2)
# S = S .- mean(S, dims = 2)
# σ1 = σ1 .- mean(σ1, dims = 2)

fig, ax = plt.subplots(3, 1, figsize = (15, 15))
v = 0.5; levels =collect(-v:0.1:v)
CM = ax[1].contourf(tecco, z[1:kmax], θ, cmap = cm.balance, vmin = -v, vmax = v, levels = levels, extend = "both"); 
CS = ax[1].contour(tecco, z[1:kmax], θ, vmin = -v, vmax = v, levels = levels, colors = "k")
# ax[1].clabel(CS, inline=1, fontsize=20, inline_spacing = 10)
fig.colorbar(CM, ax = ax[1], orientation = "vertical", fraction = 0.05, extend = "both")
ax[1].set_title("Temperature Anomaly")

v = 0.1; levels =collect(-v:0.025:v)
CM = ax[2].contourf(tecco, z[1:kmax], S, cmap = cm.delta, vmin = -v, vmax = v, levels = levels, extend = "both");
CS = ax[2].contour(tecco, z[1:kmax], S, vmin = -v, vmax = v, levels = levels, colors = "k")
# ax[2].clabel(CS, inline=1, fontsize=20, inline_spacing = 10, fmt = "%1.2f")
fig.colorbar(CM, ax = ax[2], orientation = "vertical", fraction = 0.05, extend = "both")
ax[2].set_title("Salinity Anomaly")

v = 0.5; levels =collect(-v:0.05:v)
CM = ax[3].contourf(tecco, z[1:kmax], σ1, cmap = cm.curl, vmin = -v, vmax = v, levels = levels, extend = "both");
CS = ax[3].contour(tecco, z[1:kmax], σ1, vmin = -v, vmax = v, levels = levels, colors = "k")
fig.colorbar(CM, ax = ax[3], orientation = "vertical", fraction = 0.05, extend = "both")
# ax[3].clabel(CS, inline=1, fontsize=20, inline_spacing = 10, fmt = "%1.2f")
ax[3].set_title("Salinity Anomaly")


[a.set_xlabel("time") for a in ax]
[a.set_ylabel("depth") for a in ax]
fig.tight_layout()
fig