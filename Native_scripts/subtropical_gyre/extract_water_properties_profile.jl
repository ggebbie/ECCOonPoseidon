#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
# include("../../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings,
    PyCall
# using .OHC_helper
import PyPlot as plt
@pyimport seaborn as sns;
@pyimport pandas as pd;
sns.set_theme(context = "talk", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
              
cm = pyimport("cmocean.cm");


include(srcdir("config_exp.jl"))

#using . is faster 
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
 
include("SargassoMask.jl")

cell_depths = get_cell_thickness(msk, ΔzF, Γ.hFacC); 
cell_volumes = (get_cell_volumes(area, cell_depths));

# H = OHC_helper.smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
# inv_H = 1 ./ H
tecco = 1992+1/24:1/12:2018; nz = 50

fcycle = 1 # units: yr^{-1}
# for removing seasonal cycle from monthly averaged timeseries
overtones= 4; # semi-annual, quad-annual, etc.: helps capture an asymmetric seasonal cycle
Ecycle,Fcycle = seasonal_matrices(fcycle,tecco,overtones)


function filter_waterprops_budget(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    filelist = searchdir(diagpath[expname],"state_3d_set1") 
    datafilelist_θ  = filter(x -> occursin("data",x),filelist)
 
    var_names = ["θ", "S"]
    nt = length(datafilelist_θ); nz = 50
    vars = Dict(varname => zeros(Float32, nz, nt) for varname in var_names)
    face, index, _ = findlatlon(λ, ϕ, -80, 26);

    @time for tt = 1:312
        println(tt)
        fnameθ = datafilelist_θ[tt]
        θS = γ.read(diagpath[expname]*fnameθ,MeshArray(γ,Float32,100))

        #horizontal convergences
        θ = θS[:, 1:50]
        S = θS[:, 51:100]
        for k = 1:50
            vars["θ"][k, tt]  = θ[face, k][index]
            vars["S"][k, tt]  = S[face, k][index]
        end

    end

    return vars
end

for expname in ["iter129_bulkformula"]
    savename = datadir("native/" * expname * "FtLd_WaterProp_profile_budget_z.jld2")
    θS = filter_waterprops_budget(diagpath, expname, γ)
    jldsave(savename, θS = θS)
end


fig, ax = plt.subplots(2, 1, figsize = (15, 10))
region = "FtLd"
savename = datadir("native/iter129_bulkformula" * region * "_WaterProp_profile_budget_z.jld2")
θS = load(savename)["θS"]
kmax = 20;
θ = θS["θ"][1:kmax, :]; S = θS["S"][1:kmax, :]
θ[θ .== 0] .= NaN
S[S .== 0] .= NaN

θmean = mean(θ, dims = 2); Smean = mean(S, dims = 2)

[θ[i, :] .= remove_seasonal(θ[i, :][:],Ecycle,Fcycle) for i=1:kmax]
[S[i, :] .= remove_seasonal(S[i, :][:],Ecycle,Fcycle) for i=1:kmax]

v = 0.5; levels =collect(-v:0.1:v)

CB = ax[1].pcolormesh(tecco, z[1:kmax], θ, cmap = cm.balance, vmin = -v, vmax = v); 
fig.colorbar(CB, ax = ax[1], label = "K")
ax[1].set_title("80W, 26N;Temperature Anomaly (Seasonal Cycle Removed)")

v = 0.1; levels =collect(-v:0.025:v)
CB = ax[2].pcolormesh(tecco, z[1:kmax], S, cmap = cm.delta, vmin = -v, vmax = v);
fig.colorbar(CB, ax = ax[2], label = "PSU")
# ax[2].clabel(CS, inline=1, fontsize=17, inline_spacing = 10)
ax[2].set_title("80W, 26N;Salinity Anomaly (Seasonal Cycle Removed)")
[a.set_xlabel("time") for a in ax]
[a.set_ylabel("depth") for a in ax]
fig.tight_layout()
ax[1].invert_yaxis()
ax[2].invert_yaxis()
fig
fig.savefig(plotsdir("native/subtropical_gyre/HovmollerFtLtd.png"))
fig