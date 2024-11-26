#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

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
sns.set_theme(context = "talk", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
              
cm = pyimport("cmocean.cm");


include(srcdir("config_exp.jl"))

#using . is faster 
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
 
include("SargassoMask.jl")

cell_depths = OHC_helper.get_cell_depths(msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));

H = OHC_helper.smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H
tecco = 1992+1/24:1/12:2018
#pre-allocate

ΔV = zeros(Float32, nz)
for k=1:nz
    ΔV[k] = Float32(sum(cell_volumes[:, k]))
end


function heat_flux_profile(ds::MeshArray, ΔV)
    nz = size(ds, 2)
    vol_avg = zeros(Float32, nz)

    for ff=1:5, k=1:nz
        vol_avg[k] += Float32(sum(ds[ff, k])) ./ ΔV[k]
    end
    return vol_avg
end


function filter_waterprops_budget(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    filelist = searchdir(diagpath[expname],"state_3d_set1") 
    datafilelist_θ  = filter(x -> occursin("data",x),filelist)
 
    var_names = ["θ", "S"]
    nt = length(datafilelist_θ); nz = 50
    vars = Dict(varname => zeros(Float32, nz, nt) for varname in var_names)

    @time for tt = 1:312
        println(tt)
        fnameθ = datafilelist_θ[tt]
        θS = γ.read(diagpath[expname]*fnameθ,MeshArray(γ,Float32,100))

        #horizontal convergences
        θ = θS[:, 1:50]
        S = θS[:, 51:100]

        vars["θ"][:, tt]  = heat_flux_profile(θ .* cell_volumes, ΔV); 
        vars["S"][:, tt]  = heat_flux_profile(S .* cell_volumes, ΔV); 
    end

    return vars
end

for expname in ["iter129_bulkformula"]
    savename = datadir("native/" * expname * region * "_WaterProp_budget_z.jld2")
    θS = filter_waterprops_budget(diagpath, expname, γ)
    jldsave(savename, θS = θS)
end


# vm = maximum(abs.(θ["iter129_bulkformula"] - mean(θ["iter129_bulkformula"], dims = 2)))
fig, ax = plt.subplots(2, 1, figsize = (15, 10))
savename = datadir("native/iter129_bulkformula" * region * "_WaterProp_budget_z.jld2")
θS = load(savename)["θS"]
kmax = 30; 
θ = θS["θ"][1:kmax, :]; S = θS["S"][1:kmax, :]
θ = θ .- mean(θ, dims = 2)
S = S .- mean(S, dims = 2)

v = 0.5; levels =collect(-v:0.1:v)

ax[1].contourf(tecco, z[1:kmax], θ, cmap = cm.balance, vmin = -v, vmax = v, levels = levels, extend = "both"); 
CS = ax[1].contour(tecco, z[1:kmax], θ, vmin = -v, vmax = v, levels = levels, colors = "k")
ax[1].clabel(CS, inline=1, fontsize=15, inline_spacing = 10)
ax[1].set_title("Temperature Anomaly")

v = 0.1; levels =collect(-v:0.025:v)
ax[2].contourf(tecco, z[1:kmax], S, cmap = cm.delta, vmin = -v, vmax = v, levels = levels, extend = "both");
CS = ax[2].contour(tecco, z[1:kmax], S, vmin = -v, vmax = v, levels = levels, colors = "k")
ax[2].clabel(CS, inline=1, fontsize=15, inline_spacing = 10)
ax[2].set_title("Salinity Anomaly")
[a.set_xlabel("time") for a in ax]
[a.set_ylabel("depth") for a in ax]
fig.tight_layout()
fig