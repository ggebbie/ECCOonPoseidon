#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, FLoops
using .OHC_helper
import PyPlot as plt 
using PyCall
include(srcdir("config_exp.jl"))
cmo = pyimport("cmocean.cm");

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
diagpath["seasonalclimatology_iter0_multarr0"] = "/batou/ECCOv4r4/exps/seasonalclimatology_iter0/run_multarr0/diags/"

# abbreviations for each experiment for labels, etc.
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
nexp = length(shortnames) # number of experiments
 
ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));

H = OHC_helper.smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H
tecco = 1992+1/24:1/12:2018
nz = 50

ΔV = zeros(Float32, nz)
[ΔV[k] = Float32(sum(cell_volumes[:, k])) for k=1:nz]
lvls = findall( -3000 .<= z[:].<= -2000)

masked_volume = cell_volumes[:, lvls]
sum_masked_volume = Float32(sum(masked_volume))
@pyimport seaborn as sns;

sns.set_theme(context = "notebook", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
              
fname = datadir("native/" * region * "_THETA_levels_" * ".jld2")
θ = load(fname)["θ"]
removetimemean(x) = x .- mean(x, dims = 2)
vm = 100 .* maximum(abs.(removetimemean(θ["iter129_bulkformula"][40:45, :])))
vm = round(vm)
levels = -vm:0.2:vm 
vm = 1.5*vm
titles = ["ECCO V4r4 with Full Forcing", "ECCO V4r4 with Climatological Forcing"]

fig, ax = plt.subplots(1, 2, figsize = (15, 5))
for (i, expname) in enumerate(["iter129_bulkformula", "iter0_bulkformula"])
    ax[i].contourf(tecco, abs.(z), 100 .* removetimemean(θ[expname]), levels = levels, 
    vmin = -vm, vmax = vm, cmap = cmo.balance, extend = "both")
    ax[i].contour(tecco, abs.(z[1:37]), 100 .* removetimemean(θ[expname][1:37, :]), levels = levels, 
    vmin = -vm, vmax = vm, colors = "k")
    CS = ax[i].contour(tecco, abs.(z[37:44]), 100 .* removetimemean(θ[expname][37:44, :]), levels = levels, 
    vmin = -vm, vmax = vm, colors = "k")
    ax[i].contour(tecco, abs.(z[44:end]), 100 .* removetimemean(θ[expname][44:end, :]), levels = levels, 
    vmin = -vm, vmax = vm, colors = "k")
    # ax[i].clabel(CS, inline=1, fontsize=12, inline_spacing = 10)
    ax[i].set_title(titles[i])
    ax[i].set_ylim(1000, 5000)
    ax[i].invert_yaxis()
    ax[i].set_xlabel("time"); ax[i].set_ylabel("Depth [m]")
end
fig.suptitle("North Pacific Temperature Anomaly [cK]")
fig
fig.savefig(plotsdir("native/TempAnomalies" * region * ".png"), dpi = 500)