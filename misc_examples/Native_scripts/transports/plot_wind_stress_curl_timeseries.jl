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
@pyimport seaborn as sns;
sns.set_theme(context = "notebook", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
colors =  sns.color_palette("deep")[1:4]

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
tecco = 1992+1/24:1/12:2018

ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region, extent = "not", include_bering = false)
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));
area_mask = area .* PAC_msk

ω = 2π / 86400; r = 6.378e6
beta = 2 * ω * cosd(45) / r ; rho = 1000

fname = datadir("native/" * region * "_areamean_τcurl.jld2")
curlτ_dict = load(fname)["curlτ_dict"]
fig, axes = plt.subplots(1, figsize = (10, 4), sharey = false, sharex = true)
curlτ_dict["Difference"] = curlτ_dict["iter129_bulkformula"] .- curlτ_dict["seasonalclimatology"]

ϕ_min_mask, ϕ_max_mask = OHC_helper.get_ϕ_max_min_mask(region, Γ, λ, ϕ, basins, basin_list)

for (i, expname) in enumerate(["iter129_bulkformula", "seasonalclimatology", "Difference"])
    axes.plot(tecco, 100e3 .* curlτ_dict[expname] .* inv(beta * rho), linewidth = 1, alpha =1, label = expname, c = colors[i])
    axes.set_title("Area-Averaged North Pacific Wind Stress Curl")
end
axes.legend()
fig
sum(area .* PAC_msk)
cosd(90)