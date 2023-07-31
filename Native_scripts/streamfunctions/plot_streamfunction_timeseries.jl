include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour, PyCall, 
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings
using .OHC_helper
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))
@pyimport seaborn as sns;

colors =  sns.color_palette("deep")[1:4][[1, 3, 4]]
labels = ["ECCO", "ECCO Forcing", "ECCO Init. Cond. and κ", "CTRL"]
sns.set_theme(context = "poster", style = "white", font_scale = 1.0,
              palette = sns.color_palette("deep"));

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
# _, vastdiagpath = listexperiments(vastexprootdir())
# update_paths_dict1!(diagpath, vastdiagpath)
# update_paths_dict!(diagpath, vastdiagpath)

runpath,diagpath = listexperiments(vastexprootdir())[2]
ocean_mask = OHC_helper.wet_pts(Γ)
region = "PAC"; 
tecco = collect(1992+1/24:1/12:2018)
E,F = trend_matrices(tecco)

vcat(1:36,312-36+1:312)
XΨ=collect(-89.0:89.0); Y=reverse(z); #coordinate variables

tecco = 1992+1/24:1/12:2018
nz = 50; tstart = 1; tstop = 312; nt = tstop - tstart + 1;

myexps = ["iter0_bulkformula", "iter129_bulkformula", "noinitadjust", 
"seasonalclimatology", "nosfcadjust",
"seasonalclimatology_iter0", "climatological_tau"]


myexps = ["iter129_bulkformula","noinitadjust", "iter0_bulkformula"]
label = ["i129", "i129 Forcing", "i0"]

fig,axes=plt.subplots(2, 1, figsize = (15, 7.5), sharex = true, sharey = false)
Ψ_maximum = Dict()
Ψ_NPAC_W = Dict()

for (i, expname) in enumerate(myexps)
    Ψ_exp = jldopen(datadir("ΨwBolustimeseries_"*region*"_" * expname *".jld2"))["Ψ_exp"]
    Ψ_exp = reverse(permutedims(Ψ_exp, (2, 1, 3)), dims = 1)
    Ψ_maximum[expname] = -mean(Ψ_exp[44:46, 87:100, :], dims = [1,2])[1, 1, :]
    Ψ_NPAC_W[expname] = -Ψ_exp[42, 114, :]

    # Ψmean = round(Float32(mean(Ψ_maximum[expname])), sigdigits = 2)
    # println(Ψmean)
    axes[1].plot(tecco, 1e-6 .* Ψ_maximum[expname], 
    label = label[i], 
    color = colors[i], lw = 3)
    println(mean(Ψ_NPAC_W[expname]))
    axes[2].plot(tecco, 1e-6 .* Ψ_NPAC_W[expname], 
    label = label[i], 
    color = colors[i], lw = 3)

end
axes[1].legend(frameon = false, ncols = 3)
axes[1].set_title("Pacific Streamfunction Maximum")
axes[2].set_title("NPAC Upwelling")

fig.tight_layout()
fig
# fig.savefig(plotsdir("native/StreamfunctionMaximum_" *  region * ".png"), dpi = 500)