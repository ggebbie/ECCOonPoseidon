include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings
using .OHC_helper
using ColorSchemes
import NaNMath as nm
import PyPlot as plt
using PyCall

include(srcdir("config_exp.jl"))
@pyimport cmocean.cm as cm
@pyimport seaborn as sns;

colors =  sns.color_palette("deep")[1:4]
labels = ["ECCO", "ECCO Forcing", "ECCO Init. Cond. and κ", "CTRL"]
sns.set_theme(context = "notebook", style = "white", font_scale = 1.0,
              palette = sns.color_palette("deep"));

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
 
ocean_mask = wet_pts(Γ)
region = "PAC"; 
tecco = collect(1992+1/24:1/12:2018)
E,F = trend_matrices(tecco[tstart:tstop])

vcat(1:36,312-36+1:312)
XΨ=collect(-89.0:89.0); Y=reverse(z); #coordinate variables

tecco = 1992+1/24:1/12:2018
nz = 50; tstart = 1; tstop = 312; nt = tstop - tstart + 1;

Ψ_exp = jldopen(datadir("Ψtimeseries_"*region*".jld2"))["Ψ_exp"]
Ψ_exp = Dict(key => reverse(permutedims(value, (2, 1, 3)), dims = 1) for (key, value) in Ψ_exp)
Ψ_exp = Dict(key => -mean(value[44:46, 80:90, :], dims = [1,2]) for (key, value) in Ψ_exp) #draws a circle around the maximum overturning
fig,axes=plt.subplots(1, 5, figsize = (10, 3), sharex = true, sharey = true)
for (i, expname) in enumerate(keys(shortnames))
    axes[i].plot(tecco, 1e-6 .* Ψ_exp[expname][1, 1, :], label = expname)
    axes[i].set_title(expname)
end
fig.suptitle("Pacific Streamfunction Maximum Evolution in ECCO [Sv per century]")
fig.tight_layout()
fig
fig.savefig(plotsdir("native/StreamfunctionMaximum_" *  region * ".png"), dpi = 500)