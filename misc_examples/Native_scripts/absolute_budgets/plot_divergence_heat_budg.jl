#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings,
    PyCall, ColorSchemes
using .OHC_helper
import PyPlot as plt
@pyimport seaborn as sns;
@pyimport pandas as pd;

colors =  sns.color_palette("colorblind")[1:4]
labels = ["ECCO", "ECCO Forcing", "ECCO Init. and κ", "CTRL"]

sns.set_theme(context = "poster", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));

include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());

ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
nexp = length(shortnames) # number of experiments

tecco = 1992+1/24:1/12:2018

ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
NPAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
region = "PAC56"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")

region = "PAC56noNPAC"
NoNPAC_msk = PAC_msk .- NPAC_msk

cell_depths = OHC_helper.get_cell_depths(NoNPAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));

uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)

# GTF = OHC_helper.get_geothermalheating(γ, Γ) .* PAC_msk
nz = 50

ΔV = zeros(Float32, nz)
GTF_z=  zeros(Float32, nz)
[ΔV[k] = Float32(sum(cell_volumes[:, k])) for k=lvls]
# [GTF_z[k] = sum(GTF[:, k]) for k=lvls]
# GTF_z = repeat(GTF_z, 1, 312)
V = sum(ΔV)

integrate(x) = cumsum(hcat([0.0], x .* 2.628e+6), dims=2)
volume_weight_column(x) =  sum(x .* ΔV, dims =1) ./ V 
sum_fluxes(d) = sum(d[varn][:] for varn in ["κxyθ", "uvθ", "wθ", "κzθ",  "GTF"])
anomaly(x) = x .- mean(x)

expname = "iter129_bulkformula"; i = 1
# @time for (i, expname) in enumerate(["iter129_bulkformula"])
fname = datadir(expname * region * "_THETA_budget_z.jld2")
θ_budg = load(fname)["dθ"]

θ_budg_deriv = Dict()
for var in keys(θ_budg)
    tmp = volume_weight_column(θ_budg[var])
    θ_budg_deriv[var] = tmp
end

# θ_budg_deriv["GTF"] = volume_weight_column(GTF_z)[:]
θ_budg_deriv["GTF"] = zeros(312)
θ_budg_deriv["sum"] = sum_fluxes(θ_budg_deriv)

dθ = (vec(θ_budg_deriv["θ"])[2:end] .- vec(θ_budg_deriv["θ"])[1:end-1]) ./ 2.628e+6
dθ_approx = vec(θ_budg_deriv["sum"])

diffusion = vec(θ_budg_deriv["κxyθ"] .+ θ_budg_deriv["κzθ"]) .+ θ_budg_deriv["GTF"]
advection = vec(θ_budg_deriv["uvθ"] .+ θ_budg_deriv["wθ"])
tecco_int = (tecco[2:end] .+ tecco[1:end-1])./2
fig,axs=plt.subplots(1,4, sharey = true, figsize = (25, 5))
axs[1].plot(tecco_int,  dθ, color = colors[i]); 
axs[2].plot(tecco, dθ_approx, color = colors[i], label = labels[i])
axs[3].plot(tecco, advection , color = colors[i], label = labels[i])
axs[4].plot(tecco, diffusion, color = colors[i]); 
    
# end
# fig.suptitle("North Pacific Temperature Anomaly Budget")
axs[1].set_title("Model Output \n "*L"\theta", y=1.07)
axs[2].set_title("Reconstruction \n "*L"\theta", y=1.07)
axs[3].set_title("Advection \n " * L"\int F_{adv} dt", y=1.05)
axs[4].set_title("Diffusion + Geothermal \n " *L"\int F_{diff} + F_{geo} dt", y=1.05)
axs[1].set_ylabel("[cK]")
# axs[1].set_yticks([-1, -0.5, 0.0, 0.5, 1])
# axs[1].set_ylim(-1.2, 1.2)


axs[2].legend()
sns.move_legend(axs[2], "lower center", 
                bbox_to_anchor=(1.1, -0.34), ncol=4, 
                frameon=true, borderaxespad=0.)
fig.subplots_adjust(wspace = 0.3)

fig

# fig.savefig(plotsdir("AnomalyHeatBudgetDeriv_" * region * ".png"),
#  dpi = 600, bbox_inches = "tight")
# fig