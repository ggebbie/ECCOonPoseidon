#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

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
# NPAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
# region, extent = "not")
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")

region = "NPAC"
# NoNPAC_msk = PAC_msk .- NPAC_msk

cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)

# GTF = OHC_helper.get_geothermalheating(γ, Γ) .* PAC_msk
nz = 50

ΔV = zeros(Float32, nz)
GTF_z=  zeros(Float32, nz)
[ΔV[k] = Float32(sum(cell_volumes[:, k])) for k=lvls]
# [GTF_z[k] = sum(GTF[:, k]) for k=lvls]
GTF_z = repeat(GTF_z, 1, 312)
V = sum(ΔV)

integrate(x) = cumsum(hcat([0.0], x .* 2.628e+6), dims=2)
volume_weight_column(x) =  sum(x .* ΔV, dims =1) ./ V 
sum_fluxes(d) = sum(d[varn][:] for varn in ["κxyθ", "uvθ", "wθ", "κzθ",  "GTF"])
anomaly(x) = x .- mean(x)

fig,axs=plt.subplots(1,4, sharey = true, figsize = (25, 5))
labels = ["Iteration 129", "No Init. Adjust", "No Forcing Adjust", "Iteration 0"]
alpha = [1, 0.5, 0.5, 1]
@time for (i, expname) in enumerate(keys(shortnames))

    # @time for (i, expname) in enumerate(["iter129_bulkformula"])
    fname = datadir(expname * region * "_THETA_budget_z.jld2")
    θ_budg = load(fname)["dθ"]

    θ_budg_int = Dict()
    for var in keys(θ_budg)
        tmp = volume_weight_column(θ_budg[var])
        θ_budg[var] = tmp
        if var in ["κxyθ", "uvθ", "wθ", "κzθ"] 
            θ_budg_int[var] = integrate(θ_budg[var])[1:end-1]
        end
    end

    # θ_budg_int["GTF"] = volume_weight_column(GTF_z)
    θ_budg_int["GTF"] = volume_weight_column(GTF_z)
    θ_budg_int["sum"] = sum_fluxes(θ_budg_int)

    θ = vec(θ_budg["θ"])
    θ_approx = vec(θ_budg_int["sum"])
    diffusion = vec(θ_budg_int["κxyθ"] .+ θ_budg_int["κzθ"]) .+ θ_budg_int["GTF"]
    advection = vec(θ_budg_int["uvθ"] .+ θ_budg_int["wθ"])

    axs[1].plot(tecco, 100 .* (θ .- θ[1]), color = colors[i], alpha = alpha[i]); 
    axs[2].plot(tecco, 100 .* θ_approx, color = colors[i], label = labels[i], alpha = alpha[i])
    axs[3].plot(tecco, 100 .* advection , color = colors[i], label = labels[i], alpha = alpha[i])
    axs[4].plot(tecco, 100 .* diffusion, color = colors[i], alpha = alpha[i])
end

fig
# end
# fig.suptitle("North Pacific Temperature Budget")
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

fig.savefig(plotsdir("AnomalyHeatBudgetSameAxis_" * region * ".png"),
 dpi = 900, bbox_inches = "tight")
# fig