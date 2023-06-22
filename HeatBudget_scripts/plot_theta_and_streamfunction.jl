include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings
using .OHC_helper
using Plots
using ColorSchemes
import NaNMath as nm
import PyPlot as plt
using PyCall

include(srcdir("config_exp.jl"))
@pyimport cmocean.cm as cm
@pyimport seaborn as sns;
@pyimport pandas as pd;

subtract_dict(D, subkey) = Dict(key => value .- D[subkey] for (key, value) in D)
mask_dict(D, lvls) = Dict(key => value[lvls, :] for (key, value) in D)

colors =  sns.color_palette("deep")[1:4]
sns.set_theme(context = "talk", style = "white", font_scale = 1.0,
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
#first 3 yrs tecco[1:36]
#last 3 yrs tecco[end-36+1:end]
vcat(1:36,312-36+1:312)
XΨ=collect(-89.0:89.0); Y=reverse(z); #coordinate variables
Xθ = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)

insuffix = "sfctobot"
outsuffix = "1tobot"
uplvl = -1e3; botlvl = -Inf
lvls = findall( botlvl .<= z[:].<= uplvl)
Y_mask = reverse(z[lvls])

Ψ_exp_mean = jldopen(datadir("Ψ"*region*"_" * insuffix * ".jld2"))["Ψ_exp_mean"]
Ψ_exp_mean = Dict(key => reverse(value', dims = 1) for (key, value) in Ψ_exp_mean)
Ψ_exp_mean = mask_dict(Ψ_exp_mean, lvls)
Ψ_exp_mean_anom0 = subtract_dict(Ψ_exp_mean, "iter0_bulkformula") 
Ψ_exp_labels = [L"\Psi^{\Delta F, \Delta T }", L"\Psi^{\Delta F}",
                L"\Psi^{\Delta T}", L"\Psi^{0}"]

θ_zonal = jldopen(datadir("θzonal"*region*"_" * insuffix * ".jld2"))["θ_zonal"]
θ_zonal_mean = Dict(key => mean(value, dims = 3)[:, :, 1] for (key, value) in θ_zonal)
θ_zonal_mean = mask_dict(θ_zonal_mean, lvls)
θ_zonal_mean_anom0 = subtract_dict(θ_zonal_mean, "iter0_bulkformula")
θ_labels = [L"\theta^{\Delta F, \Delta T}", L"\theta^{\Delta F}",L"\theta^{\Delta T}", L"\theta^{0}"]


θ_bounds = maximum(abs.(extrema(θ_zonal_mean_anom0)))
Ψ_bounds = 1e-6.* maximum(abs.(extrema(Ψ_exp_mean_anom0)))
Ψlevels = LinRange(-Ψ_bounds, Ψ_bounds, 13)
θlevels = LinRange(-θ_bounds, θ_bounds, 25)

fig,ax=plt.subplots(4,1, figsize = (12.5, 20), sharey = "row")
cp = Vector{Any}(missing, 1)
θ_labels_anom = [sym * L" - \theta^{0}" for sym in θ_labels]
Ψ̄_anom_labels = [exp * L"- \Psi^{0}" for exp in Ψ_exp_labels]      

@time for (i, expname) in enumerate(keys(shortnames))
    θ_zonal_exp = θ_zonal_mean_anom0[expname]
    Ψ_exp = Ψ_exp_mean_anom0[expname]
    cp[1] = ax[i].contourf(Xθ, Y_mask, reverse(θ_zonal_exp, dims = 1), 
                    levels = θlevels, cmap = cm.balance)
    cs = ax[i].contour(XΨ, Y_mask, reverse(1e-6.* Ψ_exp, dims = 1), 
                    levels = Ψlevels, colors="k")
    ax[i].clabel(cs, fontsize=15, inline=true, fmt = "%.2f")
    ax[i].set_xlim(-39, 59)
    ax[i].set_xticks(-40:10:60)
    ax[i].set_title(Ψ̄_anom_labels[i] * " (Solid); " * θ_labels_anom[i] *" (Colored)")

end
[a.set_ylabel("Depth [m]") for a in ax]
[a.set_xlabel("Latitude") for a in ax]

fig.colorbar(cp[1], ax = ax[4], orientation = "horizontal")
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/MassTransport/" * "TrspNθAnomiter0_" * region * "_" * outsuffix * ".png",
dpi = 1000)


fig,axs=plt.subplots(4,2, figsize = (7.5, 15))
θ_bounds = maximum(abs.(extrema(θ_zonal_mean)))
Ψ_bounds = 1e-6.* maximum(abs.(extrema(Ψ_exp_mean)))
Ψlevels = LinRange(-Ψ_bounds, Ψ_bounds, 13)
θlevels = LinRange(-θ_bounds, θ_bounds, 25)
cp = Vector{Any}(missing, 1)
cs = Vector{Any}(missing, 1)

@time for (i, expname) in enumerate(keys(shortnames))
    θ_zonal_exp = θ_zonal_mean[expname]
    Ψ_exp = Ψ_exp_mean[expname]
    ax = axs[i, :]
    cp[i] = ax[i].contourf(Xθ, Y_mask, reverse(θ_zonal_exp, dims = 1), 
                    levels = θlevels, cmap = cm.balance)
    cs[i] = ax[i].contourf(XΨ, Y_mask, reverse(1e-6.* Ψ_exp, dims = 1), 
                    levels = Ψlevels, cmap = cm.delta)
    # ax[i].clabel(cs, fontsize=15, inline=true, fmt = "%.2f")
    # fig.colorbar(cp) # Add a colorbar to a plot
end
[ax.set_xlim(-39, 59) for ax in axs]
[ax.set_xticks(-40:10:60) for ax in axs]
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/MassTransport/" * "TrspNextθ_" * region * "_" * outsuffix * ".png",
dpi = 1000)

using JLD2
# close("all")

