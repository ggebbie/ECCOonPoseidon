include("../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour, PyCall, 
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))
@pyimport seaborn as sns;
@pyimport matplotlib.patches as patches

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

include(srcdir("plot_and_dir_config.jl"))

@pyimport cmocean.cm as cmos


region = "PAC"; 
tecco = collect(1992+1/24:1/12:2018)

tecco = 1992+1/24:1/12:2018
nz = 50; tstart = 1; tstop = 312; nt = tstop - tstart + 1;

sns.set_theme(context = "poster", style = "white", font_scale = 1.0,
              palette = sns.color_palette("deep"));

uplvl = -1e3; botlvl = -Inf
lvls = findall( botlvl .<= z[:].<= uplvl)
Y_mask = reverse(z[lvls])

cp = Vector{Any}(missing, 1)  
# alabels = ["Iteration 129", "Iteration 0", "Initial 129", "Initial 0"]
clevels = ["-9.0", "-7.5", "-1.5"]


vars =  ["only_init", "only_kappa", "only_sfc", "iter129_bulkformula",  "iter0_bulkformula"]
Ψ_means = Dict()
for (i, expname) in enumerate(vars)
    read_file = jldopen(datadir("Ψ_EulBol_timeseries_"*region*"_" * expname *".jld2"))
    Ψ_exp = read_file["Ψ_exp_timeseries"]
    Ψ_mean = mean(Ψ_exp, dims = 3)[:, :, 1]
    Ψ_means[expname] = Ψ_mean

end

Ψ_means["Δ"] = Ψ_means["iter129_bulkformula"] .- Ψ_means["iter0_bulkformula"]

ϕ_avg = jldopen(datadir("Ψ_Eul_timeseries_"*region*"_" * "only_init" *".jld2"))["ϕ_avg"]

vars =  ["iter0_bulkformula","iter129_bulkformula", "Δ"]
plot_labels["Δ"] = "Δ"
for (i, expname) in enumerate(vars)
    Ψ_exp = 1e-6.* Ψ_means[expname]
    fig,ax=plt.subplots(1, 1, figsize = (12, 7))
    ax.set_facecolor("black")
    Ψ_bounds = 10.5
    levels = collect(-Ψ_bounds:1.5:Ψ_bounds)
    if expname == "Δ"
        Ψ_bounds = 3
        levels = -Ψ_bounds:0.5:Ψ_bounds
    end
    ax.contourf(ϕ_avg, z,  Ψ_exp, cmap=cmos.delta,levels = levels, 
    vmin = -1.5*Ψ_bounds, vmax = 1.5*Ψ_bounds, extend = "both")
    cs2 = ax.contour(ϕ_avg, z[30:46], Ψ_exp[30:46, :], colors="k",levels = levels)
    ax.contour(ϕ_avg, z[1:30], Ψ_exp[1:30, :], colors="k",levels = levels)
    ax.contour(ϕ_avg, z[46:end], Ψ_exp[46:end, :], colors="k",levels = levels)
    if expname != "Δ"
        labels = ax.clabel(cs2, levels = levels[2:8], fontsize=20, inline=true, fmt = "%.1f", 
        inline_spacing = 12, rightside_up = true, use_clabeltext = true)
    else
        labels = ax.clabel(cs2, fontsize=20, inline=true, fmt = "%.1f", 
        inline_spacing = 12, rightside_up = true, use_clabeltext = true)

    end



    for label in labels
        text = label.get_text()
        label.set_fontsize(23)  # Adjust the fontsize as desired
        label.set_fontweight("medium")
        if (text in clevels) && (expname != "Δ")
            label.set_fontweight("bold")
            label.set_fontsize(22)  # Adjust the fontsize as desired
            label.set_zorder(100)    # Set the zorder value
        end
    end
    rect = patches.Rectangle((23, 2000), 59, 1000, linewidth=3, edgecolor="red",facecolor="none", alpha = 0.7)
    ax.add_patch(rect)
    rect = patches.Rectangle((2.5, 4000), 0.1, 1500, linewidth=3, edgecolor="yellow",facecolor="none", alpha = 0.7)
    ax.add_patch(rect)

    ax.invert_yaxis()
    ax.set_xticks(-40:20:60)
    ax.set_xlim(-34, 60)
    lab = string.(abs.(collect(-40:20:60)))
    lab = lab .* ["°S", "°S", "", "°N", "°N", "°N"]
    ax.set_xticklabels(lab)
    ax.set_title(plot_labels[expname])
    ax.set_ylabel("Depth [m]", fontweight = "bold")
    ax.set_xlabel("Latitude", fontweight = "bold")

    fig.savefig(plotsdir("native/paper_figures/ΨEulBol" * expname * ".png"), bbox_inches = "tight", dpi = 400)

end
fig



i = 1
expname = "iter129_bulkformula"
Ψ_exp = 1e-6.* Ψ_means[expname]
fig,ax=plt.subplots(1, 1, figsize = (12, 7))
ax.set_facecolor("black")
Ψ_bounds = 10.5
levels = collect(-Ψ_bounds:1.0:Ψ_bounds)
if expname == "Δ"
    Ψ_bounds = 3
    levels = -Ψ_bounds:0.5:Ψ_bounds
end
ax.contourf(ϕ_avg, z,  Ψ_exp, cmap=cmos.delta,levels = levels, 
vmin = -1.5*Ψ_bounds, vmax = 1.5*Ψ_bounds, extend = "both")
cs2 = ax.contour(ϕ_avg, z[30:46], Ψ_exp[30:46, :], colors="k",levels = levels)
ax.contour(ϕ_avg, z[1:30], Ψ_exp[1:30, :], colors="k",levels = levels)
ax.contour(ϕ_avg, z[46:end], Ψ_exp[46:end, :], colors="k",levels = levels)

labels = ax.clabel(cs2, fontsize=30, inline=true, fmt = "%.1f", 
inline_spacing = 30, rightside_up = true, use_clabeltext = true)

fig

for label in labels
    text = label.get_text()
    label.set_fontsize(23)  # Adjust the fontsize as desired
    label.set_fontweight("medium")
    if (text in clevels) && (expname != "Δ")
        label.set_fontweight("bold")
        label.set_fontsize(22)  # Adjust the fontsize as desired
        label.set_zorder(100)    # Set the zorder value
    end
end
fig
rect = patches.Rectangle((23, 2000), 59, 1000, linewidth=3, edgecolor="red",facecolor="none", alpha = 0.7)
ax.add_patch(rect)
rect = patches.Rectangle((2.5, 4000), 0.1, 1500, linewidth=3, edgecolor="yellow",facecolor="none", alpha = 0.7)
ax.add_patch(rect)

ax.invert_yaxis()
ax.set_xticks(-40:20:60)
ax.set_xlim(-34, 60)
lab = string.(abs.(collect(-40:20:60)))
lab = lab .* ["°S", "°S", "", "°N", "°N", "°N"]
ax.set_xticklabels(lab)
ax.set_title(plot_labels[expname])
ax.set_ylabel("Depth [m]", fontweight = "bold")
ax.set_xlabel("Latitude", fontweight = "bold")
fig