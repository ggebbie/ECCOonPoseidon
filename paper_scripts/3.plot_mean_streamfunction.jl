include("../../../src/intro.jl")

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


vars =  ["iter0_bulkformula" "only_wind";
         "mean_tau_noadjust_redo_bf" "mean_tau_yesadjust_redo_bf";
         "only_buoyancy" "only_init"; 
         "only_kappa" "iter129_bulkformula"][:]
Ψ_means = Dict()
for (i, expname) in enumerate(vars)
    read_file = jldopen(datadir("Ψ_EulBol_timeseries_"*region*"_" * expname *".jld2"))
    Ψ_exp = read_file["Ψ_exp_timeseries"]
    Ψ_mean = mean(Ψ_exp, dims = 3)[:, :, 1]
    Ψ_means[expname] = Ψ_mean

end


ϕ_avg = jldopen(datadir("Ψ_Eul_timeseries_"*region*"_" * "only_init" *".jld2"))["ϕ_avg"]
# plot_labels_list = ["ECCO V4r4 First-Guess\n(Iteration 0)", "Adjusted Wind", 
# "Seasonal First-Guess Wind", 
# "Seasonal Adjusted Wind"]
vars =  ["iter0_bulkformula" "only_wind";
         "mean_tau_noadjust_redo_bf" "mean_tau_yesadjust_redo_bf";
         "only_buoyancy" "only_init"; 
         "only_kappa" "iter129_bulkformula"][:]

# exps =  ["iter0_bulkformula" "iter129_bulkformula", "only_buoyancy", 
# "only_init", "only_kappa", "only_wind"]
plot_labels_list = ["ECCO V4r4 First-Guess\n(Iteration 0)" "Adjusted Wind";
        "Seasonal First-Guess Wind" "Seasonal Adjusted Wind"; 
        "Adjusted Buoyancy" "Adjusted I.C."; 
        "Adjusted Mixing" "ECCO V4r4 Official Release\n(Iteration 129)"][:]

fig,axes=plt.subplots(4, 2, figsize = (17.5, 25), sharey = true, sharex = true)

# , 
# "Adjusted\nBuoyancy", "Adjusted\nI.C.", "Adjusted\nMixing", 
# "Adjusted\nWind"]

cms = []

# Plot a thicker arrow
iϕ = Base.findmin(abs.(ϕ_avg .- 22.5))[2]
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= -z[:].<= uplvl)
iz1 = lvls[1]
iz2 = lvls[end]
Ψ_bounds = 10.0
levels = collect(-Ψ_bounds:1:Ψ_bounds)

for (i, expname) in enumerate(vars)
    # ax = axes[:][i]
    # Ψ_exp = Ψ_means[expname]
    print(i)
    ax = axes[:][i]
    ax.set_title(plot_labels_list[i])
    Ψ_exp =  1e-6 .* Ψ_means[expname]
    ax.set_facecolor("black")

    cf = ax.contourf(ϕ_avg, z,  Ψ_exp, cmap=cmos.delta,levels = levels, 
    vmin = -1.5*Ψ_bounds, vmax = 1.5*Ψ_bounds, extend = "both")
    push!(cms, cf)

    Wout = -Ψ_exp[iz1, iϕ]; Wout = -round(Wout, digits = 2)
    Win = -Ψ_exp[iz2, iϕ]; Win = round(Win, digits = 2)
    Vout = -(Ψ_exp[iz1, iϕ] .- Ψ_exp[iz2, iϕ]); Vout = round(Vout, digits = 2)
    println(Wout, " ", Win, " ", Vout)
    cs2 = ax.contour(ϕ_avg, z[42:46], Ψ_exp[42:46, :], colors="k",levels = [-9.0, -2.0 ])
    ax.contour(ϕ_avg, z[1:42], Ψ_exp[1:42, :], colors="k",levels = [-9.0, -2.0 ])
    ax.contour(ϕ_avg, z[46:end], Ψ_exp[46:end, :], colors="k",levels = [-9.0, -2.0 ])

    labels = ax.clabel(cs2, fontsize=23, inline=true, fmt = "%.1f", 
    inline_spacing = 12, rightside_up = true, use_clabeltext = true)
    for label in labels
        text = label.get_text()
        label.set_fontweight("bold")
        label.set_fontsize(20)  # Adjust the fontsize as desired
        label.set_zorder(100)    # Set the zorder value

    end
    rect = patches.Rectangle((23, 2000), 59, 1000, linewidth=3, edgecolor="red",facecolor="none", alpha = 0.7)
    ax.add_patch(rect)

    # ax.annotate("Thicker Arrow", xy=(45, 1700), xytext=(45, 2000), arrowprops=arrow_params)
    arrow_params = Dict("arrowstyle"=>"->",  "color"=>"red", "linewidth"=>4, "relpos"=>(0.5, 0.0) )  # Increase the linewidth to make the arrow thicker
    ax.annotate("", xy=(40, 1000), xytext=(40, 2000), horizontalalignment="center", arrowprops=arrow_params, transform=ax.transData)
    ax.text(42, 1800, string(abs(Wout)) * " Sv", transform=ax.transData, fontsize=20, color = "red")
    ax.annotate("", xy=(40, 2900), xytext=(40, 4000), horizontalalignment="center", arrowprops=arrow_params, transform=ax.transData)
    ax.text(42, 3550, string(abs(Win)) * " Sv", transform=ax.transData, fontsize=20, color = "red")
    arrow_params = Dict("arrowstyle"=>"->",  "color"=>"red", "linewidth"=>4, "relpos"=>(0, 0.4) )  # Increase the linewidth to make the arrow thicker
    ax.annotate("", xy=(5, 2500), xytext=(23, 2500), horizontalalignment="center", arrowprops=arrow_params, transform=ax.transData)
    ax.text(4, 2200, string(abs(Vout)) * " Sv", transform=ax.transData, fontsize=20, color = "red")


    ax.set_xticks(-40:20:60)
    ax.set_xlim(-34, 60)
    ax.set_ylim(800, 5900)

    ax.set_title(plot_labels_list[i])
    ax.tick_params(which="both", bottom=true)



end
fig
axes[1].invert_yaxis()

[a.set_ylabel("Depth [m]", fontweight = "bold") for a in axes[:, 1][:]]
[a.tick_params(which="both", bottom=true, left = true) for a in axes[:, 1][:]]

# axes[2].set_ylabel("Depth [m]", fontweight = "bold")

lab = string.(abs.(collect(-40:20:60)))
lab = lab .* ["°S", "°S", "", "°N", "°N", "°N"]
axes[2].set_xticklabels(lab)
axes[4].set_xticklabels(lab)
axes[end, 1].set_xlabel("Latitude", fontweight = "bold")
axes[end, 2].set_xlabel("Latitude", fontweight = "bold")

# axes[2].tick_params(which="both", bottom=true, left = true)

fig.subplots_adjust(wspace = 0.05, hspace = .25)

fig_labs = uppercase.(["a" "b"; "c" "d"; "e" "f"; "g" "h"][:])

for (i, a) in enumerate(axes)
    a.annotate(fig_labs[i], (0.93, 0.02), fontsize = 25, color = "white", 
    xycoords="axes fraction", fontweight = "bold")
end
fig


fig.savefig(plotsdir("wind_rewrite/3.streamfunctions.png"), bbox_inches = "tight", dpi = 400)

