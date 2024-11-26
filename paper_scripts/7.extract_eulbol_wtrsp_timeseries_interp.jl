using Pkg
Pkg.activate(".")
include("../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings
using Statistics
import NaNMath as nm
using PyCall
import PyPlot as plt
using NaNStatistics
@pyimport cmocean.cm as cmos
@pyimport matplotlib.ticker as ticker
@pyimport matplotlib.colors as colors_plt
@pyimport matplotlib.tri as tri
@pyimport matplotlib.gridspec as gridspec

include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

include(srcdir("plot_and_dir_config.jl"))
 
ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "false", include_bering = true)

region = "NPAC"; 
NPAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "false", include_bering = true)

tecco = 1992+1/24:1/12:2018

vars =  ["only_init", "only_kappa", "only_sfc", "iter129_bulkformula",  "iter0_bulkformula"]
vars =  ["only_buoyancy", "only_wind"]

uplvl = -2.4e3; botlvl = -2.6e3; suffix = "2to3"
lvls = findall( botlvl .<= -z[:].<= uplvl)

diagpath["mean_tau_noadjust_redo_bf"] = vastdiagdir("seasonalclimatology", "run_only_clim_iter0_tau_bf")
diagpath["mean_tau_yesadjust_redo_bf"] = vastdiagdir("seasonalclimatology", "run_only_clim_iter129_tau_bf2")

# (lon, lat, f,i,j,w) = load(datadir("0.5deg_interpolation_factors.jld2"), "lon", "lat", "f", "i", "j", "w")

lat_mask0 = (ϕ .>= -40)
lat_mask1 =  (ϕ .<= 65)
lat_mask = lat_mask0 .* lat_mask1
mask_dict = Dict()
mask_dict["Atlantic"] = basin_mask(["Atlantic"], γ); mask_dict["Atlantic"] .= mask_dict["Atlantic"] .* lat_mask
mask_dict["Pacific"] = basin_mask(["Pacific"], γ); mask_dict["Pacific"] .= mask_dict["Pacific"] .* lat_mask
mask_dict["Indian"] = basin_mask(["Indian"], γ); mask_dict["Atlantic"] .= mask_dict["Atlantic"] .* lat_mask
mask_dict["North Pacific"] = NPAC_msk

diff_res_dict = Dict()
diff_eul_dict = Dict()
diff_bol_dict = Dict()

diff_res0_dict = Dict()
diff_eul0_dict = Dict()
diff_bol0_dict = Dict()
for basin in keys(mask_dict)

    Wres_reg = Dict()
    Weul_reg = Dict()
    Wbol_reg = Dict()

    exps = ["mean_tau_noadjust_redo_bf", "mean_tau_yesadjust_redo_bf"]

    Wres_reg[exps[1]] = zeros(312)
    Wres_reg[exps[2]] = zeros(312)

    Weul_reg[exps[1]] = zeros(312)
    Weul_reg[exps[2]] = zeros(312)

    Wbol_reg[exps[1]] = zeros(312)
    Wbol_reg[exps[2]] = zeros(312)

    area_mask = mask_dict[basin] .* area

    for expname in exps
        fname = datadir("native/" * expname * "_W_residual_2500m.jld2")
        Wres = load(fname)["Wres"]

        fname = datadir("native/" * expname * "_W_eul_2500m.jld2")
        Weul = load(fname)["Weul"]

        fname = datadir("native/" * expname * "_W_bol_2500m.jld2")
        Wbol = load(fname)["Wbol"]

        print(size(Wres))
        for tt in 1:312
            Wres_reg[expname][tt] = sum(Wres[:, tt] .* area_mask)
            Weul_reg[expname][tt] = sum(Weul[:, tt] .* area_mask)
            Wbol_reg[expname][tt] = sum(Wbol[:, tt] .* area_mask)

        end
    end
    diff_res_dict[basin] = Wres_reg[exps[2]] .- Wres_reg[exps[1]]
    diff_eul_dict[basin] = Weul_reg[exps[2]] .- Weul_reg[exps[1]]
    diff_bol_dict[basin] = Wbol_reg[exps[2]] .- Wbol_reg[exps[1]]

    diff_res0_dict[basin] = Wres_reg[exps[1]]
    diff_eul0_dict[basin] = Weul_reg[exps[1]]
    diff_bol0_dict[basin] = Wbol_reg[exps[1]]

end

function average_every_k(lst::Vector{T}, k::Int) where T
    # Determine the number of full chunks
    n_chunks = div(length(lst), k)
    
    # Initialize an empty array to hold the averages
    averages = Vector{Float64}(undef, n_chunks)
    
    # Iterate over each chunk and calculate the average
    for i in 1:n_chunks
        chunk = lst[(i-1)*k + 1:i*k]
        averages[i] = mean(chunk)
    end
    
    # If there are leftover items that don't form a full chunk, 
    # you can optionally handle them (e.g., by averaging or ignoring)
    
    return averages
end

window = 12
tecco_avg = average_every_k(collect(tecco), window)
# Wres_reg["Difference"] = (Wres_reg["mean_tau_yesadjust_redo_bf"] .- Wres_reg["mean_tau_noadjust_redo_bf"])
fig, ax = plt.subplots(1, 4, figsize = (16, 5), sharey = true)
for (i, basin) in enumerate(keys(mask_dict))
    # ax[i].plot(tecco_avg, 1e-6 .* average_every_k(diff_eul_dict[basin], 12), label = "Eulerian Transport", alpha = 0.3, c = "r")
    # ax[i].plot(tecco_avg, 1e-6 .* average_every_k(diff_bol_dict[basin], 12), label = "Bolus Transport", alpha = 0.3, c = "b", linestyle = "-")
    ax[i].plot(tecco_avg, 1e-6 .* average_every_k(diff_res_dict[basin], window), label = "Residual Transport", c = "k", alpha = 0.7)
    ax[i].plot(tecco_avg, 1e-6 .* average_every_k(diff_res0_dict[basin], window), label = "Residual Transport")
    ax[i].plot(tecco_avg, 1e-6 .* average_every_k(diff_res0_dict[basin] .+ diff_res_dict[basin], window), label = "Residual Transport")

    # ax[i].legend()
    ax[i].set_title(basin * " Ocean")
    ax[i].grid()
end
fig

# Create a figure
fig = plt.figure(figsize=(14, 10))
# Create subplots with a GridSpec layout
gs = fig.add_gridspec(5, 6)

# Add subplots in the specified positions
ax1 = fig.add_subplot(py"$gs[1:3, 0:2]")  # Top left
ax2 = fig.add_subplot(py"$gs[1:3, 2:4]")  # Top right
ax3 = fig.add_subplot(py"$gs[1:3, 4:]")  # Bottom left
ax4 = fig.add_subplot(py"$gs[3:5, 1:5]")  # Bottom left
ax = [ax1, ax2, ax3, ax4]

window = 12
tecco_avg = average_every_k(collect(tecco), window)
# Wres_reg["Difference"] = (Wres_reg["mean_tau_yesadjust_redo_bf"] .- Wres_reg["mean_tau_noadjust_redo_bf"])
# fig, ax = plt.subplots(1, 4, figsize = (16, 5), sharey = true)

plot_labels_list = [
"Climatological First-Guess Wind Stress", 
"Climatological Adjusted Wind Stress", "Response to Climatological\nWind Stress Adjustments"]

exp_colors["mean_tau_noadjust_redo_bf"] = exp_colors["iter0_bulkformula"] 
exp_colors["mean_tau_yesadjust_redo_bf"] = exp_colors["only_wind"] 
exp_colors["Difference"] = "k"


for (i, basin) in enumerate(["Pacific", "Indian", "Atlantic", "North Pacific"])
    # ax[i].plot(tecco_avg, 1e-6 .* average_every_k(diff_eul_dict[basin], 12), label = "Eulerian Transport", alpha = 0.3, c = "r")
    # ax[i].plot(tecco_avg, 1e-6 .* average_every_k(diff_bol_dict[basin], 12), label = "Bolus Transport", alpha = 0.3, c = "b", linestyle = "-")
    ax[i].plot(tecco_avg, 1e-6 .* average_every_k(diff_res0_dict[basin], window), label = plot_labels_list[1], alpha = 0.5, c = exp_colors["mean_tau_noadjust_redo_bf"])
    ax[i].plot(tecco_avg, 1e-6 .* average_every_k(diff_res0_dict[basin] .+ diff_res_dict[basin], window), label = plot_labels_list[2], alpha = 0.5, c = exp_colors["mean_tau_yesadjust_redo_bf"])
    ax[i].plot(tecco_avg, 1e-6 .* average_every_k(diff_res_dict[basin], window), label = plot_labels_list[3], c = "k", alpha = 0.7)

    # ax[i].legend()
    ax[i].set_title(basin * " Ocean")
    ax[i].set_ylabel("Upward Transport [Sv]", fontweight = "bold")
    ax[i].set_xlabel("time", fontweight = "bold")

    ax[i].grid()
end
fig.tight_layout()
# fig.subplots_adjust(wspace = 0.07)
ax[4].legend(loc = "lower center", fontsize=15,
              frameon = false, bbox_to_anchor = (0.5, -0.6), 
              ncols = 2)
fig.suptitle("Upward Vertical Transport through 2500 meters", y = 0.87, fontsize=25)
fig.savefig(plotsdir("wind_rewrite/7.vertical_transports_examples.png"), 
            bbox_inches = "tight", dpi = 400)
fig


fig