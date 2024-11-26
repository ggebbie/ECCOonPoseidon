include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall, 
    RollingFunctions, LinearAlgebra, Statistics, FFTW, LaTeXStrings
import PyPlot as plt 
import NaNMath as nm
include(srcdir("config_exp.jl"))
@pyimport matplotlib.patches as patches

include(srcdir("MeshArraysPlots.jl"))

tecco = 1992+1/24:1/12:2018
nz = 50

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
λ_wrap = wrap_λ(λ)
ocean_mask = wet_pts(Γ)
@pyimport cmocean.cm as cmo

region = "PAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

region = "NPAC"
NPAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
pad_arr(x) = vcat(reverse(x), x, reverse(x))
unpad_arr(x, arr_len) = x[arr_len+1:2*arr_len]
function low_pass(signal)
    nt = length(signal)
    signal_mean = mean(signal)
    ff = digitalfilter(Lowpass(1/24, fs = 1),Butterworth(2))
    filtered_signal = filtfilt(ff, pad_arr(signal .- signal_mean))
    return unpad_arr(filtered_signal, nt) .+ signal_mean
end


include(srcdir("plot_and_dir_config.jl"))

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); wherenan = (!isnan).(ϕ_avg .* 1)

normalize(x, dims) = (x .- mean(x, dims = dims)) ./ std(x, dims = dims)
normalize(x) = (x[:] .- mean(x[:])) ./ std(x[:])
square_error(x, y) = sum((x .- mean(x) .* (y .- mean(y))))
corr(x, y) = cov(x, y) / (std(x) * std(y))
# corr(x, y, dims) = cov(x, , normalize(y, dims); dims = dims) / (var(normalize(x); dims = dims) * var(normalize(y); dims = dims))

vars =  ["only_kappa",  "iter0_bulkformula"]
ϕ_avg = jldopen(datadir("Ψ_Eul_timeseries_PAC_only_init" *".jld2"))["ϕ_avg"]
Ws = Dict(); ΨEul = Dict(); ΨBol = Dict(); ΨEulBol = Dict()
z_ref = findall( -3300 .<= -z[:].<= -2000)
wθ_resid = Dict(); wθ_eul = Dict(); wθ_bol = Dict()

for (i, expname) in enumerate(vars)
    read_file = jldopen(datadir("Ψ_Eul_timeseries_PAC_" * expname *".jld2"))
    Ψ_exp = read_file["Ψ_exp_timeseries"]
    ΨEul[expname] = 1 .* Ψ_exp

    read_file = jldopen(datadir("Ψ_EulBol_timeseries_PAC_" * expname *".jld2"))
    Ψ_exp = read_file["Ψ_exp_timeseries"]
    ΨBol[expname] = Ψ_exp .- ΨEul[expname]
    ΨEulBol[expname] = 1 .* Ψ_exp
end

ΨEulBol["only_kappa"] .-= ΨEulBol["iter0_bulkformula"]
ΨEul["only_kappa"] .-= ΨEul["iter0_bulkformula"]
ΨBol["only_kappa"] .-= ΨBol["iter0_bulkformula"]

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 

Wres = -(ΨEulBol["only_kappa"][:, 1:end-2, :] .- ΨEulBol["only_kappa"][:, 3:end, :])
WEul = -(ΨEul["only_kappa"][:, 1:end-2, :] .- ΨEul["only_kappa"][:, 3:end, :])
WBol = -(ΨBol["only_kappa"][:, 1:end-2, :] .- ΨBol["only_kappa"][:, 3:end, :])

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 
area_avg = zonal_sum(PAC_msk .*area)
area_avg = area_avg[(!isnan).(ϕ_avg), :][:]
ϕ_avg = ϕ_avg[(!isnan).(ϕ_avg), :][:]
labels = ["ΔWres", "ΔWEul", "ΔWBol"]


Compensation = WEul .+ WBol ./ (abs.(WEul).+ abs.(WBol))
Compensation = mean(Compensation, dims = 3)[:, :, 1]

fig,axs=plt.subplots(figsize = (8, 6), sharex = true)
axs.set_facecolor("black")
levels = -1:0.1:1
CM = axs.contourf(ϕ_avg[2:end-1], z, Compensation,levels = levels, cmap=cmo.balance, 
vmin = -1, vmax = 1, extend = "both")
axs.invert_yaxis()
fig.colorbar(CM, orientation = "horizontal", fraction  =0.04)
axs.set_xticks(-40:20:60)
axs.set_xlim(-34, 60)
axs.set_ylabel("Depth [m]", fontweight = "bold")
lab = string.(abs.(collect(-40:20:60)))
lab = lab .* ["°S", "°S", "", "°N", "°N", "°N"]
axs.set_xticklabels(lab)
# axs.set_title(labels[i])
rect = patches.Rectangle((23, 2000), 59, 1000, linewidth=3, edgecolor="black",facecolor="none", alpha = 0.7)
axs.add_patch(rect)
fig
fig.savefig(plotsdir("native/paper_figures/ΔW_timemean_corr_EulnBol_mix.png"), bbox_inches = "tight", dpi = 400)




ΨCorr = zeros(size(Wres)[1:2])
for i in 1:size(Wres, 1), j in 1:size(Wres, 2)
    x = low_pass(WEul[i, j, :][:])
    y = low_pass(WBol[i, j, :][:])
    ΨCorr[i, j] = corr(x, y) 
end
    
fig,axs=plt.subplots(figsize = (8, 6), sharex = true)
ax2.set_facecolor("black")
levels = -1:0.2:1
CM = axs.contourf(ϕ_avg[2:end-1], z, ΨCorr,levels = levels, cmap=cmo.balance, 
vmin = -1, vmax = 1, extend = "both")
axs.invert_yaxis()
fig.colorbar(CM, orientation = "horizontal", fraction  =0.04)
axs.set_xticks(-40:20:60)
axs.set_xlim(-34, 60)
axs.set_ylabel("Depth [m]", fontweight = "bold")
lab = string.(abs.(collect(-40:20:60)))
lab = lab .* ["°S", "°S", "", "°N", "°N", "°N"]
axs.set_xticklabels(lab)
# axs.set_title(labels[i])
rect = patches.Rectangle((23, 2000), 59, 1000, linewidth=3, edgecolor="black",facecolor="none", alpha = 0.7)
axs.add_patch(rect)
fig
fig.savefig(plotsdir("native/paper_figures/ΔW_timemean_corr_EulnBol_mix_lowpass.png"), bbox_inches = "tight", dpi = 400)


