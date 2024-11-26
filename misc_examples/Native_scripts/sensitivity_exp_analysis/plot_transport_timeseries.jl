
include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    LaTeXStrings, PyCall, RollingFunctions, DSP
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");
@pyimport seaborn as sns;

colors =  sns.color_palette("deep")[1:4]
sns.set_theme(context = "talk", font_scale = 1.0,
              palette = sns.color_palette("deep"));
sns.set_style("white")

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
tecco = 1992+1/24:1/12:2018

ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "false", include_bering = true)

ϕ_avg = zonal_average(ϕ, area .* PAC_msk)
ϕ_avg = ϕ_avg[isfinite.(ϕ_avg)]

region = "NPAC"; 
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall(botlvl .< -z[:] .< uplvl)

NPAC_boundidx = Base.findmin(abs.(ϕ_avg .- 23))[2]

pad_arr(x) = vcat(reverse(x), x, reverse(x))
unpad_arr(x, arr_len) = x[arr_len+1:2*arr_len]
function low_pass_filter(signal)
    ff = digitalfilter(Lowpass(1/36, fs = 1),Butterworth(2))
    filtered_signal = filtfilt(ff, pad_arr(signal))
    return unpad_arr(filtered_signal, 312)
end

obtain_Vin(Ψ, NPAC_boundidx) = Ψ[lvls[end] + 1, NPAC_boundidx, :] .- Ψ[lvls[1], NPAC_boundidx, :]
function obtain_Wtop(Ψ, NPAC_boundidx)
    Win = Ψ[lvls[1], NPAC_boundidx, :]
    return -Win
end

function obtain_Wbot(Ψ, NPAC_boundidx)
    Wout = Ψ[lvls[end] + 1, NPAC_boundidx, :]
    return -Wout
end
open_streamfunction(expname, type) = jldopen(datadir("Ψ_" *type * "_timeseries_PAC_" * expname *".jld2"))["Ψ_exp_timeseries"]

get_face_transports(Ψ, NPAC_boundidx) = (obtain_Vin(Ψ, NPAC_boundidx), obtain_Wtop(Ψ, NPAC_boundidx), obtain_Wbot(Ψ, NPAC_boundidx))

window_size = 2
window_pad = Int(window_size/2)
extract_interannual(x) =  vcat(fill(NaN32, window_pad), 
                            rollmean(x[:], window_size), fill(NaN32, window_pad-1))
round_mean(x) = round(mean(x), digits = 2)
include(srcdir("plot_and_dir_config.jl"))

vars =  ["iter0_bulkformula", "only_init", "iter129_bulkformula", "only_kappa", "only_sfc"]
fig, ax = plt.subplots(1, 3, figsize = (20, 10))
for expname in vars 
    println(expname)
    ΨEul = open_streamfunction(expname, "Eul")
    Vin, Wout, Win = 1.f-6 .* (get_face_transports(ΨEul, NPAC_boundidx) )
    println("Mean Values: Vin ", round_mean(Vin), " Win ", round_mean(Win), " Wout ", round_mean(Wout))
    ax[1].plot(tecco, low_pass_filter(Vin), label = plot_labels[expname], c = exp_colors[expname])
    ax[2].plot(tecco, low_pass_filter(Win), label = plot_labels[expname], c = exp_colors[expname])
    ax[3].plot(tecco, -low_pass_filter(Wout), label = plot_labels[expname], c = exp_colors[expname])
end
[a.legend(frameon = false) for a in ax]
fig

fig, ax = plt.subplots(1, 3, figsize = (20, 10))
for expname in vars 
    println(expname)
    ΨEulBol = open_streamfunction(expname, "EulBol")
    Vin, Wout, Win = 1.f-6 .* (get_face_transports(ΨEulBol, NPAC_boundidx) )
    println("Mean Values: Vin ", round_mean(Vin), " Win ", round_mean(Win), " Wout ", round_mean(Wout))
    println("Mean Filtered Values: Vin ", round_mean(low_pass_filter(Vin)), " Win ", round_mean(low_pass_filter(Win)), " Wout ", round_mean(low_pass_filter(Wout)))
    ax[1].plot(tecco, low_pass_filter(Vin), label = plot_labels[expname], c = exp_colors[expname])
    ax[2].plot(tecco, low_pass_filter(Win), label = plot_labels[expname], c = exp_colors[expname])
    ax[3].plot(tecco, -low_pass_filter(Wout), label = plot_labels[expname], c = exp_colors[expname])
end
[a.legend(frameon = false) for a in ax]
fig


fig, ax = plt.subplots(1, 3, figsize = (20, 10))
for expname in vars 
    println(expname)
    ΨEulBol = open_streamfunction(expname, "EulBol")
    ΨEul = open_streamfunction(expname, "Eul")
    ΨBol = ΨEulBol .- ΨEul
    Vin, Wout, Win = 1.f-6 .* (get_face_transports(ΨBol, NPAC_boundidx) )
    println("Mean Values: Vin ", round_mean(Vin), " Win ", round_mean(Win), " Wout ", round_mean(Wout))
    println("Mean Filtered Values: Vin ", round_mean(low_pass_filter(Vin)), " Win ", round_mean(low_pass_filter(Win)), " Wout ", round_mean(low_pass_filter(Wout)))
    ax[1].plot(tecco, low_pass_filter(Vin), label = plot_labels[expname], c = exp_colors[expname])
    ax[2].plot(tecco, low_pass_filter(Win), label = plot_labels[expname], c = exp_colors[expname])
    ax[3].plot(tecco, -low_pass_filter(Wout), label = plot_labels[expname], c = exp_colors[expname])
end
[a.legend(frameon = false) for a in ax]

fig