include("../../src/intro.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LinearAlgebra, LaTeXStrings, PyCall, DataFrames, CSV, DSP
import PyPlot as plt
import NaNMath as nm

@pyimport matplotlib.animation as anim

cm = pyimport("cmocean.cm");
@pyimport seaborn as sns;

sns.set_theme(context = "notebook", style = "ticks",
              palette = sns.color_palette("colorblind"));
              
include(srcdir("config_exp.jl"))
runpath,diagpath = listexperiments(exprootdir());
include(srcdir("plot_and_dir_config.jl"))
tecco = 1992+1/24:1/12:2018; nz = 50

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "PAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
reg_mask = LLCcropC(PAC_msk,γ)

τ = jldopen(datadir("native/oceτ_timeseries_regular_grid.jld2"))["τ"]

(τx_129, τy_129, curlτ129) = τ["iter129_bulkformula"]
(τx_0, τy_0, curlτ0) = τ["iter0_bulkformula"]

nx, ny, nt = size(curlτ129)

pad_arr(x) = vcat(reverse(x), x, reverse(x))
unpad_arr(x, arr_len) = x[arr_len+1:2*arr_len]
function low_pass_filter(signal)
    # Create a low-pass filter
    signal_nomean = signal .- mean(signal)
    ff = digitalfilter(Lowpass(1/36, fs = 1),Butterworth(2))
    # Apply the filter to the signal
    filtered_signal = filtfilt(ff, pad_arr(signal_nomean))

    return unpad_arr(filtered_signal, 312)
end

function compute_EOF(data)
    Y = reshape(data .* reg_mask, nx*ny, nt); 
    Y .-= mean(Y, dims = 2); wherenan = isnan.(Y)
    Y[wherenan] .= 0.0
    nl = size(Y, 1)
    for ll in 1:nl
        Y[ll, :] .= low_pass_filter(Y[ll, :]) #apply low-pass filter before EOF
    end

    return Y, svd(Y; full = false)
end

data_129, EOF_129 = compute_EOF(τx_129)
data_0, EOF_0 = compute_EOF(τx_0)
EOFS = [EOF_0, EOF_129]

labels = ["Iteration 0", "Iteration 129"]

reg_λ, reg_ϕ = LLCcropC(λ, γ), LLCcropC(ϕ, γ)
fig, axs = plt.subplots(2, figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
signs = [1, -1]
for i = 1:2
    ax = axs[i]
    arr = reshape(EOFS[i].U[:, 1], nx, ny)
    bounds = nm.maximum(abs.(arr)) 


    cf = ax.pcolormesh(reg_λ, reg_ϕ,  signs[i] * arr, 
                    vmin = -bounds, vmax = bounds, shading="nearest", 
                    transform=projPC, rasterized = true, cmap = cm.curl)

    ax.set_title(labels[i])
    gl = ax.gridlines(crs=projPC, draw_labels=true, linewidth=2, 
                    color="gray", alpha=0, linestyle="--")
    gl.top_labels = false; gl.right_labels = false
    ax.coastlines()
end
fig

fig.colorbar(cf, ax=ax, orientation = "horizontal", 
fraction=0.04, pad = 0.1, extend = "both", label = "N per m²")
#do an EOF

ϕ_avg = zonal_average(ϕ, area .* PAC_msk)
ϕ_avg = ϕ_avg[isfinite.(ϕ_avg)]

NPAC_boundidx = Base.findmin(abs.(ϕ_avg .- 23))[2]
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall(botlvl .< -z[:] .< uplvl)

function obtain_Wbot(Ψ, NPAC_boundidx)
    Wout = Ψ[lvls[end] + 1, NPAC_boundidx, :]
    return -Wout
end
open_streamfunction(expname, type) = jldopen(datadir("Ψ_" *type * "_timeseries_PAC_" * expname *".jld2"))["Ψ_exp_timeseries"]

expname = "only_sfc"
ΨEul = open_streamfunction(expname, "Eul")
Wbot = obtain_Wbot(ΨEul, NPAC_boundidx)

expname = "iter0_bulkformula"
ΨEul = open_streamfunction(expname, "Eul")
Wbot0 = obtain_Wbot(ΨEul, NPAC_boundidx)


fig, ax = plt.subplots(2, 1, figsize = (10, 5))
PC = vec(EOFS[1].U[:, 1]' * data_0); println(mean(PC))
ax[1].plot(tecco, PC, label = "PC 1 " * labels[1])
PC = vec(EOFS[2].U[:, 1]' * data_129); println(mean(PC))
ax[1].plot(tecco, PC, label = "PC 1 " * labels[2])
ax[1].legend()

ax[2].plot(tecco[12*3:end], (low_pass_filter(Wbot0) .+ mean(Wbot0))[12*3:end], label = "CTRL")
ax[2].plot(tecco[12*3:end], (low_pass_filter(Wbot) .+ mean(Wbot))[12*3:end], label = "FORCING")


ax[2].legend()
fig