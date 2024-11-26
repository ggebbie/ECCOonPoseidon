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

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "NPAC"; 
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
    ff = digitalfilter(Lowpass(1/24, fs = 1),Butterworth(2))
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

fig, ax = plt.subplots(figsize = (10, 5))
PC = vec(EOFS[1].U[:, 1]' * data_0); println(mean(PC))
ax.plot(tecco, PC, label = "PC 1 " * labels[1])
PC = vec(EOFS[2].U[:, 1]' * data_129); println(mean(PC))
ax.plot(tecco, PC, label = "PC 1 " * labels[2])
ax.legend()
fig

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
