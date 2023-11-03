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


reg_λ, reg_ϕ = LLCcropC(λ, γ), LLCcropC(ϕ, γ)
fig, axs = plt.subplots(1, figsize=(7,10), subplot_kw=Dict("projection"=> proj0))
axs.set_extent((120, 285, -40, 70),crs=projPC)

arr = [curlτ129 .- curlτ0]

ax = axs
data = mean(arr[1], dims = 3)[:, :, 1] .* 1e12
# data = data .* reg_mask
data[data .== 0.0] .= NaN
bounds = 4


cf = ax.pcolormesh(reg_λ, reg_ϕ,  data, 
                vmin = -bounds, vmax = bounds, shading="nearest", 
                transform=projPC, rasterized = true, cmap = cm.curl)

ax.set_title("FORCING Effect")
gl = ax.gridlines(crs=projPC, draw_labels=true, linewidth=2, 
                color="gray", alpha=0, linestyle="--")
gl.top_labels = false; gl.right_labels = false
ax.coastlines()

fig.colorbar(cf, ax=ax, orientation = "horizontal", 
fraction=0.04, pad = 0.1, extend = "both", label = "N per (10 km)³")
fig.savefig(plotsdir("native/sensitivity_exps/Windstress_curl_change.png"), dpi = 400)
fig



reg_λ, reg_ϕ = LLCcropC(λ, γ), LLCcropC(ϕ, γ)
fig, axs = plt.subplots(1, figsize=(7,10), subplot_kw=Dict("projection"=> proj0))
axs.set_extent((120, 285, -40, 70),crs=projPC)

arr = [curlτ129]

ax = axs
data = mean(arr[1], dims = 3)[:, :, 1] .* 1e12
# data = data .* reg_mask
data[data .== 0.0] .= NaN
bounds = 4


cf = ax.pcolormesh(reg_λ, reg_ϕ,  data, 
                vmin = -bounds, vmax = bounds, shading="nearest", 
                transform=projPC, rasterized = true, cmap = cm.curl)

ax.set_title("FORCING Effect")
gl = ax.gridlines(crs=projPC, draw_labels=true, linewidth=2, 
                color="gray", alpha=0, linestyle="--")
gl.top_labels = false; gl.right_labels = false
ax.coastlines()

fig.colorbar(cf, ax=ax, orientation = "horizontal", 
fraction=0.04, pad = 0.1, extend = "both", label = "N per (10 km)³")
fig.savefig(plotsdir("native/sensitivity_exps/Windstress_curl_change.png"), dpi = 400)
fig

#do an EOF


fig, axs = plt.subplots(1, figsize=(7,10), subplot_kw=Dict("projection"=> proj0))
axs.set_extent((120, 285, -40, 70),crs=projPC)

arr = [τx_129 .- τx_0]

ax = axs
data = mean(arr[1], dims = 3)[:, :, 1] .* 1e2
data = data .* reg_mask
data[data .== 0.0] .= NaN
bounds = 1


cf = ax.pcolormesh(reg_λ, reg_ϕ,  data, 
                vmin = -bounds, vmax = bounds, shading="nearest", 
                transform=projPC, rasterized = true, cmap = cm.curl)

ax.set_title("FORCING Effect")
gl = ax.gridlines(crs=projPC, draw_labels=true, linewidth=2, 
                color="gray", alpha=0, linestyle="--")
gl.top_labels = false; gl.right_labels = false
ax.coastlines()

fig.colorbar(cf, ax=ax, orientation = "horizontal", 
fraction=0.04, pad = 0.1, extend = "both", label = "N per (10 km)³")
fig.savefig(plotsdir("native/sensitivity_exps/Windstress_curl_change.png"), dpi = 400)
fig