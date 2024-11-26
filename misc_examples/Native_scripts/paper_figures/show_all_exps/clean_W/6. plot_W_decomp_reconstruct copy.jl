include("../../../../src/intro.jl")

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

pad_arr(x) = vcat(reverse(x), x, reverse(x))
unpad_arr(x, arr_len) = x[arr_len+1:2*arr_len]
function low_pass(signal)
    nt = length(signal)
    signal_mean = mean(signal)
    ff = digitalfilter(Lowpass(1/(8*12 + 1), fs = 1),Butterworth(4))
    filtered_signal = filtfilt(ff, pad_arr(signal .- signal_mean))
    return unpad_arr(filtered_signal, nt) .+ signal_mean
end

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
λ_wrap = wrap_λ(λ)
ocean_mask = wet_pts(Γ)
cmo = pyimport("cmocean.cm")

region = "PAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)

region = "NPAC"
NPAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

maskW = 1 .* NPAC_msk
[maskW.f[i][iszero.(maskW.f[i])] .= NaN for i in 1:5]
[maskW.f[i][isfinite.(maskW.f[i])] .= 1 for i in 1:5]

include(srcdir("plot_and_dir_config.jl"))

W_mean = jldopen(datadir("W_2000_all_times.jld2"))["W"]

W_mean["diff_wind"] = W_mean["only_wind"] .- W_mean["iter0_bulkformula"]
W_mean["diff_init"] = W_mean["only_init"] .- W_mean["iter0_bulkformula"]
W_mean["diff_kappa"] = W_mean["only_kappa"] .- W_mean["iter0_bulkformula"]
years = [12*5, 12*15, 12*25]; nyears = length(years)
filt_W_dict =  Dict();
for expname in ["diff_wind", "diff_init", "diff_kappa"]
    filt_W = 1 .* W_mean[expname]
    for i in 1:5
        println(i)
        tmp = filt_W[i, :]
        tmp = cat(tmp..., dims = 3)
        nx, ny = size(tmp[:, :, 1])
        for ix in 1:nx
            for iy in 1:ny
                low_pass_tmp = low_pass(tmp[ix, iy, :])
                tmp[ix, iy, :] .= low_pass_tmp
            end
        end
        tmp1 = [tmp[:, :, i] for i in 1:312]
        filt_W.f[i, :] .= tmp1
    end
    
    mean_W = W_mean[expname][:, 1:nyears] .* 0.0
    fill!(mean_W, 0.0)
    nmax = 12
    for (j, y) in enumerate(years)
        TMP = filt_W[:, y:(y+nmax)]
        for i in 1:nmax
            tmp = (mean_W.f[:, j] .+ TMP.f[:, i]) ./ nmax
            mean_W.f[:, j] .= tmp
        end
    end
    filt_W_dict[expname] = mean_W
end

vmaxs = [1e3, 500, 1e3]
# vmaxs = 1e-6 .* [1e3, 1e2, 1e3]

fig, axs = plt.subplots(3, 3, figsize=(14,10), subplot_kw=Dict("projection"=> proj0))
CF = Any[]
for (i, expname) in enumerate(["diff_wind", "diff_init", "diff_kappa"])
    ax = axs[:, i]
    mean_W = filt_W_dict[expname]
    vmax = vmaxs[i]
    for (j, y) in enumerate(years)
        data = mean_W[:, j] .* maskW
        data = data .* area
        print(vmax)
        for ff = 1:5
            cf = ax[j].pcolormesh(λ[ff], ϕ[ff],  data[ff], transform=projPC, 
            cmap = cmo.balance, vmin = -vmax, vmax = vmax, rasterized = true, shading = "nearest") 
            push!(CF, cf)
        end
        ax[j].set_extent((120, 260, 23, 65),crs=projPC)
        ax[j].coastlines()
    end
end
fig
fig.colorbar(CF[1], ax = axs[end, 1], orientation  ="horizontal", label = L"\Delta w_{2000}" * " " *L" [m^3 / s]")
fig.colorbar(CF[16], ax = axs[end, 2], orientation  ="horizontal", label = L"\Delta w_{2000}" * " " *L" [m^3 / s]")
fig.colorbar(CF[44], ax = axs[end, 3], orientation  ="horizontal", label = L"\Delta w_{2000}" * " " *L" [m^3 / s]")

fig.tight_layout()

[a.set_extent((120, 255, 23, 65),crs=projPC) for a in axs[:]]
[a.set_aspect(1.3) for a in axs[:]]

fig.subplots_adjust(hspace = -0.6)

for a in axs[:, 1][:]
    gl = a.gridlines(crs=projPC, draw_labels=true,
                        linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.top_labels = false
    gl.right_labels = false
    gl.bottom_labels = false
end

for a in axs[end, :][:]
    gl = a.gridlines(crs=projPC, draw_labels=true,
                        linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.left_labels = false
    gl.right_labels = false
    gl.top_labels = false
end

years_text = ["1996", "2006", "2016"]
exp_text = ["Wind Stress\nAdjustment Effect", "Initial Condition\nAdjustment Effect", "Mixing Parameter\nAdjustment Effect"]
for i in 1:3
    a = axs[i, 1]
    a.annotate(years_text[i], xy=[-0.3, 0.5], 
    transform=a.transAxes,
    xycoords = a.transAxes,
    rotation = 90, 
    fontweight = "bold", 
    ha="center")
    a = axs[1, i]
    a.annotate(exp_text[i], xy=[0.5, 1.1], 
    transform=a.transAxes,
    xycoords = a.transAxes,
    rotation = 0, 
    fontweight = "bold", 
    ha="center")
end

fig.savefig(plotsdir("native/paper_figures/ΔW_maps.png"), bbox_inches = "tight", dpi = 400)
