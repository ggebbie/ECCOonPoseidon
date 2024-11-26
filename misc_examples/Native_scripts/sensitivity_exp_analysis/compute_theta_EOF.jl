include("../../src/intro.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LinearAlgebra, LaTeXStrings, PyCall, DataFrames, CSV, DSP
import PyPlot as plt
import NaNMath as nm

using CSV

PDO = CSV.read("/home/ameza/ECCOonPoseidon/Native_scripts/winds/annpdo.csv",  DataFrame,
header=false)

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
function extract_tau(expname, γ)
    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(τdatafilelist)
    tauX = zeros(size(reg_mask)..., nt)
    tauY = zeros(size(reg_mask)..., nt)

    for tt = 1:nt
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        Tname = τdatafilelist[tt]
        
        τx, τy = extract_ocnTAU(diagpath, expname , τdatafilelist[tt], γ)
        τxC, τyC = velocity2center(τx, τy, Γ)
        τE, τN = rotate_uv(τxC, τyC, Γ)
        tauX[:, :, tt] .= LLCcropC(τE, γ)
        tauY[:, :, tt] .= LLCcropC(τN, γ)
    end
    return tauX, tauY
end

function low_pass_filter(signal)
    # Create a low-pass filter
    ff = digitalfilter(Lowpass(1/36, fs = 1),Butterworth(3))
    # Apply the filter to the signal
    filtered_signal = filt(ff, signal)

    return filtered_signal
end


expname = "iter129_bulkformula"
τx_129, τy_129 = extract_tau(expname, γ);

expname = "iter0_bulkformula"
τx_0, τy_0 = extract_tau(expname, γ);

nx, ny, nt = size(τx)

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

tecco_yearly = 1992:2017
yearly_129, τx_129_EOF = compute_EOF(τx)
yearly_0, τx_0_EOF = compute_EOF(τx_0)
labels = ["Iteration 0", "Iteration 129"]
normalize = (x .- mean(x)) ./ std(x)
fig, ax = plt.subplots(figsize = (10, 5))
EOF = τx_0_EOF.U[:, 1]
PC = vec(EOF' * yearly_0); println(mean(PC))
ax.plot(tecco_yearly, -PC, label = "PC 1 " * labels[1])

EOF = τx_129_EOF.U[:, 1]
PC = vec(EOF' * yearly_129); println(mean(PC))
ax.plot(tecco_yearly, -PC, label = "PC 1 " * labels[2])

ax.legend()
fig


fig, axs = plt.subplots(2, figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
EOFS = [τx_0_EOF, τx_129_EOF]
for i = 1:2
    ax = axs[i]
    arr = reshape(EOFS[i].U[:, 1], nx, ny)
    bounds = nm.maximum(abs.(arr)) 


    cf = ax.pcolormesh(reg_λ, reg_ϕ,  -arr, 
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
