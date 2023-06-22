#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
using .OHC_helper
import NaNMath as nm
include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));

expname = "iter129_bulkformula"
runpath,diagpath = listexperiments(exprootdir());

filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"trsp_3d_set3")
datafilelist_R  = filter(x -> occursin("data",x),filelist)

convR = []
approx_convR = []
θFr = MeshArray(γ,Float32,50); fill!(θFr, 0.0)
θFr_conv = MeshArray(γ,Float32,50); fill!(θFr_conv, 0.0)

for tt = 1:312
    print(tt)
    dθr = γ.read(diagpath[expname]*datafilelist_R[tt],MeshArray(γ,Float32,150))
    wθ_conv = dθr[:, 1:50];
    κzθ_conv = dθr[ :, 51:100] .+ dθr[ :, 101:150]; dθr = nothing 

    θ = γ.read(diagpath[expname]*datafilelist_θ[tt],MeshArray(γ,Float32, 50))
    UV = γ.read(diagpath[expname]*datafilelist_uvw[tt],MeshArray(γ,Float32,150))
    u = UV[:, 1:50]; v = UV[:, 51:100]; w = UV[:, 101:150]
    W = w .* area
    for k = 1:49 
        θFr[:, k+1] = (θ[:, k] + θ[:, k+1]) ./2 #linear interpolate to faces
        θFr_conv[:, k+1] = θFr[:, k+1] .* W[:, k+1]
    end
    
    for ff=1:5, k=1:49
        wθ_conv.f[ff, k] .= wθ_conv.f[ff, k+1] .- wθ_conv.f[ff, k] #in - out 
        θFr_conv.f[ff, k] .= θFr_conv.f[ff, k+1] .- θFr_conv.f[ff, k] 
    end

    push!(approx_convR, sum(θFr_conv[:, 39] .* PAC_msk))
    push!(convR, sum(wθ_conv[:, 39] .* PAC_msk))
end

ctrl_vol = sum(cell_volumes[:, 39])
tecco = 1992+1/24:1/12:2018

fig, axes = plt.subplots(2, 1, figsize = (9, 5))
mean_residual = 0
axes[1].plot(tecco, convR ./ ctrl_vol, label = L"\overline{\nabla (W \theta)}", c = "blue")
axes[1].plot(tecco, approx_convR ./ ctrl_vol, label = L"\nabla (\overline{W} \overline{\theta})}", c = "orange")
axes[1].legend()
axes[1].set_ylabel("[°C / s]")
axes[1].set_title("Vertical Advective Convergence")
axes[2].plot(tecco, cumsum(convR ./ ctrl_vol) .* 2.628e+6, label = L"\int{\overline{\nabla (W \theta)} dt}", c = "blue")
axes[2].plot(tecco, cumsum((approx_convR .+ mean_residual)./ ctrl_vol) .* 2.628e+6, label = L"\int{\nabla (\overline{W} \overline{\theta})} dt}", c = "orange")
axes[2].legend()
axes[2].set_title("Vertical Advective Convergence")
axes[2].set_ylabel("[°C]")
[ax.set_xlabel("time") for ax in axes]
fig.tight_layout()
fig

