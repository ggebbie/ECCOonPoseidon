include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall
import PyPlot as plt 
include(srcdir("config_exp.jl"))
include(srcdir("MeshArraysPlots.jl"))

tecco = 1992+1/24:1/12:2018
nz = 50

@pyimport cmocean.cm as cmo


(ϕ,λ) = latlonC(γ)
area = readarea(γ)
λ_wrap = wrap_λ(λ)
ocean_mask = wet_pts(Γ)

region = "NPAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

include(srcdir("plot_and_dir_config.jl"))

area_mask = area .* PAC_msk

cs, sn = get_cs_and_sn(γ)

function get_transports(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    fileroot = "trsp_3d_set1"
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(datafilelist);
    We_PAC = zeros(50, nt)
    Web_PAC = zeros(50, nt)

    @time for tt = 1:nt
        println(tt)
        fname = datafilelist[tt]
        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        u, v, w = extract_eulerian_velocities(diagpath, expname, fname, γ)
        We_PAC[:, tt] .= lateral_sum(w .* area_mask)
    end

    return We_PAC
end


expname = "only_sfc"
W_NE_sfc = get_transports(diagpath, expname, γ); 
expname = "iter0_bulkformula"
W_NE_0= get_transports(diagpath, expname, γ); 

expname = "only_kappa"
W_NE_kap = get_transports(diagpath, expname, γ); 

W_kap_eff = W_NE_kap .- W_NE_0
W_sfc_eff = W_NE_sfc .- W_NE_0

fig, ax = plt.subplots(3, 1, figsize = (21,12), sharey = true, sharex = true)
ax[1].pcolormesh(tecco, -z,1e-6 .*  W_NE_0, cmap = cmo.balance, vmin = -5, vmax = 5)
ax[1].set_title("CTRL")
ax[2].pcolormesh(tecco, -z,1e-6 .*  W_kap_eff, cmap = cmo.balance, vmin = -5, vmax = 5)
ax[2].set_title("MIXING Effect")
ax[ 3].pcolormesh(tecco, -z,1e-6 .*  W_sfc_eff, cmap = cmo.balance, vmin = -5, vmax = 5)
ax[3].set_title("SFC Effect")
# [a.set_ylim(-5000, 0) for a in ax]
[a.set_ylabel("Depth [m]") for a in ax]
fig.suptitle("North Pacific Vertical Transport")
fig