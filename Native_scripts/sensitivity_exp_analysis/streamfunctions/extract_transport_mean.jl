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
    w_mean = MeshArray(γ, Float32, 50); fill!(w_mean, 0.0)

    @time for tt = 1:nt
        println(tt)
        fname = datafilelist[tt]
        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        u, v, w = extract_eulerian_velocities(diagpath, expname, fname, γ)
        w_mean .+= w ./nt
    end

    return w_mean
end


expname = "only_sfc"
W_PAC_sfc= get_transports(diagpath, expname, γ); 

expname = "iter0_bulkformula"
W_PAC_0 = get_transports(diagpath, expname, γ); 

expname = "only_kappa"
W_PAC_mix = get_transports(diagpath, expname, γ); 

expname = "only_init"
W_PAC_init = get_transports(diagpath, expname, γ); 



fig, ax = plt.subplots(3, 2, figsize = (21,12), sharey = true, sharex = true)
ax[1, 1].pcolormesh(tecco, -z,1e-6 .*  W_NE_0, cmap = cmo.balance, vmin = -5, vmax = 5)
ax[1, 2].pcolormesh(tecco, -z,1e-6 .*  W_NW_0, cmap = cmo.balance, vmin = -5, vmax = 5)

ax[2, 1].pcolormesh(tecco, -z,1e-6 .*  W_NE_129, cmap = cmo.balance, vmin = -5, vmax = 5)
ax[2, 2].pcolormesh(tecco, -z,1e-6 .*  W_NW_129, cmap = cmo.balance, vmin = -5, vmax = 5)

ax[3, 1].pcolormesh(tecco, -z,1e-6 .*  (W_NE_129 .- W_NE_0), cmap = cmo.balance, vmin = -5, vmax = 5)
ax[3, 2].pcolormesh(tecco, -z,1e-6 .*  (W_NW_129 .- W_NW_0), cmap = cmo.balance, vmin = -5, vmax = 5)


fig, ax = plt.subplots(3, 2, figsize = (21,12), sharey = true, sharex = true)
ax[1, 1].pcolormesh(tecco, -z,1e-6 .*  W_NE_0, cmap = cmo.balance, vmin = -5, vmax = 5)
ax[1, 2].pcolormesh(tecco, -z,1e-6 .*  W_NW_0, cmap = cmo.balance, vmin = -5, vmax = 5)

ax[2, 1].pcolormesh(tecco, -z,1e-6 .*  W_NE_kap, cmap = cmo.balance, vmin = -5, vmax = 5)
ax[2, 2].pcolormesh(tecco, -z,1e-6 .*  W_NW_kap, cmap = cmo.balance, vmin = -5, vmax = 5)

ax[3, 1].pcolormesh(tecco, -z,1e-6 .*  (W_NE_129 .- W_NE_0), cmap = cmo.balance, vmin = -5, vmax = 5)
ax[3, 2].pcolormesh(tecco, -z,1e-6 .*  (W_NW_129 .- W_NW_0), cmap = cmo.balance, vmin = -5, vmax = 5)



fig, ax = plt.subplots(3, 1, figsize = (21,12), sharey = true, sharex = true)
W_PAC_0 = W_NE_0 .+ W_NW_0
W_PAC_129 = W_NE_129 .+ W_NW_129
ax[1, 1].pcolormesh(tecco, -z,1e-6 .*  W_PAC_0, cmap = cmo.balance, vmin = -7, vmax = 7)
ax[2, 1].pcolormesh(tecco, -z,1e-6 .*  W_PAC_129, cmap = cmo.balance, vmin = -7, vmax = 7)
ax[3, 1].pcolormesh(tecco, -z,1e-6 .*  (W_PAC_129 .- W_PAC_0), cmap = cmo.balance, vmin = -7, vmax = 7)
fig
