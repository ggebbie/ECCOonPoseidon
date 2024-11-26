include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall
import PyPlot as plt 
include(srcdir("config_exp.jl"))

include(srcdir("MeshArraysPlots.jl"))

tecco = 1992+1/24:1/12:2018
nz = 50

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
λ_wrap = wrap_λ(λ)
ocean_mask = wet_pts(Γ)

region = "NPAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

include(srcdir("plot_and_dir_config.jl"))

NW_PAC = region_mask(ocean_mask, λ_wrap, ϕ, (0, 90), (0, 190))
NE_PAC = region_mask(ocean_mask, λ_wrap, ϕ, (0, 90), (190, 360))

cs, sn = get_cs_and_sn(γ)


function get_transports(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    fileroot = "trsp_3d_set1"
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(datafilelist);
    w_mean = MeshArray(γ, Float32, 50); fill!(w_mean, 0.0)
    N_mean = MeshArray(γ, Float32, 50); fill!(N_mean, 0.0)

    @time for tt = 1:nt
        println(tt)
        fname = datafilelist[tt]
        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        u, v, w = extract_eulerian_velocities(diagpath, expname, fname, γ)
        w_mean .+= w    / nt

        N_mean .+= Nθ    / nt
        
    end

    return w_mean, N_mean
end


# utrsp, vtrsp = UVtoTransport(u, v, Γ)
# Eθ, Nθ = rotate_UV_native(utrsp, vtrsp, cs, sn) 
        
expname = "only_sfc"
W_mean_sfc, V_mean_sfc= get_transports(diagpath, expname, γ); 

expname = "only_kappa"
W_mean_mix, V_mean_mix= get_transports(diagpath, expname, γ); 

expname = "only_init"
W_mean_init, V_mean_init= get_transports(diagpath, expname, γ); 

expname = "iter0_bulkformula"
W_mean_0, V_mean_0 = get_transports(diagpath, expname, γ); 

fig, axs = plt.subplots(1, 3, figsize=(25,10), subplot_kw=Dict("projection"=> proj0))

force_eff = W_mean_sfc .- W_mean_0; force_eff = force_eff[:, 38] .* area
mix_eff = W_mean_mix .- W_mean_0; mix_eff = mix_eff[:, 38] .* area
init_eff = W_mean_init .- W_mean_0; init_eff = init_eff[:, 38] .* area

vmax = maximum(mix_eff .* PAC_msk) /4

for ff = 1:5
    force_eff.f[ff][force_eff.f[ff] .== 0.0] .= NaN
    mix_eff.f[ff][mix_eff.f[ff] .== 0.0] .= NaN
    init_eff.f[ff][init_eff.f[ff] .== 0.0] .= NaN
end
[a.set_extent((120, 285, -70, 56),crs=projPC) for a in axs]

DXCsm=2*Γ.DXC; DYCsm=2*Γ.DYC;


pcolormesh_ma(axs[1], λ, ϕ, MeshArrays.smooth(init_eff,DXCsm,DYCsm,Γ), 
cmo.balance; bounds = (-vmax, vmax))
pcolormesh_ma(axs[2], λ, ϕ, MeshArrays.smooth(mix_eff,DXCsm,DYCsm,Γ), 
cmo.balance; bounds = (-vmax, vmax))
pcolormesh_ma(axs[3], λ, ϕ, MeshArrays.smooth(force_eff,DXCsm,DYCsm,Γ), 
cmo.balance; bounds = (-vmax, vmax))

fig