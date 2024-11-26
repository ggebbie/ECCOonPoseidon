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
    w_mean = MeshArray(γ, Float32); fill!(w_mean, 0.0)

    @time for tt = 1:nt
        println(tt)
        fname = datafilelist[tt]
        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        u, v, w = extract_eulerian_velocities(diagpath, expname, fname, γ, 43)
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

expname = "iter129_bulkformula"
W_PAC_129 = get_transports(diagpath, expname, γ); 

# proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
W_Dict = Dict() 
W_Dict["iter0_bulkformula"] = W_PAC_0 .* 1
W_Dict["iter129_bulkformula"] = W_PAC_129 .* 1
W_Dict["only_init"] = W_PAC_init .* 1
W_Dict["only_kappa"] = W_PAC_mix .* 1
W_Dict["only_sfc"] = W_PAC_sfc .* 1

proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
# λ_wrap = wrap_λ(λ)

eff_exps = ["only_init", "only_kappa", "only_sfc"]
W_dicteff = deepcopy(W_Dict)
[W_dicteff[expt] .-= W_Dict["iter0_bulkformula"] for expt in eff_exps]
W_dicteff["SUM"] = 1. * W_Dict["iter0_bulkformula"] 
[W_dicteff["SUM"] .+= W_dicteff[expt] for expt in eff_exps]

#setting up plotting stuff      
#plot τx
fig, ax = plt.subplots(2, 3, figsize=(20,7.5), 
                       subplot_kw=Dict("projection"=> proj0), sharey = true)
# [a.coastlines() for a in ax]
# [a.set_extent((120, 240, 20, 56),crs=projPC) for a in ax]
CF = Any[]

for (i, var) in enumerate(keys(W_dicteff))
    for ff = 1:5
        data = W_dicteff[var].f[ff] .* 100 * 86400 ; data[data .== 0.0] .= NaN
        cf = ax[i].pcolormesh(λ_wrap.f[ff], ϕ.f[ff], data, shading="nearest",
    cmap = cm.balance, transform = projPC, vmin = -15, vmax = 15, rasterized = true)
    push!(CF, cf)
    end
    ax[i].set_title(var)
    ax[i].coastlines(resolution="110m")
    ax[i].set_extent((120, 260, 17, 56),crs=projPC)
    # ax[i].set_title(plot_labels_effect[var])
    gl = ax[i].gridlines(crs=projPC, draw_labels=true,
                        linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.top_labels = false
    gl.right_labels = false
end
fig
# ax.set_extent((-180, 180, -70, 56),crs=projPC)
fig.colorbar(CF[1], ax = ax[:], orientation = "horizontal", 
fraction = 0.04, label = "cm per day")
fig
fig.savefig(plotsdir("native/sensitivity_exps/velocity3000_adjusts.png"), bbox_inches = "tight", dpi = 400);
sum(curlτ_dict["only_init"])
sum(curlτ_dict["only_sfc"])
sum(curlτ_dict["only_kappa"])

#points of interest for only_kappa
iff, ixxiyy, _ = findlatlon(λ, ϕ, -154, 45)
ax[4].scatter(λ[iff][ixxiyy], ϕ[iff][ixxiyy], transform = projPC, c= "green")
iff, ixxiyy, _ = findlatlon(λ, ϕ, -150, 47)
ax[4].scatter(λ[iff][ixxiyy], ϕ[iff][ixxiyy], transform = projPC, c= "green")
iff, ixxiyy, _ = findlatlon(λ, ϕ, 150, 37)
ax[4].scatter(λ[iff][ixxiyy], ϕ[iff][ixxiyy], transform = projPC, c= "green")

fig