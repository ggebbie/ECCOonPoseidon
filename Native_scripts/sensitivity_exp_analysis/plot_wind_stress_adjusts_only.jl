#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../../src/intro.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LinearAlgebra, LaTeXStrings, PyCall
import PyPlot as plt
import NaNMath as nm
@pyimport matplotlib.animation as anim

cm = pyimport("cmocean.cm");
@pyimport seaborn as sns;
include(srcdir("MeshArraysPlots.jl"))

sns.set_theme(context = "notebook", style = "ticks",
              palette = sns.color_palette("colorblind"));
              
include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)

region = "PAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
ocean_mask = wet_pts(Γ)

tecco= 1992+1/24:1/12:2018 # ecco years
nt = length(tecco)

runpath,diagpath = listexperiments(exprootdir());
include(srcdir("plot_and_dir_config.jl"))

τx_dict = Dict(); τy_dict = Dict(); curlτ_dict = Dict(); ekup_dict = Dict()
vars = ["iter0_bulkformula", "iter129_bulkformula", "only_init", "only_kappa", "only_sfc"]

for expname in vars
    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    τx_dict[expname] = MeshArray(γ,Float32); fill!(τx_dict[expname], 0.0)
    τy_dict[expname] = MeshArray(γ,Float32); fill!(τy_dict[expname], 0.0)
    curlτ_dict[expname] = MeshArray(γ,Float32); fill!(curlτ_dict[expname], 0.0)
    for tt = 1:nt
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        Tname = τdatafilelist[tt]
        
        τx, τy = extract_ocnTAU(diagpath, expname , τdatafilelist[tt], γ)
        τcurl = MeshArrays.curl(τx, τy, Γ)
        curlτ_dict[expname].+= τcurl ./ nt
    end
    # curlτ_dict[expname] = MeshArrays.curl(τx_dict[expname], τy_dict[expname], Γ)
    # τxC, τyC = velocity2center(τx_dict[expname], τy_dict[expname], Γ)
    # τE, τN = rotate_uv(τxC, τyC, Γ)
    # τx_dict[expname] = τE
    # τy_dict[expname] = τN
end


# proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
proj0 = ECCOonPoseidon.cartopy.crs.Robinson(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
# λ_wrap = wrap_λ(λ)

eff_exps = ["only_init", "only_kappa", "only_sfc"]
curlτ_dicteff = deepcopy(curlτ_dict)
[curlτ_dicteff[expt] .-= curlτ_dict["iter0_bulkformula"] for expt in eff_exps]

#setting up plotting stuff      
#plot τx
fig, ax = plt.subplots(5, figsize=(15,20), 
                       subplot_kw=Dict("projection"=> proj0), sharey = true)
[a.coastlines() for a in ax]
[a.set_extent((120, 240, 20, 56),crs=projPC) for a in ax]'
CF = Any[]
for (i, var) in enumerate(keys(curlτ_dicteff))
    for ff = 1:5
        data = curlτ_dicteff[var].f[ff]; data[data .== 0.0] .= NaN
        cf = ax[i].pcolormesh(λ_wrap.f[ff], ϕ.f[ff], data, 
    cmap = cm.curl, transform = projPC, vmin = -1e-12, vmax = 1e-12)
    push!(CF, cf)
    end
    ax[i].set_title(var)

end
# ax.set_extent((-180, 180, -70, 56),crs=projPC)
fig.colorbar(CF[1], ax = ax[:], orientation = "horizontal")
fig.savefig(plotsdir("native/sensitivity_exps/wind_adjusts.png"), bbox_inches = "tight", dpi = 400);
sum(curlτ_dict["only_init"])
sum(curlτ_dict["only_sfc"])
sum(curlτ_dict["only_kappa"])