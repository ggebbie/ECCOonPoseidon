#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LinearAlgebra, LaTeXStrings, PyCall
import PyPlot as plt
import NaNMath as nm
@pyimport matplotlib.animation as anim

cm = pyimport("cmocean.cm");
@pyimport seaborn as sns;

sns.set_theme(context = "notebook", style = "ticks",
              palette = sns.color_palette("colorblind"));
              
include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)

ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

region = "PAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
ocean_mask = OHC_helper.wet_pts(Γ)

tecco= 1992+1/24:1/12:2018 # ecco years
nt = length(tecco)

runpath,diagpath = listexperiments(exprootdir());
diagpath["clim_tau_iter0"] = "/vast/ECCOv4r4/exps/clim_tau_iter0/run/diags/"

τx_dict = Dict(); τy_dict = Dict(); 
curlτ_dict = Dict(); ekup_dict = Dict()
Qnet_dict = Dict()
vars = ["iter129_bulkformula", "iter0_bulkformula"]

println(collect(keys(test)))
for expname in vars
    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    forcingdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    τx_dict[expname] = MeshArray(γ,Float32); fill!(τx_dict[expname], 0.0)
    τy_dict[expname] = MeshArray(γ,Float32); fill!(τy_dict[expname], 0.0)
    Qnet_dict[expname] = MeshArray(γ,Float32); fill!(Qnet_dict[expname], 0.0)

    for tt = 1:nt
        forcename = forcingdatafilelist[tt]
        forcings = mdsio2dict(diagpath[expname],forcename[1:end-5],γ)
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)

        for ff = 1:5
            Qnet_dict[expname].f[ff] .+= forcings["oceQnet"].f[ff] ./ nt
            τx_dict[expname].f[ff] .+= forcings["oceTAUX"].f[ff] ./ nt
            τy_dict[expname].f[ff] .+= forcings["oceTAUY"].f[ff] ./ nt
        end
    end
    curlτ_dict[expname] = MeshArrays.curl(τx_dict[expname], τy_dict[expname], Γ)
    τxC, τyC = velocity2center(τx_dict[expname], τy_dict[expname], Γ)
    τE, τN = rotate_uv(τxC, τyC, Γ)
    τx_dict[expname] = τE
    τy_dict[expname] = τN
end

diff_dict(d, x, y) = d[x] .- d[y]
div_dict(d, x, y) = d[x] ./ d[y]
# shortnames["climatological_tau"] = "mean_wind_i129"
# shortnames["clim_tau_iter0"] = "mean_wind_i0"

diff_label = shortnames[vars[1]] * " minus " * shortnames[vars[2]]
curlτ_dict[diff_label] = diff_dict(curlτ_dict, vars[1], vars[2])
τx_dict[diff_label] = diff_dict(τx_dict, vars[1], vars[2])
τy_dict[diff_label] = diff_dict(τy_dict, vars[1], vars[2])
Qnet_dict[diff_label] = diff_dict(Qnet_dict, vars[1], vars[2])
push!(vars, diff_label) 
Qnet_dict[diff_label] = 100 .* div_dict(Qnet_dict, vars[3], vars[2])
curlτ_dict[diff_label] = 100 .* div_dict(curlτ_dict, vars[3], vars[2])
for ff in 1:5
    Qnet_dict[diff_label][ff][(!isfinite).(Qnet_dict[diff_label][ff])] .= NaN
    curlτ_dict[diff_label][ff][(!isfinite).(curlτ_dict[diff_label][ff])] .= NaN
end
# proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
proj0 = ECCOonPoseidon.cartopy.crs.Robinson(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
λ_wrap = OHC_helper.wrap_λ(λ)

#setting up plotting stuff       
bounds = [0.15, 0.15, 0.04]
CF = Any[1, 1, 1]

#plot curl
fig, ax = plt.subplots(1, 3, figsize=(15,20), 
                       subplot_kw=Dict("projection"=> proj0), sharey = true)
[a.coastlines() for a in ax]
[a.set_extent((120, 285, -70, 56),crs=projPC) for a in ax]

bounds = [1e-12, 1e-12, 200]
for (i, var) in enumerate(vars)
    CF[i] = OHC_helper.pcolormesh_ma(ax[i], λ_wrap, ϕ, curlτ_dict[var], 
                        bounds[i], cm.curl, projPC, true)
    ax[i].set_title(var)
end

fig.colorbar(CF[1], ax = ax[1:2], fraction = 0.01, pad = 0.03, orientation = "horizontal", label = "[N/m³]")
fig.colorbar(CF[3], ax = ax[3], fraction = 0.01,  pad = 0.03, orientation = "horizontal", label = "[N/m³]")
fig.suptitle("Mean Wind Stress Curl [1992 - 2017]", y = 0.31)
# fig.savefig(plotsdir("native/mean_curlτ_" * shortnames[vars[1]] * "_" * shortnames[vars[2]] * ".png"), bbox_inches = "tight", dpi = 400);
fig


#plot Qnet
fig, ax = plt.subplots(1, 3, figsize=(15,20), 
                       subplot_kw=Dict("projection"=> proj0), sharey = true)
[a.coastlines() for a in ax]
[a.set_extent((120, 285, -70, 56),crs=projPC) for a in ax]
bounds = [50, 50, 200]

for (i, var) in enumerate(vars)
    CF[i] = OHC_helper.pcolormesh_ma(ax[i], λ_wrap, ϕ, Qnet_dict[var], 
                        bounds[i], cm.balance, projPC, true)
    ax[i].set_title(var)
end

fig.colorbar(CF[1], ax = ax[1:2], fraction = 0.01, pad = 0.03, orientation = "horizontal", label = "[W/m²]")
fig.colorbar(CF[3], ax = ax[3], fraction = 0.01,  pad = 0.03, orientation = "horizontal", label = "% Difference")
fig.suptitle("Mean Qnet [1992 - 2017]", y = 0.31)
# fig.savefig(plotsdir("native/mean_curlτ_" * shortnames[vars[1]] * "_" * shortnames[vars[2]] * ".png"), bbox_inches = "tight", dpi = 400);
fig