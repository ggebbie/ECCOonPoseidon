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
area = readarea(γ)

ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

region = "PAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
ocean_mask = OHC_helper.wet_pts(Γ)
cell_depths = OHC_helper.get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));
H = OHC_helper.smush(cell_depths[:, 1:38], γ); H[findall(H .== 0.0)] = Inf

vertical_average(ds) = OHC_helper.depth_average(ds[:, 1:38], cell_depths[:, 1:38], H, γ)

tecco= 1992+1/24:1/12:2018 # ecco years
nt = length(tecco)

runpath,diagpath = listexperiments(exprootdir());
τx_dict = Dict(); τy_dict = Dict(); curlτ_dict = Dict(); ekup_dict = Dict()
vars = ["iter129_bulkformula", "iter0_bulkformula"]

#plot ekman pumping velocity
f(ϕ) = abs(ϕ) > 1 ? (2 * (2π / 86400) * sind(ϕ)) : Inf32
fs = MeshArray(γ,Float32); fs.f .= map.(x -> f.(x), ϕ.f)
finv = 1 ./ fs; finv = Float32.(finv) 
σ2 = jldopen(datadir("native/native_sigma2_timeavg_1992_2017.jld2"))["ρθSavg"]
σ2 = Dict(k => vertical_average(v["σ2"]) for (k, v) in σ2)
inv_σ2 = Dict()
for (k, v) in σ2
    inv_σ2[k] = deepcopy(σ2[k]);
    for ff in eachindex(inv_σ2[k])
        inv_σ2[k].f[ff][inv_σ2[k].f[ff] .== 0.0] .= Inf
        inv_σ2[k].f[ff] .= 1 ./ inv_σ2[k].f[ff]
    end
end

for expname in vars
    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    τx_dict[expname] = MeshArray(γ,Float32); fill!(τx_dict[expname], 0.0)
    τy_dict[expname] = MeshArray(γ,Float32); fill!(τy_dict[expname], 0.0)

    for tt = 1:nt
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        Tname = τdatafilelist[tt]
        
        τx, τy = OHC_helper.extract_ocnTAU(diagpath, expname , τdatafilelist, tt, γ)

        for ff = 1:5
            τx_dict[expname].f[ff] .+= τx.f[ff] ./ nt
            τy_dict[expname].f[ff] .+= τy.f[ff] ./ nt
        end
    end
    curlτ_dict[expname] = MeshArrays.curl(τx_dict[expname], τy_dict[expname], Γ)
    fσinv = inv_σ2[expname] .* finv
    ekup_dict[expname] =  MeshArrays.curl(τx_dict[expname] .* fσinv, τy_dict[expname] .* fσinv, Γ)

    τxC, τyC = velocity2center(τx_dict[expname], τy_dict[expname], Γ)
    τE, τN = rotate_uv(τxC, τyC, Γ)
    τx_dict[expname] = τE
    τy_dict[expname] = τN
end

diff_dict(d, x, y) = d[x] .- d[y]
diff_label = shortnames[vars[1]] * " minus " * shortnames[vars[2]]
curlτ_dict[diff_label] = diff_dict(curlτ_dict, vars[1], vars[2])
τx_dict[diff_label] = diff_dict(τx_dict, vars[1], vars[2])
τy_dict[diff_label] = diff_dict(τy_dict, vars[1], vars[2])
ekup_dict[diff_label] = diff_dict(ekup_dict, vars[1], vars[2])

push!(vars, diff_label) 

# proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
proj0 = ECCOonPoseidon.cartopy.crs.Robinson(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
λ_wrap = OHC_helper.wrap_λ(λ)

#setting up plotting stuff       
bounds = [0.15, 0.15, 0.04]
CF = Any[1, 1, 1]

#plot τx
fig, ax = plt.subplots(1, 3, figsize=(15,20), 
                       subplot_kw=Dict("projection"=> proj0), sharey = true)

for (i, var) in enumerate(vars)
    CF[i] = OHC_helper.pcolormesh_ma(ax[i], λ_wrap, ϕ, τx_dict[var], 
                        bounds[i], cm.curl, projPC, true)
    ax[i].set_title(var)
end
# ax.set_extent((-180, 180, -70, 56),crs=projPC)

[a.coastlines() for a in ax]
[a.set_extent((-180, 180, -70, 56),crs=projPC) for a in ax]
fig.colorbar(CF[1], ax = ax[1:2], fraction = 0.01, pad = 0.03, orientation = "horizontal", label = "[N/m²]")
fig.colorbar(CF[3], ax = ax[3], fraction = 0.01,  pad = 0.03, orientation = "horizontal", label = "[N/m²]")
fig.suptitle("Mean Zonal Wind Stress [1993 - 2017]", y = 0.25)
fig.savefig(plotsdir("native/mean_τx_" * shortnames[vars[1]] * "_" * shortnames[vars[2]] * ".png"), bbox_inches = "tight", dpi = 400);

#plot τy
fig, ax = plt.subplots(1, 3, figsize=(15,20), 
                       subplot_kw=Dict("projection"=> proj0), sharey = true)
                       
for (i, var) in enumerate(vars)
    CF[i] = OHC_helper.pcolormesh_ma(ax[i], λ_wrap, ϕ, τy_dict[var], 
                        bounds[i], cm.curl, projPC, true)
    ax[i].set_title(var)
end

[a.coastlines() for a in ax]
[a.set_extent((-180, 180, -70, 56),crs=projPC) for a in ax]
fig.colorbar(CF[1], ax = ax[1:2], fraction = 0.01, pad = 0.03, orientation = "horizontal", label = "[N/m²]")
fig.colorbar(CF[3], ax = ax[3], fraction = 0.01,  pad = 0.03, orientation = "horizontal", label = "[N/m²]")
fig.suptitle("Mean Meridional Wind Stress [1993 - 2017]", y = 0.25)
fig.savefig(plotsdir("native/mean_τy_" * shortnames[vars[1]] * "_" * shortnames[vars[2]] * ".png"), bbox_inches = "tight", dpi = 400);


#plot curl
fig, ax = plt.subplots(1, 3, figsize=(15,20), 
                       subplot_kw=Dict("projection"=> proj0), sharey = true)
                       maximum(curlτ_dict["iter129_bulkformula"]      )
[a.coastlines() for a in ax]
[a.set_extent((120, 285, -70, 56),crs=projPC) for a in ax]
labels = ["Iteration 129", "Iteration 0", "Iteration 129 minus Iteration 0"]

for (i, var) in enumerate(vars)
    CF[i] = OHC_helper.pcolormesh_ma(ax[i], λ_wrap, ϕ, curlτ_dict[var], 
                        1e-12, cm.curl, projPC, true)
    ax[i].set_title(labels[i])
end

fig.colorbar(CF[1], ax = ax[1:2], fraction = 0.01, pad = 0.03, orientation = "horizontal", label = "[N/m³]")
fig.colorbar(CF[3], ax = ax[3], fraction = 0.01,  pad = 0.03, orientation = "horizontal", label = "[N/m³]")
fig.suptitle("Mean Wind Stress Curl [1992 - 2017]", y = 0.31)
fig.savefig(plotsdir("native/mean_curlτ_" * shortnames[vars[1]] * "_" * shortnames[vars[2]] * ".png"), bbox_inches = "tight", dpi = 400);
fig

bounds = [1e-9, 1e-9, 0.5e-9]
fig, ax = plt.subplots(1, 3, figsize=(15,20), 
                       subplot_kw=Dict("projection"=> proj0), sharey = true)
[a.coastlines() for a in ax]
[a.set_extent((120, 285, -70, 56),crs=projPC) for a in ax]
for (i, var) in enumerate(vars)
    CF[i] = OHC_helper.pcolormesh_ma(ax[i], λ_wrap, ϕ, ekup_dict[var], 
    bounds[i], cm.curl, projPC, true)
    ax[i].set_title(var)
end

fig.colorbar(CF[1], ax = ax[1:2], fraction = 0.01, pad = 0.03, orientation = "horizontal", label = "[N/m²]")
fig.colorbar(CF[3], ax = ax[3], fraction = 0.01,  pad = 0.03, orientation = "horizontal", label = "[N/m²]")
fig.suptitle("Mean Ekman upwelling Curl [1992 - 2017]", y = 0.31)
fig.savefig(plotsdir("native/mean_ekup_" * shortnames[vars[1]] * "_" * shortnames[vars[2]] * ".png"), bbox_inches = "tight", dpi = 400);
fig