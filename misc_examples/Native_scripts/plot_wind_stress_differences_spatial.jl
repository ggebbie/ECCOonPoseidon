#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LinearAlgebra, LaTeXStrings, PyCall
import PyPlot as plt
import NaNMath as nm
@pyimport matplotlib.animation as anim

cm = pyimport("cmocean.cm");
@pyimport seaborn as sns;

sns.set_theme(context = "poster", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
              
include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "NPAC"; 
ocean_mask = OHC_helper.wet_pts(Γ)
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
sns.set_theme(context = "poster", style = "white", font_scale = 1.0,
              palette = sns.color_palette("deep"));

#computed on standard pressures, 
#may have a different value if using pressures diagnosed from ECCO
#https://ecco-v4-python-tutorial.readthedocs.io/VectorCalculus_ECCO_barotropicVorticity.html
runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

tecco= 1992+1/24:1/12:2018 # ecco years
nt = length(tecco)

reg_λ = LLCcropC(λ,γ); reg_ϕ = LLCcropC(ϕ,γ)

reg_mask = LLCcropC(ocean_mask, γ); reg_mask[reg_mask .== 0] .= NaN
reg_PAC = LLCcropC(PAC_msk, γ); reg_PAC[reg_PAC .== 0.0] .= NaN

τxs = Dict(key => zeros(size(reg_PAC)..., nt) for key in keys(shortnames)); 
τys = Dict(key => zeros(size(τxs[key])...) for key in keys(shortnames)); 

for expname in keys(shortnames)
        filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
        τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    for tt = 1:nt

        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        Tname = τdatafilelist[tt]
        
        @time EXF = γ.read(diagpath[expname]*τdatafilelist[tt],MeshArray(γ,Float32,15))
        τx = EXF[:, 14]; τy = EXF[:, 15]; 
        # curlτ = curl(τx,τy,Γ)
        τE, τN = rotate_uv(τx, τy, Γ)

        τx_reg = LLCcropC(τE,γ); τy_reg = LLCcropC(τN,γ); 
        # curlτ_reg = LLCcropC(curlτ,γ);
        
        τxs[expname][:, :, tt] = τx_reg .* reg_mask
        τys[expname][:, :, tt] = τy_reg .* reg_mask
        # curlτs .+= curlτ_reg ./ nt
    end
end

proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

τ_adj = mean(τxs["iter129_bulkformula"], dims = 3)[:, :, 1]
τ_unadj = mean(τxs["iter0_bulkformula"], dims = 3)[:, :, 1]
Δτ = τ_adj .- τ_unadj
rotate_velocity!
fig, ax = plt.subplots(ncols = 2, figsize=(15,20), 
subplot_kw=Dict("projection"=> proj0), sharey = true)

bounds = maximum([nm.maximum(abs.(var .* reg_PAC)) for var in [τ_adj, τ_unadj]])
var_labels = ["Adjusted τ", "Unadjusted τ"]
levels = collect(-bounds:0.01:bounds)
cfs = []
for (i, var) in enumerate([τ_adj, τ_unadj])
    cf = ax[i].contourf(reg_λ, reg_ϕ,  var, vmin = -bounds, vmax = bounds, levels = levels,
    transform=projPC, cmap = cm.curl, extend = "both")
    push!(cfs, cf)
    gl = ax[i].gridlines(crs=projPC, draw_labels=true, linewidth=2, 
                    color="gray", alpha=0, linestyle="--")
    gl.top_labels = false; gl.right_labels = false
    
    ax[i].coastlines()
    ax[i].set_title(var_labels[i])
    ax[i].set_extent((-326, -60, -60, 70),crs=projPC)
end
fig.colorbar(cfs[1], ax=ax, orientation = "horizontal", 
fraction=0.01, pad = 0.05, extend = "both", label = L"Sv")
fig

fig, ax = plt.subplots(figsize=(15,20), 
subplot_kw=Dict("projection"=> proj0), sharey = true)
bounds = round(nm.maximum(abs.(Δτ .* reg_PAC)), digits = 2)
levels = collect(-bounds:0.005:bounds)
cf = ax.contourf(reg_λ, reg_ϕ,  Δτ, vmin = -bounds, vmax = bounds, 
    levels = levels,
    transform=projPC, cmap = cm.curl, extend = "both")

gl = ax.gridlines(crs=projPC, draw_labels=true, linewidth=2, 
                color="gray", alpha=0, linestyle="--")
gl.top_labels = false; gl.right_labels = false

ax.coastlines()
ax.set_title("Adjustment to τ")
ax.set_extent((-326, -60, -60, 70),crs=projPC)

fig.colorbar(cf, ax=ax, orientation = "horizontal", 
fraction=0.03, pad = 0.05, extend = "both", label = L"Sv")
fig