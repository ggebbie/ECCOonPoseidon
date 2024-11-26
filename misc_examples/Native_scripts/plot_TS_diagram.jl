#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour, MeshArrays, 
MITgcmTools, JLD2, DrWatson, Statistics, JLD2, 
Printf, PyCall, LaTeXStrings
import NaNMath as nm, PyPlot as plt
@pyimport cmocean.cm as cmo
@pyimport seaborn as sns;
sns.set_theme(context = "poster", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
include(srcdir("config_exp.jl"))'

(ϕ,λ) = latlonC(γ);
area = readarea(γ);
meshgrid(x, y) = (x' .* ones(length(y)), ones(length(x))' .* y);

ocean_mask = OHC_helper.wet_pts(Γ);
cell_depths = OHC_helper.get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);

runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
tecco= 1992+1/24:1/12:2018 # ecco years
nt = length(tecco)
#read in the first time step of S and θ
expname = "iter129_bulkformula"

filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

#select time and location to do the transect 
tt = 1
latlons = ((-90, -20), (-120, -10), (-179, -30), (-120, 30), (-150, -30))
llcolors = ["red", "blue", "green", "purple", "orange"]
fig, ax = plt.subplots(1,1, figsize = (10, 10));

transects = Dict(); p₀ = 2000
for (ipair, llpair) in enumerate(latlons)
    face, index, _ = OHC_helper.findlatlon(λ, ϕ, llpair[1], llpair[2]);
    θS = γ.read(diagpath[expname]*datafilelist_θ[tt],MeshArray(γ,Float32,100));

    θ_list = Float32[]; S_list = Float32[]; σ_list = Float32[]; vol_list = Float64[];
    θ = θS[:, 1:50]; S = θS[:, 51:end]; 

    for k = 20:50 #k = 20:43 means from z ∈ (-3000, -300)
        σ2 = OHC_helper.densityJMD95(θ.f[face, k][index],S.f[face, k][index], NaN, p₀) .- 1000. #EOS from MITGCM 
        push!(θ_list, θ.f[face, k][index])
        push!(S_list, S.f[face, k][index])
        push!(σ_list, σ2)
        push!(vol_list, cell_volumes[face, k][index])
    end
    [(var_list[vol_list .== 0.0] .= NaN) for var_list in [θ_list, S_list, σ_list]];
    transects[llpair] = Dict{String}{Tuple{Float32, Float32}}()
    transects[llpair]["θ extrema"] = nm.extrema(θ_list);
    transects[llpair]["S extrema"] = nm.extrema(S_list)
    
    ax.scatter(S_list, θ_list, c = llcolors[ipair], alpha = 0.3, label = "(lon, lat): " *  string(latlons[ipair]));
end
S_min = minimum([value["S extrema"][1] for value in values(transects)])
S_max = maximum([value["S extrema"][2] for value in values(transects)])
θ_min = minimum([value["θ extrema"][1] for value in values(transects)])
θ_max = maximum([value["θ extrema"][2] for value in values(transects)])

S_levels = S_min:0.005:S_max; θ_levels = θ_min:0.1:θ_max;
S_grid, θ_grid = meshgrid(S_levels, θ_levels);
σgrid = OHC_helper.densityJMD95.(θ_grid,S_grid, 0.0, p₀) .- 1000.; #EOS from MITGCM does not require ref pressure

sns.set_theme(context = "poster", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
CS = ax.contour(S_grid, θ_grid, σgrid, levels = sigma2grid(), colors = "black", linewidths = 0.75);
ax.clabel(CS, fontsize=15, inline=true, fmt = "%.2f");
ax.set_xlabel("Practical Salinity"); ax.set_ylabel("Potential Temperature");
ax.set_title("Pacific Ocean T-S Diagram (t = 1992)"); ax.legend()
ax.set_ylim(1.0, 10)
fig.savefig(plotsdir("native/TSDiagram_t0.png"), bbox_inches = "tight")

proj0 = ECCOonPoseidon.cartopy.crs.Robinson(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

fig2, axs = plt.subplots(1, 1, figsize=(10,10), subplot_kw=Dict("projection"=> proj0))
for (i, pairs) in enumerate(latlons)
    axs.scatter(x=pairs[1], y=pairs[2], c = llcolors[i], s = 20, transform=projPC, label = "(lon, lat): " *  string(latlons[i]))
end
axs.legend()
axs.set_extent((-190, 180, -70, 56),crs=projPC)
gl = axs.gridlines(crs=projPC, draw_labels=true,
                linewidth=2, color="gray", alpha=0.2, linestyle="--")
gl.top_labels = false; gl.bottom_labels = true
gl.right_labels = false
axs.coastlines(resolution="110m")
fig2
fig2.savefig(plotsdir("native/TSDiagramLocations_t0.png"), bbox_inches = "tight")
