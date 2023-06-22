#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LinearAlgebra, LaTeXStrings, PyCall
import PyPlot as plt, NaNMath as nm
@pyimport matplotlib.animation as anim

cm = pyimport("cmocean.cm");
@pyimport seaborn as sns
sns.set_theme(context = "paper", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "PAC56"; 
ocean_mask = OHC_helper.wet_pts(Γ)
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")

runpath,diagpath = listexperiments(exprootdir());

#read in the first time step of S and θ
expname = "seasonalclimatology"

uplvl = -2.4e3; botlvl = -2.6e3;
lvls = findall( botlvl .<= z[:].<= uplvl)

filelist = searchdir(diagpath[expname],"state_3d_set1") 
datafilelist_θ  = filter(x -> occursin("data",x),filelist)

filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

θ = Float32[]
τE = Float32[]
tecco= 1992+1/24:1/12:2018; # ecco years
nt = length(tecco); nz = length(z);

fcycle = 1 # units: yr^{-1}
overtones= 4; # semi-annual, quad-annual, etc.: helps capture an asymmetric seasonal cycle
Ecycle,Fcycle = seasonal_matrices(fcycle,tecco,overtones)
tt = 1
diagpath["nointerannual"]*τdatafilelist[tt]
moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]
standardize(x) = (x .- mean(x)) ./ std(x)

latlons = ((-90, -20), (-175, -30), (-179, +30), (-145, 30), (157, 30))
llcolors = ["red", "blue", "green", "purple", "orange"]
θ_dict = Dict(); τ_dict = Dict()
fidict = Dict()
for (ipair, llpair) in enumerate(latlons)
    face, index, _ = OHC_helper.findlatlon(λ, ϕ, llpair[1], llpair[2]);
    fidict[llpair] = (face, index)
    θ_dict[llpair] = []
    τ_dict[llpair] = []
end
for tt = 1:nt
    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    @time θσ1 = γ.read(diagpath[expname]*datafilelist_θ[tt],MeshArray(γ,Float32,nz))
    #if bulkformula τx and τy are at 14 and 15. 
    #if flux-forced, τx and τy are at 9 and 10. 

    @time EXF = γ.read(diagpath[expname]*τdatafilelist[tt],MeshArray(γ,Float32,10))

    τx = EXF[:, 9]; τy = EXF[:, 10]; curlτ = curl(τx,τy,Γ)
    τxC, τyC = velocity2centers(τx, τy, Γ)
    τE_, τN_ = rotate_uv(τxC, τyC, Γ)
    for (ipair, llpair) in enumerate(latlons)
        face, index = fidict[llpair]
        push!(θ_dict[llpair], θσ1[face, lvls[1]][index])
        push!(τ_dict[llpair], τE_[face][index])
    end
end

fig, axes = plt.subplots(nrows = 2, ncols = 5, figsize=(15,3), 
                        sharex = true, sharey = true)

for (ipair, llpair) in enumerate(latlons)
    axes1 = axes[1, :]
    axes2 = axes[2, :]
    θ  = θ_dict[llpair]
    τE = -τ_dict[llpair]

    τE_rolling = moving_average(τE,13)
    axes1[ipair].set_title("(Lon, Lat) = " * string(latlons[ipair]) * "\n Temperature")
    axes1[ipair].plot(tecco, standardize(θ), color = llcolors[ipair], linewidth = 1.0)
    axes2[ipair].set_title("negative Zonal Wind Stress")
    axes2[ipair].plot(tecco[7:end-6], standardize(τE_rolling),  
                     color = llcolors[ipair], alpha = 0.7, 
                     linewidth = 1.0, zorder = 2.5)
                     
end
[ax.set_ylim(-3, 3) for ax in axes]
fig
fig.savefig(plotsdir("native/THETA_TAU_Coherence_" * expname *".png"), bbox_inches = "tight")

proj0 = ECCOonPoseidon.cartopy.crs.Robinson(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
fig2, axs = plt.subplots(1, 1, figsize=(5,5), subplot_kw=Dict("projection"=> proj0))
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
sns.move_legend(axs, "lower center", 
                bbox_to_anchor=(0.5, -0.5), ncol=3, 
                frameon=true, borderaxespad=0.)

fig2
fig2.savefig(plotsdir("native/TSDiagramLocations_t0.png"), bbox_inches = "tight")


# function spatial_correlation(bs, Us)
#     nt = length(bs)
#     Up = Us .- mean(Us, dims = 3)
#     bp = bs .- mean(bs)
#     σU = std(Up, dims = 3)[:, :, 1]; σb = std(bp)
#     σUσb = σU .* σb
    
#     c = zeros(size(Us[:, :, 1])...)
#     [c[:, :] .+= (Up[:, :, tt] .* bp[tt] ./ nt) ./ σUσb for tt in 1:nt]
#     return c
# end


# proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
# projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

# fig, axes = plt.subplots(ncols= 3, figsize=(15,20), 
# subplot_kw=Dict("projection"=> proj0), sharey = true)
# vars = [τxs, τys, curlτs ]; var_labels = ["τ_x", "τ_y", "curl τ"]
# bounds = maximum([nm.maximum(abs.(var)) for var in vars])
# cfs = Any[]

# dtbs =  (bs[2:end] .- bs[1:end-1]) ./ 2.628e+6 #time derivative

# for (i, ax) in enumerate(axes)
#     # linear interpolate 
#     var = (vars[i][:, :, 1:end-1] .+ vars[i][:, :, 2:end]) ./ 2
#     cs = spatial_correlation(dtbs, var)
#     cf = ax.pcolormesh(reg_λ, reg_ϕ,  cs, 
#                    vmin = -bounds, vmax = bounds, shading="nearest", 
#                    transform=projPC, rasterized = true, cmap = cm.curl)
#     push!(cfs, cf)
#     ax.scatter(reg_λ[ilonilat], reg_ϕ[ilonilat], s=100, marker = "*",c = "blue")
#     gl = ax.gridlines(crs=projPC, draw_labels=true, linewidth=2, 
#                   color="gray", alpha=0, linestyle="--")
#     gl.top_labels = false; gl.right_labels = false
#     if i > 1
#         gl.left_labels = false
#     end
#     # ax.coastlines()
#     ax.set_title(var_labels[i])
#     ax.set_extent((-326, -60, -60, 70),crs=projPC)

# end
# # fig.suptitle("Correlation of dθ/dt with atmospheric variables")
# # fig.tight_layout()
# fig.colorbar(cfs[1], ax=axes, orientation = "horizontal", 
# fraction=0.01, pad = 0.05, extend = "both", label = L"\rho")
# fig