#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour, MeshArrays, 
MITgcmTools, JLD2, DrWatson, Statistics, JLD2, 
Printf, PyCall, LaTeXStrings
import NaNMath as nm, PyPlot as plt
# @pyimport cmocean.cm as cmo
@pyimport seaborn as sns;
sns.set_theme(context = "poster", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
include(srcdir("config_exp.jl"))'

(ϕ,λ) = latlonC(γ);
area = readarea(γ);
meshgrid(x, y) = (x' .* ones(length(y)), ones(length(x))' .* y);

region = "NPAC"
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);

runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
tecco= 1992+1/24:1/12:2018 # ecco years
nt = length(tecco); nz = length(z)

expname = "iter129_bulkformula"
filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"

λ_wrap = deepcopy(λ); 
for ff ∈ [3,4] 
    λ_wrap[ff][λ_wrap[ff] .<= 0 ] = λ_wrap[ff][λ_wrap[ff] .<= 0 ] .+ 360
end

U=0.0*Γ.hFacW; fill!(U, 0.0) #varmeta not quite correct
V=0.0*Γ.hFacS; fill!(V, 0.0)#varmeta not quite correct

for tt in 1:nt
    fname = datafilelist[tt]
    u, v, w = OHC_helper.extract_velocities(diagpath, expname , datafilelist_uvw[tt], γ)
    println(tt)
    for i in eachindex(U)
        U[i]=U[i]+(u[i]/nt)
        V[i]=V[i]+(v[i]/nt)
    end
end

(Utr,Vtr)=UVtoTransport(U,V,Γ)

cs, sn = OHC_helper.get_cs_and_sn(γ)
Evel, Nvel = OHC_helper.rotate_UV_native(Float32.(Utr), Float32.(Vtr), cs, sn);
ϕ_min_mask, ϕ_max_mask = OHC_helper.get_ϕ_max_min_mask(region, Γ, λ, ϕ, basins, basin_list)

latitudecircle(x, ϕ_min_mask) = vcat([x.f[ff][ϕ_min_mask.f[ff] .!= 0.0] for ff = 1:5]...)

lons = latitudecircle(λ, ϕ_min_mask)
lons[lons .< 0.0] .+= 360
Ntrsps = zeros(50, length(lons))
[Ntrsps[k, :] .= latitudecircle(Nvel[:, k], ϕ_min_mask) for k = 1:50]
zonal_trsp = sum(Ntrsps, dims = 2)
zonal_trsp_noWB = sum(Ntrsps[:, 6:end], dims = 2)

Ntrsps[Ntrsps .== 0.0] .= NaN
fig, ax = plt.subplots()
levels = -0.1:0.02:0.1
CS = ax.contourf(lons, z, 1e-6 .* Ntrsps, cmap = cmo.balance, 
vmin = -0.1, vmax = 0.1, levels = levels, extend = "both");
ax.set_xlabel("Longitude"); ax.set_ylabel("Depth");
ax.set_title("Northward Velocity"); ax.set_xlim(120, 270)
fig

fig, ax = plt.subplots()
CS = ax.plot(1e-6 .* zonal_trsp[25:end], z[25:end])
CS = ax.plot(1e-6 .* zonal_trsp_noWB[25:end], z[25:end])

# ax.set_ylim(1.0, 10)
fig


zonal_trsp
# fig.savefig(plotsdir("native/TSDiagram_t0.png"), bbox_inches = "tight")