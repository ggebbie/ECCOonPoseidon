# this script will compile and save a basinwide average of OHC 
# for a particular basin (basin_name) 

include("../src/intro.jl")
include("../src/OHC_helper.jl")
using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools,
PyPlot, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf
using .OHC_helper

pygui(false)
standardize(x) = (x .- mean(x))/std(x)

include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir())
shortnames = expnames()
marks = expsymbols()
nexp = length(shortnames) # number of experiments

nz = 50
ocean_mask = Γ.hFacC[:,1] 
ocean_mask[findall(ocean_mask.>0.0)].=1.0
ocean_mask[findall(ocean_mask.<=0.0)].=0.0
uplvl = -2000; botlvl = -3000
lvls = findall( botlvl .<= z[:].<= uplvl)

""" comparing different surface temperatures """
cell_depths = get_cell_depths(ocean_mask, Δz, Γ.hFacC)
cell_volumes = get_cell_volumes(area, cell_depths)
area_vec = vec(area .* ocean_mask)
vol_sfcvec = vec(cell_volumes[:, 1])

θ_sfc = get_mnth1_θ(shortnames, diagpath, γ)
keys_ = keys(θ_sfc)
ignore_list = ["noIA", "129ff"]
new_keys = setdiff(keys_, ignore_list)

fig, ax = plt.subplots(4)
fig.suptitle("ECCO SST distribution for Jan. 15, 1992")
i = 0 
lvl = 1 
plot_ocn_field_histogram!(θ_sfc, area_vec, new_keys, lvl, ax)
fig.savefig(plotsdir("SST_hist"))

""" comparing different internal global temperatures 2km -3km """
vol_sfcvec = vec(cell_volumes[:, lvls])
θ_sfc = get_mnth1_θ(shortnames, diagpath, γ)

fig, ax = plt.subplots(4)
fig.suptitle("ECCO T distribution for Jan. 15, 1992")
plot_ocn_field_histogram!(θ_sfc, vol_sfcvec, new_keys, lvls, ax)
fig.savefig(plotsdir("T_hist"))

"""comparing different internal global temperatures 2km -3km 
 in the Pacific """
basin_name = "Pacific"
basinID=findall(basin_list.==basin_name)[1]
basin_mask=similar(basins)
for ff in 1:length(area)
    above_SO = (ϕ[ff] .> -56.0) #removes southern ocean 
    basin_mask[ff] .= ocean_mask[ff].*(basins[ff].==basinID) .* above_SO
end
cell_depths = get_cell_depths(basin_mask, Δz, Γ.hFacC)
cell_volumes = get_cell_volumes(area, cell_depths)
vol_sfcvec = vec(cell_volumes[:, lvls])
fig, ax = plt.subplots(4)
fig.suptitle("ECCO T distribution for Jan. 15, 1992")
plot_ocn_field_histogram!(θ_sfc, vol_sfcvec, new_keys, lvls, ax)
fig.savefig(plotsdir("T_hist_pac"))

""" comparing total diffusivities (in pacific) at 2km-3km"""
basin_name = "Pacific"
basinID=findall(basin_list.==basin_name)[1]
basin_mask=similar(basins)
for ff in 1:length(area)
    above_SO = (ϕ[ff] .> -56.0) #removes southern ocean 
    basin_mask[ff] .= ocean_mask[ff].*(basins[ff].==basinID) .* above_SO
end
cell_depths = get_cell_depths(basin_mask, Δz, Γ.hFacC)
cell_volumes = get_cell_volumes(area, cell_depths)
vol_sfcvec = vec(cell_volumes[:, lvls])
fig, ax = plt.subplots(3, figsize = (7,6))
fig.suptitle("ECCO Diffusion coefs. for Jan. 15, 1992")
diffs, diffs_keys = get_diff_arrs(γ)
plot_ocn_init_histogram!(diffs, vol_sfcvec, diffs_keys, lvls, ax)
for a in ax
    a.set_ylabel("Frequency")
    a.set_xlabel(" m^2/s")
end
fig.tight_layout()
fig.savefig(plotsdir("diff_hist"),bbox_inches="tight")
""" comparing surface energy globally"""
sfc_adjust = get_sfcfrc_fnames("adjust")
