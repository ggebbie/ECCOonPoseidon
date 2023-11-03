#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../../../src/intro.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LaTeXStrings, PyCall, GibbsSeaWater, Dierckx, Interpolations


import NaNMath as nm
@pyimport seaborn as sns
@pyimport pandas as pd
@pyimport cmocean.cm as cmo

include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());


include(srcdir("plot_and_dir_config.jl"))

region = "PAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "full")
cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths);

fname = "/batou/eccodrive/files/Version4/Release4/input_init/xx_kapgm.0000000129.data";
# fname = "/batou/eccodrive/files/Version4/Release4/input_init/total_kapgm_r009bit11.bin";

# fname = "/batou/eccodrive/files/Version4/Release4/input_init/xx_kapredi.0000000129.data";


κGM = read_bin(fname,Float32,γ)
κGM_zonal = zonal_average(κGM, cell_volumes)
vmax = nm.maximum(κGM_zonal) 
ϕ_avg = zonal_average(ϕ, PAC_msk .* area); wherenan = (!isnan).(ϕ_avg .* 1)

fig,axs=plt.subplots(figsize = (15, 7.5))
cb = axs.pcolormesh(ϕ_avg[wherenan], z,  κGM_zonal[:, wherenan], cmap=cmo.delta, 
vmin = -vmax, vmax = vmax)
fig.colorbar(cb)
axs.invert_yaxis()
fig