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

sns.set_theme(context = "paper", style = "ticks",
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
curlτ_dict = Dict(); ekup_dict = Dict()
vars = ["iter129_bulkformula", "seasonalclimatology"]


area_weight(x, w) = sum(x .* w)
area_average(x, w) = sum(x .* w) / (sum(w))
p₀ = 2000; nt = 312 
for expname in vars
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    for tt = 1:nt
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)

        θname = datafilelist_θ[tt]
        Tname = τdatafilelist[tt]

        τx, τy = OHC_helper.extract_ocnTAU(diagpath, expname , τdatafilelist, tt, γ)
        curlτ_dict[expname][tt] = area_average(MeshArrays.curl(τx, τy, Γ), area_mask)

    end
end

diff_dict(d, x, y) = d[x] .- d[y]
diff_label = shortnames[vars[1]] * " minus " * shortnames[vars[2]]
curlτ_dict[diff_label] = diff_dict(curlτ_dict, vars[1], vars[2])
# ekup_dict[diff_label] = diff_dict(ekup_dict, vars[1], vars[2])
push!(vars, diff_label) 

#setting up plotting stuff       
bounds = [0.15, 0.15, 0.04]
CF = Any[1, 1, 1]


#plot curlτ
fig, ax = plt.subplots(1, 3, figsize=(15,3), sharey = false)
fig.suptitle("Area Averaged curl(τ) [1993 - 2017]")
for (i, var) in enumerate(vars)
    ax[i].plot(tecco[1:nt], curlτ_dict[var][1:nt]) #m/s -> cm/century 
    ax[i].set_title(var)
end
ax[1].set_ylabel("cm per century")
fig

# #plot ekman upwelling
# fig, ax = plt.subplots(1, 3, figsize=(15,5), sharey = true)
# fig.suptitle("Net Ekman Upwelling Transport ρ⁻¹curl(τ/f) [1993 - 2017]")

# for (i, var) in enumerate(vars[1:3])
#     ax[i].plot(tecco[1:nt], ekup_dict[var][1:nt]) #m/s -> cm/century 
#     ax[i].set_title(var)
# end
# ax[1].set_ylabel("cm per century")
# fig


