#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings

import NaNMath as nm
using .OHC_helper
using PyCall

@pyimport seaborn as sns;
@pyimport pandas as pd;
colors =  sns.color_palette("colorblind")[1:4]
labels = ["ECCO", "ECCO Forcing", "ECCO w/out Forcing Adjust.", "CTRL"]
sns.set_theme(context = "poster", style = "ticks", font_scale = 1.2,
              palette = sns.color_palette("colorblind"));
pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;
include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
tecco = 1992+1/24:1/12:2018

runpath,diagpath = listexperiments(exprootdir());
diagpath["seasonalclimatology_iter0"]
diagpath["seasonalclimatology_iter0_multarr0"] = "/batou/ECCOv4r4/exps/seasonalclimatology_iter0/run_run_multarr0/diags/"
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
#load in latitude mask 
ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not", include_bering = false)

#create volume mask
area = readarea(γ)
cell_depths = get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
H = smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);

function volume_average_by_depth(ds::MeshArray, ΔV::MeshArray, γ)
    nz = size(ds, 2)
    vol_avg = zeros(Float32, nz)
    V = zeros(Float32, nz)

    for k=1:nz
        V[k] = Float32(sum(ΔV[:, k]))
    end

    for ff=1:5, k=1:nz
        vol_avg[k] += sum(ds[ff, k] .* ΔV[ff, k]) / V[k]
    end

    return vol_avg
end

#load trend trend_matrices, F is LS estimator
E,F = trend_matrices(tecco)
dθ = Dict()
#need to save these fields for faster plotting! They will not change anymore 
@time for (i, expname) in enumerate(["iter129_bulkformula"])#enumerate(keys(shortnames))
    filelist = searchdir(diagpath[expname], "state_2d_set1") # first filter for state_3d_set1
    datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    β = zeros(50)
    nt = length(datafilelist_S)
    # F =  F
    @time for tt = 1:nt
        println(tt)
        fnameS = datafilelist_S[tt]
        fnameθ = datafilelist_θ[tt]
        sθ = OHC_helper.extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_H)
        sθavg = volume_average_by_depth(sθ, cell_volumes, γ)
        for k in 1:50
            β[k] += F[2,tt] * sθavg[k] 
        end
    end
    dθ[expname] = β
    # ax.scatter(1e2.*β, -z.*1e-3, color = colors[i], label = labels[i])
end

savename = datadir("native/" * region * "_THETA_depth_trends_" * ".jld2")
jldsave(savename, dθ = dθ)