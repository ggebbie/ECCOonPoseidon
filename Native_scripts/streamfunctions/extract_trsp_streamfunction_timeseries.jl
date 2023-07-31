include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings
using .OHC_helper
import NaNMath as nm
using PyCall

subtract_dict(D, subkey) = Dict(key => value .- D[subkey] for (key, value) in D)
mask_dict(D, lvls) = Dict(key => value[lvls, :] for (key, value) in D)

include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
 
ocean_mask = OHC_helper.wet_pts(Γ)
region = "PAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "full", include_bering = true)
msk = PAC_msk;
tecco = 1992+1/24:1/12:2018

X=collect(-89.0:89.0); Y=reverse(z); #coordinate variables
area_mask = area .* PAC_msk
Xϕ = OHC_helper.ma_zonal_sum(ϕ .* area_mask) ./ OHC_helper.ma_zonal_sum(area_mask)
Xϕ = Xϕ[isfinite.(Xϕ)]
@time for expname in ["iter0_bulkformula", "iter129_bulkformula", "noinitadjust", 
    "seasonalclimatology", "nosfcadjust",
    "seasonalclimatology_iter0", "climatological_tau"]
    println(expname)
    Ψ_exp = OHC_helper.extract_meridionalΨ̄wBolus_timeseries(expname,diagpath, Γ, γ, PAC_msk)
    jldsave(datadir("ΨwBolustimeseries_"*region*"_" * expname *".jld2"); Ψ_exp)
end

fig, ax = plt.subplots()
ax.plot(Xϕ)
fig

LC=LatitudeCircles(round.(Xϕ, digits = 3),Γ)
LC=LatitudeCircles(-89:80,Γ)

LC[80].tabW[1, :, :, :]

Γ.AngleSN[4][121, 1]

unique(LC[80].tabW[:, 1])
unique(LC[80].tabW[:, end])

unique(LC[80].tabS[:, 1])
unique(LC[80].tabS[:, end])
