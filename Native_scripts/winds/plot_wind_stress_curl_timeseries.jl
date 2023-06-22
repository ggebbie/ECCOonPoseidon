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

[curlτ_dict[var] = zeros(Float32, nt) for var in vars]
[ekup_dict[var] = zeros(Float32, nt) for var in vars]

#plot ekman pumping velocity
f(ϕ) = abs(ϕ) > 1 ? (2 * (2π / 86400) * sind(ϕ)) : Inf32
fs = MeshArray(γ,Float32); fs.f .= map.(x -> f.(x), ϕ.f)
finv = 1 ./ fs; finv = Float32.(finv) 
σ2 = jldopen(datadir("native/native_sigma2_timeavg_1992_2017.jld2"))["ρθSavg"]
σ2 = Dict(k => vertical_average(v["σ2"] .+ 1000) for (k, v) in σ2)
inv_σ2 = Dict()
function inv_ma(ma::MeshArray)
    ma_copy = deepcopy(ma)
    for ff in eachindex(ma)
        ma_copy.f[ff][ma.f[ff] .== 0.0] .= Inf
        ma_copy.f[ff] .= 1 ./ ma_copy.f[ff]
    end
    return ma_copy
end

#time mean 
for (k, v) in σ2
    inv_σ2[k] = inv_ma(σ2[k]);
end

area_mask = area .* PAC_msk; 
area_weight(x, w) = sum(x .* w)
area_average(x, w) = sum(x .* w) / (sum(w))
p₀ = 2000; nt = 312 
for expname in vars
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    σ = MeshArray(γ, Float32, 50)
    for tt = 1:nt
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)

        θname = datafilelist_θ[tt]
        Tname = τdatafilelist[tt]
        #     # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        # @time θSz = γ.read(diagpath[expname]*θname,MeshArray(γ,Float32,100))
        # θz = θSz[:, 1:50]; Sz = θSz[:, 51:end]
        # for ijk in eachindex(θz)
        #     σ.f[ijk] .= OHC_helper.densityJMD95.(θz.f[ijk],Sz.f[ijk], 0.0, p₀) #EOS from MITGCM 
        # end 
        τx, τy = OHC_helper.extract_ocnTAU(diagpath, expname , τdatafilelist, tt, γ)
        curlτ_dict[expname][tt] = area_average(MeshArrays.curl(τx, τy, Γ), area_mask)
        # inv_σ = inv_ma(vertical_average(σ))
        # println(maximum(inv_σ))
        # fσinv = inv_σ .* finv
        # ek = MeshArrays.curl(τx .* fσinv, τy .* fσinv, Γ)
        # ekup_dict[expname][tt] =  area_weight(ek, Γ.DXG .* PAC_msk)
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


