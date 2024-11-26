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

region = "NPAC"; 
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
curlτ_dict = Dict();
svtrsp_dict = Dict();
wek_dict = Dict();

vars = ["iter129_bulkformula", "iter0_bulkformula", "seasonalclimatology"]

[wek_dict[var] = zeros(Float32, nt) for var in vars]
[curlτ_dict[var] = zeros(Float32, nt) for var in vars]
[svtrsp_dict[var] = zeros(Float32, nt) for var in vars]


area_mask = area .* PAC_msk; 
area_weight(x, w) = sum(x .* w)
area_average(x, w) = sum(x .* w) / (sum(w))
cs, sn = OHC_helper.get_cs_and_sn(γ)

ω = 2π / 86400; r = 6.387e6
beta = 2 * ω * cosd(45) / r ; rho = 1000
ϕ_min_mask, ϕ_max_mask = OHC_helper.get_ϕ_max_min_mask(region, Γ, λ, ϕ, basins, basin_list)
dE, dN =  OHC_helper.rotate_UV_native(Float32.(Γ.DXG), Float32.(Γ.DYG), cs, sn);
dE = dE .* ϕ_min_mask; dE = abs.(dE) 


#plot ekman pumping velocity
f0(ϕ) = abs(ϕ) > 1 ? (2 * ω * sind(ϕ)) : Inf32
f1(ϕ) =  (2 * ω * sind(ϕ))

ϕ = Float32.(ϕ)
ϕU, ϕV = OHC_helper.tofaces(ϕ, ϕ, Γ)

finv = MeshArray(γ,Float32); finv.f .= map.(x -> inv.(f0.(x)), ϕ.f)
fUinv = MeshArray(γ,Float32); fUinv.f .= map.(x -> inv.(f0.(x)), ϕU.f)
fVinv = MeshArray(γ,Float32); fVinv.f .= map.(x -> inv.(f0.(x)), ϕV.f)

#plot ekman pumping velocity
β0(ϕ) = abs(ϕ) < 80 ? (2 * ω * cosd(ϕ) /  6.387e6) : Inf32
invβ = MeshArray(γ,Float32); invβ.f .= map.(x -> inv.(β0.(x)), ϕ.f)
f = MeshArray(γ,Float32); f.f .= map.(x -> f1.(x), ϕ.f)

X = sum(ϕ_min_mask) * 100e3
area_lat_mask = area_mask .* ϕ_min_mask

for expname in vars
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    sv = MeshArray(γ,Float32); 

    for tt = 1:312
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)

        τx, τy = OHC_helper.extract_ocnTAU(diagpath, expname , τdatafilelist, tt, γ)
        curlτ = MeshArrays.curl(τx, τy, Γ)
        curlτ_dict[expname][tt] = area_average(Float32.(curlτ), area_mask)

        wek =  MeshArrays.curl(τx .* fUinv, τy .* fVinv, Γ) 
        for ff = 1:5
            sv.f[ff] .= wek.f[ff] .* invβ.f[ff] .* f.f[ff]
            wek.f[ff] .= wek.f[ff] .* area_mask.f[ff]

        end
        svtrsp_dict[expname][tt] = area_average(Float32.(sv), area_mask) .* X
        wek_dict[expname][tt] = sum(wek)

    end
end
mean(curlτ_dict["iter129_bulkformula"])
savename = datadir("native/" * region * "_areamean_τcurl.jld2")
jldsave(savename, curlτ_dict = curlτ_dict)

savename = datadir("native/" * region * "_zonalint_svtrsp.jld2")
jldsave(savename, svtrsp_dict = svtrsp_dict)


fig, ax, = plt.subplots()
ax.plot(svtrsp_dict["iter129_bulkformula"] )
fig

fig, ax, = plt.subplots()
ax.plot(curlτ_dict["iter129_bulkformula"] )
fig

fig, ax, = plt.subplots()
ax.plot(1e-6 .* svtrsp_dict["iter129_bulkformula"] )
ax.plot(1e-6 .* svtrsp_dict["seasonalclimatology"] )
fig

fig, ax, = plt.subplots()
ax.plot(1e-6 .* (wek_dict["seasonalclimatology"]))
mean(1e-6 .* wek_dict["seasonalclimatology"])

fig