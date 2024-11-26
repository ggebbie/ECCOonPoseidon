include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
area32 = Float32.(area)

ocean_mask = wet_pts(Γ)
region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
H = vertical_sum(cell_depths); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H
cell_volumes = get_cell_volumes(area, cell_depths);

runpath,diagpath = listexperiments(exprootdir());

include(srcdir("plot_and_dir_config.jl"))

ϕ_min_mask, ϕ_max_mask = get_ϕ_max_min_mask(region, Γ, λ, ϕ, basins, basin_list)
mskC, mskW, mskS = get_msk(Γ)

#define control volume 
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= -z[:].<= uplvl)
ctrl_vol = sum(cell_volumes[:, lvls])

function interpolate_to_vertical_faces!(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    ma_interp::MeshArrays.gcmarray{T, 2, Matrix{T}})  where T<:Real
    for k = 1:49 
        ma_interp[:, k+1] = (ma[:, k] .+ ma[:, k+1]) ./2 #linear interpolate to faces
    end
end


expname = "iter129_bulkformula"

filelist = searchdir(diagpath[expname], "state_2d_set1") # first filter for state_3d_set1
datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"trsp_3d_set2")
tt = 1
fnameS = datafilelist_S[tt]
fnameθ = datafilelist_θ[tt]
sθ = extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_H)


sθinterp= MeshArray(γ,Float32,50); 
sθinterp2= MeshArray(γ,Float32,50); 

nt = 312;


interpolate_to_vertical_faces!(sθ, sθinterp)
function interpolate_to_vertical_faces!(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    ma_interp::MeshArrays.gcmarray{T, 2, Matrix{T}})  where T<:Real
    for k = 1:49 
        ma_interp[:, k+1] = (ma[:, k] .+ ma[:, k+1]) ./2 #linear interpolate to faces
    end
end

test = zeros(50, 270, 90)
test_interp = zeros(50, 270, 90)
test_interp2 = zeros(50, 270, 90)

[test[i, :, :] .= view(sθinterp.f[4, i], :, :) for i in 1:50]
[test_interp[i, :, :] .= view(sθ.f[4, i], :, :) for i in 1:50]

x = test[:, 100, 2]; x[x .== 0.0] .= NaN 
y = test_interp[:, 100, 2]; y[y .== 0.0] .= NaN 

zF = cumsum(vcat(0, Γ.DRF))[2:end-1]
mult(x) = Float32.((x .- z[1:end-1]) ./ (z[2:end] .- z[1:end-1]))
mult(zF)[lvls]

function interpolate_to_vertical_faces!(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    ma_interp::MeshArrays.gcmarray{T, 2, Matrix{T}})  where T<:Real
    for k = 1:48
        ma_interp[:, k+1] = (ma[:, k] .* (1 - mult(zF)[k+1])) .+ (ma[:, k+1].* mult(zF)[k+1]) #linear interpolate to faces
    end
end
interpolate_to_vertical_faces!(sθ, sθinterp2)
[test_interp2[i, :, :] .= view(sθinterp2.f[4, i], :, :) for i in 1:50]
zs = test_interp2[:, 100, 2]; y[y .== 0.0] .= NaN 

fig, ax = plt.subplots()
ax.scatter(x, z, label = "x")
ax.scatter(y[1:49], zF, label = "y")
# ax.scatter(zs[1:49], zF, label = "z")

ax.set_ylim(1000, 3000)
ax.set_xlim(0, 3)

ax.invert_yaxis()
ax.legend()
fig