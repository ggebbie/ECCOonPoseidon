#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
using .OHC_helper
import NaNMath as nm
include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");

(ϕ,λ) = latlonC(γ)

function tofaces(uFLD::MeshArrays.gcmarray{T, 1, Matrix{T}},
    vFLD::MeshArrays.gcmarray{T, 1, Matrix{T}}) where T<:Real
    Ugrid = T.(similar(uFLD))
    Vgrid = T.(similar(uFLD))
    (tmpU,tmpV)=exch_UV(uFLD,vFLD)

    for a in 1:5
        (s1,s2)=size(uFLD.f[a])
        tmpU1=view(tmpU.f[a],1:s1,1:s2)
        tmpU2=view(tmpU.f[a],2:s1+1,1:s2)
        tmpV1=view(tmpV.f[a],1:s1,1:s2)
        tmpV2=view(tmpV.f[a],1:s1,2:s2+1)
        Ugrid.f[a] =(tmpU1+tmpU2) / 2
        Vgrid.f[a] =(tmpV1+tmpV2) / 2
    end

return Ugrid, Vgrid
end

function UVtoTrsp(U::MeshArray,V::MeshArray,G::NamedTuple)
    uTr=deepcopy(U)
    vTr=deepcopy(V)
    for i in eachindex(U)
        uTr.f[i] .=G.DRF[i[2]]*U.f[i].*G.DYG.f[i[1]]
        vTr.f[i] .=G.DRF[i[2]]*V.f[i].*G.DXG.f[i[1]]
    end
    return uTr,vTr
end

  
@time xFU, xFV = tofaces(Γ.XC, Γ.XC)
@time yFU, yFV = tofaces(Γ.YC, Γ.YC)

expname = "iter129_bulkformula"
runpath,diagpath = listexperiments(exprootdir());

filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"trsp_3d_set2")
datafilelist_H  = filter(x -> occursin("data",x),filelist) 

tt = 100

uθ_conv3D = MeshArray(γ,Float32,50);
dθλ = γ.read(diagpath[expname]*datafilelist_H[tt],MeshArray(γ,Float32,200))
κθx = dθλ[ :, 1:50]; κθy = dθλ[:, 51:100]
uθx = dθλ[:, 101:150]; uθy = dθλ[:, 151:200]; dθλ = nothing 
OHC_helper.calc_UV_conv3D!(uθx, uθy, uθ_conv3D);

θ = γ.read(diagpath[expname]*datafilelist_θ[tt],MeshArray(γ,Float32, 50))
UV = γ.read(diagpath[expname]*datafilelist_uvw[tt],MeshArray(γ,Float32,100))
u = UV[:, 1:50]; v = UV[:, 51:100]
U, V = UVtoTrsp(u, v, Γ)

θFU, θFV = MeshArray(γ,Float32,50), MeshArray(γ,Float32,50)
θU, θV = MeshArray(γ,Float32,50), MeshArray(γ,Float32,50)
uθ_conv3D_approx = MeshArray(γ,Float32,50);

for k = 1:50 
    tmp = θ[:, k]
    tmp1, tmp2 = tofaces(tmp, tmp)
    θFU[:, k] = deepcopy(tmp1)
    θFV[:, k] = deepcopy(tmp2)
    θU[:, k] = U[:, k] .* tmp1
    θV[:, k] = V[:, k] .* tmp2
end
OHC_helper.calc_UV_conv3D!(θU, θV, uθ_conv3D_approx);


import PyPlot as plt
fig, ax = plt.subplots(figsize = (7, 7))
ff = 1
ax.scatter(vec(Γ.XC[ff][40:60,150:180]), vec(Γ.YC[ff][40:60, 150:180]), c = "k", label = "center")
ax.scatter(vec(xFU[ff][40:60,150:180]), vec(yFU[ff][40:60, 150:180]), c = "r", label = "U Faces")
ax.scatter(vec(xFV[ff][40:60, 150:180]), vec(yFV[ff][40:60, 150:180]), c = "b", label = "V Faces")
ax.legend()
ax.set_xlabel("longitude")
ax.set_ylabel("latitude")
fig

ff = 1
for ff in eachindex(θ)
    θ[ff][θ[ff] .== 0] .= NaN
    θFU[ff][θFU[ff] .== 0] .= NaN
    θFV[ff][θFV[ff] .== 0] .= NaN
end

fig, axes = plt.subplots(3, figsize = (7, 7))
vmin, vmax = nm.extrema(θ[ff][30:60,150:180])
axes[1].pcolormesh(Γ.XC[ff][30:60,150:180], Γ.YC[ff][30:60,150:180], θ[ff][30:60,150:180], 
cmap = cmo.thermal, vmin = vmin, vmax = vmax)
axes[1].set_title("Tracer Grid")
axes[2].pcolormesh(xFU[ff][30:60,150:180], yFU[ff][30:60,150:180], θFU[ff][30:60,150:180], 
cmap = cmo.thermal, vmin = vmin, vmax = vmax)
axes[2].set_title("U Grid")
axes[3].pcolormesh(xFV[ff][30:60,150:180], yFV[ff][30:60,150:180], θFV[ff][30:60,150:180], 
cmap = cmo.thermal, vmin = vmin, vmax = vmax)
axes[3].set_title("V Grid")
fig


ff = 4; k = 35
for ff in eachindex(θ)
    uθ_conv3D_approx[ff][uθ_conv3D_approx[ff] .== 0] .= NaN
    uθ_conv3D[ff][uθ_conv3D[ff] .== 0] .= NaN
end

fig, axes = plt.subplots(2, figsize = (5, 10), sharex = true)
vmin, vmax = nm.extrema(uθ_conv3D_approx[ff, k][70:110, 60:80])
axes[1].pcolormesh(xFU[ff][70:110, 60:80], yFU[ff][70:110, 60:80], uθ_conv3D_approx[ff, k][70:110, 60:80], 
cmap = cmo.thermal, vmin = vmin, vmax = vmax)
axes[1].set_title(L"\nabla (\overline{U} \overline{\theta})}")
cf = axes[2].pcolormesh(xFU[ff][70:110, 60:80], yFU[ff][70:110, 60:80], uθ_conv3D[ff, k][70:110, 60:80], 
cmap = cmo.thermal, vmin = vmin, vmax = vmax)
axes[2].set_title(L"\overline{\nabla (U \theta)}")
[ax.set_ylabel("Latitude") for ax in axes]
axes[2].set_xlabel("Longitude")
fig.colorbar(cf, ax = axes, orientation = "horizontal", pad = 0.08)
fig

println("exact flux", nm.sum(uθ_conv3D[ff, k][70:110, 60:80]))
println("approx flux", nm.sum(uθ_conv3D_approx[ff, k][70:110, 60:80]))

fig, axes = plt.subplots(1, figsize = (5, 10), sharex = true)
perc_diff = 100 * (uθ_conv3D_approx[ff, k][70:110, 60:90] .- uθ_conv3D[ff, k][70:110, 60:90]) ./ uθ_conv3D[ff, k][70:110, 60:90]
vmin, vmax = -100, 100
cf = axes.pcolormesh(xFU[ff][70:110, 60:90], yFU[ff][70:110, 60:90], perc_diff, 
cmap = cmo.thermal, vmin = vmin, vmax = vmax)
axes.set_title(L"\nabla (\overline{U} \overline{\theta})}")
fig.colorbar(cf, ax = axes, orientation = "horizontal", pad = 0.08)
fig

println("exact flux", nm.sum(uθ_conv3D[ff, k][70:110, 60:80]))
println("approx flux", nm.sum(uθ_conv3D_approx[ff, k][70:110, 60:80]))

