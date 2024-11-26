include("../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
using NaNMath
import PyPlot as plt

include(srcdir("plot_and_dir_config.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)

tecco= 1992+1/24:1/12:2018 # ecco years
nz = length(z)
nt = length(tecco)
P = MeshArray(γ,Float32, nz)

for ijk in eachindex(P)
    P[ijk] .= -pstdz[ijk[2]]
end

p₀ = 2000

expname = "iter129_bulkformula"
filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
datafilelist_η = filter(x -> occursin("data",x),filelist) # second filter for "data"


zonal_avg = Dict()
P_hyd = MeshArray(γ,Float32, nz); fill!(P_hyd, 0.0)
θz_mean = MeshArray(γ,Float32, nz); fill!(θz_mean, 0.0)
σz_mean = MeshArray(γ,Float32, nz); fill!(σz_mean, 0.0)
Sz_mean = MeshArray(γ,Float32, nz); fill!(Sz_mean, 0.0)
η_mean = MeshArray(γ,Float32); fill!(η_mean, 0.0)


nt = 36
for tt in 1:nt
    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    θname = datafilelist_θ[tt]
    ηname = datafilelist_η[tt]
    # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
    @time θSz = γ.read(diagpath[expname]*θname,MeshArray(γ,Float32,2*nz))
    @time ηz = γ.read(diagpath[expname]*ηname,MeshArray(γ,Float32))
    η_mean.f[:] .+= (ηz.f[:] ./ nt)
    θz = θSz[:, 1:nz]; Sz = θSz[:, nz+1:end]
    for ijk in eachindex(P)
        θz_mean[ijk] .+= θz.f[ijk] ./ nt
        Sz_mean[ijk] .+= Sz.f[ijk] ./ nt

        σ = densityJMD95.(θz.f[ijk],Sz.f[ijk], P[ijk], p₀) #EOS from MITGCM 

        σz_mean[ijk] .+= σ ./ nt
        P_hyd[ijk] .+= (σ .* 9.81 .* Γ.DRC[ijk[2]]) ./ nt

        P_hyd[ijk][iszero.(Γ.hFacC[ijk])] .= 0.0

    end 
end

fig, ax = plt.subplots()
ax.pcolormesh(1:270, -z, zonal_average(σz_mean, cell_volumes))
fig


P_hyd_int= MeshArray(γ,Float32, nz); fill!(P_hyd_int, 0.0)
P_hyd_int.f[:, 1] .= 1 .* P_hyd.f[:, 1]
# η_mean .* 
for k in 2:50
    P_hyd_int.f[:, k] = P_hyd.f[:, k] .+ P_hyd_int.f[:, k-1]
    
end
P_hyd_arr = zonal_average(P_hyd, cell_volumes)
P_hyd_arr = cumsum(P_hyd_arr, dims = 1)
P_hyd_arr[:, 130]

fig, ax = plt.subplots()
ax.pcolormesh(1:270, -z, zonal_average(P_hyd_int, cell_volumes))
fig

fig, ax = plt.subplots(1, 2, sharey = true)
ax[1].plot(zonal_average(P_hyd_int, cell_volumes)[:, 130] .* 1e-4, -z)
ax[1].plot(pstdz, -z)
ax[2].plot(zonal_average(P_hyd_int, cell_volumes)[:, 150] .* 1e-4, -z)
ax[2].plot(pstdz, -z)
fig


fig, ax = plt.subplots(1, 2, sharey = true)
ax[1].plot(P_hyd_arr[:, 130] .* 1e-4, -z)
ax[1].plot(pstdz, -z)
ax[2].plot(P_hyd_arr[:, 135] .* 1e-4, -z)
ax[2].plot(pstdz, -z)
fig

zonal_average(P_hyd_int, cell_volumes)[ :, 140:end]

fig, ax = plt.subplots(sharey = true)
ax.pcolormesh(λ[5], ϕ[5], P_hyd_int[5, 46] .* PAC_msk[5])
fig

fig, ax = plt.subplots(sharey = true)
ax.pcolormesh(λ[5], ϕ[5], η_mean[5] .* PAC_msk[5])
fig


zonal_average(σz_mean, cell_volumes)[:, 150]
zonal_average(P_hyd_int, cell_volumes)[:, 150]


P_hyd_int= MeshArray(γ,Float32, nz); fill!(P_hyd_int, 0.0)
P_hyd_int.f[:, 1] .= 1 .* P_hyd.f[:, 1]
for ff = 1:5
    P_hyd_int.f[ff, 1] .+= (η_mean.f[ff] .* σz_mean.f[ff, 1] .* 9.81)
end
for k in 2:50
    P_hyd_int.f[:, k] = P_hyd.f[:, k] .+ P_hyd_int.f[:, k-1]
end


fig, ax = plt.subplots(sharey = true)
ax.pcolormesh(λ[5], ϕ[5], P_hyd_int[5, 46] .* PAC_msk[5])
fig

fig, ax = plt.subplots(1, 2, sharey = true)
ax[1].plot(zonal_average(P_hyd_int, cell_volumes)[:, 130] .* 1e-4, -z)
ax[1].plot(pstdz, -z)
ax[2].plot(zonal_average(P_hyd_int, cell_volumes)[:, 120] .* 1e-4, -z)
ax[2].plot(pstdz, -z)
fig

