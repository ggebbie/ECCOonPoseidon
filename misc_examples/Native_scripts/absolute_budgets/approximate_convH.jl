#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
using .OHC_helper
import NaNMath as nm
include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));
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

uθ_conv3D = MeshArray(γ,Float32,50);
uθ_conv3D_approx = MeshArray(γ,Float32,50);
convH = []
approx_convH = []
for tt = 1:312
    print(tt)
    dθλ = γ.read(diagpath[expname]*datafilelist_H[tt],MeshArray(γ,Float32,200))
    κθx = dθλ[ :, 1:50]; κθy = dθλ[:, 51:100]
    uθx = dθλ[:, 101:150]; uθy = dθλ[:, 151:200]; dθλ = nothing 
    OHC_helper.calc_UV_conv3D!(uθx, uθy, uθ_conv3D);

    θ = γ.read(diagpath[expname]*datafilelist_θ[tt],MeshArray(γ,Float32, 50))
    UV = γ.read(diagpath[expname]*datafilelist_uvw[tt],MeshArray(γ,Float32,100))
    u = UV[:, 1:50]; v = UV[:, 51:100]
    U, V = UVtoTrsp(u, v, Γ)

    θU, θV = MeshArray(γ,Float32,50), MeshArray(γ,Float32,50)

    for k = 1:50 
        tmp = θ[:, k]
        tmp1, tmp2 = tofaces(tmp, tmp)
        θU[:, k] = U[:, k] .* tmp1
        θV[:, k] = V[:, k] .* tmp2
    end
    OHC_helper.calc_UV_conv3D!(θU, θV, uθ_conv3D_approx);
    push!(approx_convH, sum(uθ_conv3D_approx[:, 39] .* PAC_msk))
    push!(convH, sum(uθ_conv3D[:, 39] .* PAC_msk))
end

ctrl_vol = sum(cell_volumes[:, 39])
tecco = 1992+1/24:1/12:2018
import PyPlot as plt
fig, axes = plt.subplots(2, 1, figsize = (9, 5))
# mean_residual = mean((convH .- approx_convH))
mean_residual = 0
axes[1].plot(tecco, convH ./ ctrl_vol, label = L"\overline{\nabla (U \theta)}", c = "blue")
axes[1].plot(tecco, approx_convH ./ ctrl_vol, label = L"\nabla (\overline{U} \overline{\theta})}", c = "orange")
axes[1].legend()
axes[1].set_ylabel("[°C / s]")
axes[1].set_title("Horizontal Advective Convergence")
axes[2].plot(tecco, cumsum(convH ./ ctrl_vol) .* 2.628e+6, label = L"\int{\overline{\nabla (U \theta)} dt}", c = "blue")
axes[2].plot(tecco, cumsum((approx_convH .+ mean_residual)./ ctrl_vol) .* 2.628e+6, label = L"\int{\nabla (\overline{U} \overline{\theta})} dt}", c = "orange")
axes[2].legend()
axes[2].set_title("Integrated Horizontal Advective Convergence")
axes[2].set_ylabel("[°C]")
[ax.set_xlabel("time") for ax in axes]
fig.tight_layout()
fig

fig, axes = plt.subplots(figsize = (9, 5))
zF = (Γ.RC[1:end-1] .+ Γ.RC[2:end]) ./ 2
zFtrue = Γ.RF[2:end-1]
resid = zF .- zFtrue
axes.scatter(zF .- zFtrue, 1:length(zF))
axes.invert_yaxis()
fig

