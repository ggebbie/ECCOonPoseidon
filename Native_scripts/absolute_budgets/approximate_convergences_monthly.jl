#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
using .OHC_helper
import NaNMath as nm
include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));
H = OHC_helper.smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H

expname = "iter129_bulkformula"
runpath,diagpath = listexperiments(exprootdir());

filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"trsp_3d_set2")
datafilelist_H  = filter(x -> occursin("data",x),filelist) 
filelist = searchdir(diagpath[expname],"trsp_3d_set3")
datafilelist_R  = filter(x -> occursin("data",x),filelist)
filelist = searchdir(diagpath[expname], "state_2d_set1") # first filter for state_3d_set1
datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"



lvls = findall( -3000 .<= z[:].<= -2000)

θFr = MeshArray(γ,Float32,50); fill!(θFr, 0.0)
wθ_conv_approx = MeshArray(γ,Float32,50); fill!(wθ_conv_approx, 0.0)
ctrl_vol = sum(cell_volumes[:, lvls])
nt = 35
nl = length(lvls)
uθ_conv3D = MeshArray(γ,Float32,nl);
uθ_conv3D_approx = MeshArray(γ,Float32,nl);
θs = [];
convH = []; convR = []
approx_convH = []; approx_convR = []


mskC, mskW, mskS = get_msk(Γ)

for tt = 1:nt
    println(tt)
    dθλ = γ.read(diagpath[expname]*datafilelist_H[tt],MeshArray(γ,Float32,200))
    κθx = dθλ[ :, 1:50]; κθy = dθλ[:, 51:100]
    uθx = dθλ[:, 101:150]; uθy = dθλ[:, 151:200]; dθλ = nothing 
    uθx = uθx[:, lvls]; uθy = uθy[:, lvls]
    dθr = γ.read(diagpath[expname]*datafilelist_R[tt],MeshArray(γ,Float32,150))
    wθ_conv = dθr[:, 1:50];
    κzθ_conv = dθr[ :, 51:100] .+ dθr[ :, 101:150]; dθr = nothing 
    wθ_conv_approx = MeshArray(γ,Float32,50); fill!(wθ_conv_approx, 0.0)

    # θ = OHC_helper.extract_sθ(expname,diagpath, γ, datafilelist_S[tt], datafilelist_θ[tt], inv_H)
    θ = γ.read(diagpath[expname]*datafilelist_θ[tt],MeshArray(γ,Float32, 50))
    u, v, w, Ub, Vb, Wb = OHC_helper.extract_velocities_and_bolus(diagpath, expname , 
    datafilelist_uvw[tt], γ, Γ, mskC, mskW, mskS)

    for a in eachindex(u)
        u.f[a] .= u.f[a] .+ Ub.f[a]
        v.f[a] .= v.f[a] .+ Vb.f[a]
        w.f[a] .= w.f[a] .+ Wb.f[a]
    end

    U, V = UVtoTrsp(u, v, Γ); 
    W = w .* area
    θ_lvls = θ[:, lvls] 
    push!(θs, sum(θ_lvls .* cell_volumes[:, lvls]) / ctrl_vol)
    θU, θV = OHC_helper.tofaces(θ_lvls, θ_lvls, Γ)
    U = U[:, lvls]; V = V[:, lvls]
    for a in eachindex(θU)
        θU.f[a] .= U.f[a] .* θU.f[a]
        θV.f[a] .= V.f[a] .* θV.f[a]
    end

    OHC_helper.calc_UV_conv3D!(uθx, uθy, uθ_conv3D);
    OHC_helper.calc_UV_conv3D!(θU, θV, uθ_conv3D_approx);

    for k = 1:49 
        θFr[:, k+1] = (θ[:, k] + θ[:, k+1]) ./2 #linear interpolate to faces
        wθ_conv_approx[:, k+1] = θFr[:, k+1] .* W[:, k+1]
    end

    wθ_conv = wθ_conv[:, lvls[end]+1] .- wθ_conv[:, lvls[1]] #in - out 
    wθ_conv_approx = wθ_conv_approx[:, lvls[end]+1] .- wθ_conv_approx[:, lvls[1]] 

    push!(approx_convR, sum(wθ_conv_approx .* PAC_msk) / ctrl_vol)
    push!(convR, sum(wθ_conv .* PAC_msk) / ctrl_vol)
    push!(approx_convH, sum(uθ_conv3D_approx .* PAC_msk) / ctrl_vol)
    push!(convH, sum(uθ_conv3D .* PAC_msk) / ctrl_vol)
end

tecco = 1992+1/24:1/12:2018
tecco = tecco[1:nt]
import PyPlot as plt
fig, axes = plt.subplots(2, 3, figsize = (15, 7))

axes[1].plot(tecco, convH, label = L"\overline{\nabla (U \theta)}", c = "blue")
axes[1].plot(tecco, approx_convH, label = L"\nabla (\overline{U} \overline{\theta})}", c = "orange")
axes[1].plot(tecco, (convH .- approx_convH), label = L"Cov(U, \theta)", c = "green")

axes[1].set_ylabel("[°C / s]")
axes[1].set_title("Horizontal Advective Convergence")
axes[2].plot(tecco, cumsum(convH) .* 2.628e+6, label = L"\int{\nabla \overline{(U \theta)} dt}", c = "blue")
axes[2].plot(tecco, cumsum(approx_convH) .* 2.628e+6, label = L"\int{\nabla (\overline{U} ~ \overline{\theta})} dt}", c = "orange")
axes[2].set_title("Integrated Horizontal Advective Convergence")
axes[2].set_ylabel("[°C]")

axes[3].plot(tecco, convR, label = L"\nabla \overline{(W \theta)}", c = "blue")
axes[3].plot(tecco, approx_convR, label = L"\nabla (\overline{W} ~ \overline{\theta})}", c = "orange")
axes[3].plot(tecco, convR .- approx_convR, label = L"Cov(W, \theta)", c = "green")

axes[3].set_ylabel("[°C / s]")
axes[3].set_title("Vertical Advective Convergence")
axes[4].plot(tecco, cumsum(convR) .* 2.628e+6, label = L"\int{\nabla \overline{(W \theta)} dt}", c = "blue")
axes[4].plot(tecco, cumsum(approx_convR) .* 2.628e+6, label = L"\int{\nabla (\overline{W} ~ \overline{\theta})} dt}", c = "orange")
axes[4].set_title("Integrated Vertical Advective Convergence")
axes[4].set_ylabel("[°C]")

axes[5].plot(tecco, (convH .+ convR), label = L"\nabla \overline{((W + U) \theta)}", c = "blue")
axes[5].plot(tecco, (approx_convH .+ approx_convR), 
label = L"\nabla ((\overline{W} + \overline{U}) \overline{\theta})}", c = "orange")
axes[5].plot(tecco, ((convH .- approx_convH) .+ (convR .- approx_convR)), label = L"Cov(W + U, \theta)", c = "green")
axes[5].set_ylabel("[°C / s]")
axes[5].set_title("Advective Convergence")

axes[6].plot(tecco, cumsum( (convH .+ convR)) .* 2.628e+6, 
label = L"\int{\nabla \overline{(W \theta + U \theta)} dt}", c = "blue")
axes[6].plot(tecco, cumsum((approx_convH .+ approx_convR)) .* 2.628e+6, 
label = L"\int{ \nabla ((\overline{W} + \overline{U}) \overline{\theta}) dt }", c = "orange")
axes[6].set_title("Integrated Advective Convergence")
axes[6].set_ylabel("[°C]")

[ax.legend(ncol = 3, fontsize = "small",handlelength = 0.3, loc = "lower left") for ax in axes]
axes[4].legend(ncol = 3, fontsize = "small",handlelength = 0.3, loc = "lower right")
[ax.set_ylim(-5e-10, 5e-10) for ax in axes[1, :]]
[ax.set_xlabel("time") for ax in axes]
fig.tight_layout()
# fig.savefig("TestConvergences_notfull.png")
fig