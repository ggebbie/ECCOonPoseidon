include("../../src/intro.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
JLD2, Printf, LaTeXStrings, PyPlot
using Dierckx, Interpolations

include(srcdir("config_exp.jl"))
# include("./IsopycnalHelpers.jl")

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "NPAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region)

area_mask = area .* PAC_msk
runpath,diagpath = listexperiments(exprootdir());

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)
sig1grid = sigma2grid("NPAC")
nσ = length(sig1grid)
TSroot = "p_on_sigma2" 

sig2dirs = Dict()
sig2dirs["iter0_bulkformula"] = rundir("iter0_bulkformula")*"sigma2/"
sig2dirs["only_init"] = vastrundir("nosfcadjust_exps", "run_adjust_init")*"sigma2/"
sig2dirs["only_kappa"] = vastrundir("nosfcadjust_exps", "run_adjust_kappa")*"sigma2/"
sig2dirs["only_sfc"] = vastrundir("nooceanadjust", "run_noadjusts")*"sigma2/"
sig2dirs["only_buoyancy"] = vastrundir("nooceanadjust", "run_noadjusts_nowind")*"sigma2/"
sig2dirs["only_wind"] = vastrundir("nooceanadjust", "run_noadjusts_nobuoyancy")*"sigma2/"
#read in the first time step of S and θ
for expname in ["only_buoyancy", "only_wind"]
    println(expname)
    # Get list of files for salinity on sigma1
    
    filelist = searchdir(sig2dirs[expname], TSroot)
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(datafilelist)

    #load trend trend_matrices, F is LS estimator
    P = zeros(Float32, nσ, nt);

    @time for tt = 1:nt
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        Tname = datafilelist[tt]
        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        @time Pσ1 = γ.read(sig2dirs[expname]*Tname,MeshArray(γ,Float32,nσ))
        area_3D = MeshArray(γ,Float32,nσ); fill!(area_3D, 0.0)
        for a in eachindex(Pσ1)
            nan_pts = isnan.(Pσ1.f[a] .* 1)
            area_3D.f[a] .= area_mask.f[a[1]] .* (!).(nan_pts)
            Pσ1.f[a][nan_pts] .= 0.0
        end
        P[:, tt] .= lateral_sum(Pσ1 .* area_3D) ./ lateral_sum(area_3D)
    end

    jldsave(datadir(expname * region * "_AVG_P_sigma2.jld2"), P= P)
end

fig, ax = plt.subplots()
cm = ax.contourf(tecco, sig2grid, P, vmin = 0, vmax = 5500, cmap = "Spectral", levels = 35)
fig.colorbar(cm)
ax.set_ylim(36.75, 37); ax.invert_yaxis()
fig