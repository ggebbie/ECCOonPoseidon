include("../../src/intro.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
JLD2, Printf, LaTeXStrings, PyPlot

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
sig1grid = ECCOtour.sigma2grid()
nσ = length(sig1grid)
TSroot = "THETA_on_sigma2" 

sig2dirs = Dict()
sig2dirs["iter0_bulkformula"] = vastrundir("iter0_bulkformula")*"sigma2/"
sig2dirs["iter129_bulkformula"] = vastrundir("iter129_bulkformula")*"sigma2/"
sig2dirs["only_init"] = vastrundir("nosfcadjust_exps", "run_adjust_init")*"sigma2/"
sig2dirs["only_kappa"] = vastrundir("nosfcadjust_exps", "run_adjust_kappa")*"sigma2/"
sig2dirs["only_buoyancy"] = vastrundir("nooceanadjust", "run_noadjusts_nowind")*"sigma2/"
sig2dirs["only_wind"] = vastrundir("nooceanadjust", "run_noadjusts_nobuoyancy")*"sigma2/"
#read in the first time step of S and θ
for expname in keys(sig2dirs)
    println(expname)
    # Get list of files for salinity on sigma1
    filelist = searchdir(sig2dirs[expname], TSroot)
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(datafilelist)

    #load trend trend_matrices, F is LS estimator
    θ = zeros(Float32, nσ, nt);

    @time for tt = 1:312
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        Tname = datafilelist[tt]
        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        @time θσ2 = γ.read(sig2dirs[expname]*Tname,MeshArray(γ,Float32,nσ))
        area_3D = MeshArray(γ,Float32,nσ); fill!(area_3D, 0.0)
        for a in eachindex(θσ2)
            nan_pts = isnan.(θσ2.f[a] .* 1)
            area_3D.f[a] .= area_mask.f[a[1]] .* (!).(nan_pts)
            θσ2.f[a][nan_pts] .= 0.0
        end
        θ[:, tt] .= lateral_sum(θσ2 .* area_3D) ./ lateral_sum(area_3D)
    end

    jldsave(datadir(expname * region * "_AVG_THETA_sigma2.jld2"), θ= θ)
end
fig, ax = plt.subplots()
cm = ax.pcolormesh(tecco, sig2grid, θ .- mean(θ, dims = 2), 
cmap = "bwr", vmin = -0.02, vmax = 0.02)
fig.colorbar(cm)

ax.set_ylim(36.75, 37); ax.invert_yaxis()
fig

fig, ax = plt.subplots()
ax.plot(tecco, θ[sig2grid .== 36.75, :][:])
fig
fig.colorbar(cm)


ax.set_ylim(36.75, 37); ax.invert_yaxis()
fig
