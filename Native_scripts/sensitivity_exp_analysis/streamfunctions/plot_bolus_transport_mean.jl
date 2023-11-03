include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall
import PyPlot as plt 
import NaNMath as nm
include(srcdir("config_exp.jl"))

tecco = 1992+1/24:1/12:2018
nz = 50

@pyimport cmocean.cm as cmo


(ϕ,λ) = latlonC(γ)
area = readarea(γ)
λ_wrap = wrap_λ(λ)
ocean_mask = wet_pts(Γ)

nanmean(x) = mean(filter(!isnan,x))
nanmean(x,y) = mapslices(nanmean,x,dims=y)

function velocity2center3D(u,v)

    uC=similar(u); vC=similar(v)

    (tmpU,tmpV)=exch_UV_llc90(u,v); 


    for a in eachindex(uC)
        (s1, s2) = size(u.f[a])
        tmpU1 = view(tmpU.f[a], 1:s1, 1:s2)
        tmpU2 = view(tmpU.f[a], 2:s1 + 1, 1:s2)
        tmpV1 = view(tmpV.f[a], 1:s1, 1:s2)
        tmpV2 = view(tmpV.f[a], 1:s1, 2:s2 + 1)

        uC.f[a] .= tmpU1 .+ tmpU2
        vC.f[a] .= reshape(nanmean([tmpV1[:] tmpV2[:]],2),size(tmpV1))
    end

    return uC, vC
end

include(srcdir("plot_and_dir_config.jl"))
include(srcdir("MeshArraysPlots.jl"))

open_velocities(fname) = (jldopen(fname)["Eul_Vel"], jldopen(fname)["Bol_Vel"])
vars =  ["only_init", "only_kappa", "only_sfc", "iter129_bulkformula",  "iter0_bulkformula"]

EKE_dict = Dict()

lvls = findall(2000 .<= z[:] .<= 3000)

@time for expname in vars
    fname = datadir("native/"* expname * "_mean_Eul_Bol_Trsp.jld2")
    Eul_Vel, Bol_Vel = open_velocities(fname)
    Ubol_interp, Vbol_interp = interpolate_to_lateral_faces(abs.(Bol_Vel["U"]), abs.(Bol_Vel["V"]), Γ)

    EKE = similar(Ubol_interp); fill!(EKE, 0.0)
    for a in eachindex(EKE)
        EKE.f[a] .= (Ubol_interp.f[a] .* (86400 * 365)).^2 .+ (Vbol_interp.f[a] .* (86400* 365)).^2
        EKE.f[a] .= 0.5 * EKE.f[a] .* Γ.DRF[a[2]]
    end
    EKE_dict[expname] = vertical_sum(EKE[:, lvls])
end

[EKE_dict[ex] .= EKE_dict[ex] .- EKE_dict["iter0_bulkformula"] for ex in ["only_init", "only_kappa", "only_sfc"]]
EKE_dict["SUM"] = EKE_dict["iter0_bulkformula"] .+ EKE_dict["only_init"] .+ EKE_dict["only_kappa"] .+ EKE_dict["only_sfc"]
vars =  ["only_init", "iter0_bulkformula", "only_kappa", "iter129_bulkformula", "only_sfc",  "SUM"]

fig, axs = plt.subplots(2, 3, figsize=(20,10), subplot_kw=Dict("projection"=> proj0))
[a.set_extent((120, 285, -40, 70),crs=projPC) for a in axs]
[a.coastlines() for a in axs]
vmax = 40
plot_labels_effects["SUM"] = "SUM"
CB = Any[]
for (i, expname) in enumerate(vars)
    data = (1040  * 1000 * 1e-15) .* EKE_dict[expname] 
    cb = pcolormesh_ma(axs[i], λ, ϕ, data, 
    cmo.curl; bounds = (-vmax, vmax), add_gridlines = true)
    axs[i].set_title(expname); push!(CB, cb)
    axs[i].set_title(plot_labels_effects[expname])
end

fig.colorbar(CB[1], ax = axs[:], orientation = "horizontal", fraction = 0.05, label = "petagrams per year²"); 
fig.savefig(plotsdir("native/sensitivity_exps/Vertically_Int_MidDepth_EKE.png"))
