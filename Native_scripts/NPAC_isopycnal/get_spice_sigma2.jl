include("../../src/intro.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
JLD2, Printf, LaTeXStrings, PyPlot
using Dierckx, Interpolations
include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "NPAC"

runpath,diagpath = listexperiments(exprootdir());

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)

sig2grid = sigma2grid("NPAC")
nσ = length(sig2grid)
# σlvls = findall( σtop .<= sig1grid .<= σbot)

θ_dict = Dict(); z_dict = Dict()
for expname in ["only_init", "only_kappa", "only_sfc",  "iter0_bulkformula", "iter129_bulkformula"]
    Pσ = jldopen(datadir(expname * region * "_AVG_P_sigma2.jld2"))["P"]
    θσ = jldopen(datadir(expname * region * "_AVG_THETA_sigma2.jld2"))["θ"]

    avg_Pσ = mean(Pσ, dims = 2)[:]; wherenan = isfinite.(avg_Pσ .* 1)
    avg_Pσ = Pσ[:, 312][wherenan]

    p_to_z = Spline1D(pstdz, z;k=3)
    θ_dict[expname] = θσ[wherenan, :]
    z_dict[expname] = p_to_z.(avg_Pσ)
end

for expname in ["iter0_bulkformula", "only_init", "only_kappa", "only_sfc", "iter129_bulkformula"]
    println(expname)
    new_θ = zeros(length(z), 312)
    
    for it in 1:312 
        θ_interp = Spline1D(z_dict[expname][5:end], θ_dict[expname][5:end, it][:];k=3)
        new_θ[:, it] .= θ_interp.(z)[:]
    end
    θ_dict[expname] = new_θ
    z_dict[expname] = z
end

jldsave(datadir(region * "_temp_sigma2_to_z.jld2"), θ_dict= θ_dict, z_dict = z_dict)

# [θ_dict[expname] .-= θ_dict["iter0_bulkformula"] for expname in ["only_init", "only_kappa", "only_sfc"]]
# θ_dict["SUM"] = 1 .* θ_dict["iter0_bulkformula"]
# [θ_dict["SUM"] .+= θ_dict[expname] for expname in ["only_init", "only_kappa", "only_sfc"]]
# z_dict["SUM"] = z_dict["iter0_bulkformula"]

fig, axs = plt.subplots(2, 3, sharey = true, figsize = (15, 10))
for (i, expname) in enumerate(["iter0_bulkformula", "only_init", "only_kappa", "only_sfc", "iter129_bulkformula", "SUM"])
    z_levs = 1500 .< z_dict[expname] .< 3500
    θ = 100 .* θ_dict[expname][z_levs, :]
    cm = axs[i].contourf(tecco, z_dict[expname][z_levs], θ .- mean(θ, dims = 2), 
    cmap = cmo.balance, vmin = -0.5, vmax = 0.5, levels = 25, extend = "both")
    axs[i].set_title(expname)
end
axs[1].set_ylim(1750, 3250); 

axs[1].invert_yaxis()

fig

fig, axs = plt.subplots(sharey = true, figsize = (15, 10))
exp_colors["SUM"] = "orange"
for (i, expname) in enumerate(["iter0_bulkformula", "only_init", "only_kappa", "only_sfc", "iter129_bulkformula"])
    z_levs = 2000 .< z_dict[expname] .< 3000
    θ = mean(θ_dict[expname][z_levs, :], dims = 1)
    axs.plot(tecco, θ[:],label = expname, color= exp_colors[expname])
    # axs[i].set_title(expname)
end
axs.legend()
fig

mean(θ_dict["only_init"][z_levs, :], dims = 1)