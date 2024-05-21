include("../../../../src/intro.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
JLD2, Printf, LaTeXStrings, PyPlot, PyCall
using Dierckx, Interpolations
using GibbsSeaWater
import NaNMath as nm
include(srcdir("config_exp.jl"))
@pyimport cmocean.cm as cmo
(ϕ,λ) = latlonC(γ)
area = readarea(γ)
z
region = "NPAC"

runpath,diagpath = listexperiments(exprootdir());

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)

sig2grid = ECCOtour.sigma2grid()
nσ = length(sig2grid)
# σlvls = findall( σtop .<= sig1grid .<= σbot)

θ_dict = Dict(); z_dict = Dict()
for expname in ["iter0_bulkformula", "iter129_bulkformula"]
    Pσ = jldopen(datadir(expname * region * "_AVG_P_sigma2.jld2"))["P"]
    θσ = jldopen(datadir(expname * region * "_AVG_THETA_sigma2.jld2"))["θ"]

    avg_Pσ = mean(1 .* Pσ[:, 1:12], dims = 2)[:]; wherenan = isfinite.(avg_Pσ .* 1)
    # avg_Pσ = avg_Pσ[wherenan]

    p_to_z = Spline1D(pstdz, z;k=3)
    θ_dict[expname] = θσ[wherenan, :]
    z_dict[expname] = abs.(gsw_z_from_p.(avg_Pσ[wherenan], 45., 0., 0.))
    # z_dict[expname] = p_to_z.(Pσ[wherenan, :])
    # z_dict[expname] = p_to_z.(avg_Pσ[wherenan])

end

θ_dict_new = Dict();
for expname in ["iter0_bulkformula", "iter129_bulkformula"]
    println(expname)
    new_θ = zeros(length(z), 312)
    
    for it in 1:312 
        θ_interp = Spline1D(z_dict[expname][10:end], θ_dict[expname][10:end, it][:];k=3)
        new_θ[:, it] .= θ_interp.(z)[:]
    end
    θ_dict_new[expname] = new_θ
    z_dict[expname] = z
end
# θσ = jldopen(datadir("iter0_bulkformula" * region * "_AVG_P_sigma2.jld2"))["P"]
# nm.maximum(θσ)
jldsave(datadir(region * "_temp_sigma2_to_z_reference_middle.jld2"), θ_dict= θ_dict_new, z_dict = z_dict)

# [θ_dict[expname] .-= θ_dict["iter0_bulkformula"] for expname in ["only_init", "only_kappa", "only_sfc"]]
# θ_dict["SUM"] = 1 .* θ_dict["iter0_bulkformula"]
# [θ_dict["SUM"] .+= θ_dict[expname] for expname in ["only_init", "only_kappa", "only_sfc"]]
# z_dict["SUM"] = z_dict["iter0_bulkformula"]

fig, ax = plt.subplots()
ax.plot(θ_dict_new["iter0_bulkformula"][40, :])
fig