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
vars = ["only_init", "only_kappa", "only_sfc",  "iter0_bulkformula", "iter129_bulkformula", "only_buoyancy", "only_wind"]
for expname in vars
    Pσ = jldopen(datadir(expname * region * "_AVG_P_sigma2.jld2"))["P"]
    θσ = jldopen(datadir(expname * region * "_AVG_THETA_sigma2.jld2"))["θ"]

    avg_Pσ = mean(Pσ, dims = 2)[:]; wherenan = isfinite.(avg_Pσ .* 1)
    avg_Pσ = Pσ[:, 312][wherenan]

    p_to_z = Spline1D(pstdz, z;k=3)
    avg_zσ = p_to_z.(avg_Pσ)

    θ_dict[expname] = θσ[wherenan, :]
    z_dict[expname] = avg_zσ
end

for expname in ["only_init", "only_kappa", "only_sfc", "iter129_bulkformula", "only_buoyancy", "only_wind"]

    avg_zσ = z_dict["iter0_bulkformula"]
    new_θ = zeros(length(avg_zσ), 312)
    
    for it in 1:312 
        θ_interp = Spline1D(z_dict[expname], θ_dict[expname][:, it];k=3)
        new_θ[:, it] .= θ_interp.(avg_zσ[:])[:]
    end
    θ_dict[expname] = new_θ
    z_dict[expname] = avg_zσ
end

E, F = trend_matrices(Float32.(tecco))

fig, axs = plt.subplots( sharey = true, figsize = (15, 10))
trends_dict = Dict()
for (i, expname) in enumerate(vars)
    trends  = θ_dict[expname][:, end] .- θ_dict[expname][:, 1]
    trends = trends[:]
    trends_interp = Spline1D(z_dict[expname][3:end], trends[3:end];k=3)
    trends_dict[expname] = zeros(length(z)); 
    trends_dict[expname][15:end] .= trends_interp.(z[15:end])
    cm = axs.plot(trends_dict[expname] .* 100, z, label = expname)
    # axs.set_title(expname)
    axs.set_xlim(0, 5); 
end
axs.set_ylim(1500, 5000); 
axs.invert_yaxis()
fig

jldsave(datadir(region * "_temp_diffs_sigma2_to_z.jld2"), diffs= trends_dict)

trends_dict["iter129_bulkformula"][lvls]

trends_dict["iter0_bulkformula"]
fig