include("../../src/intro.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
JLD2, Printf, LaTeXStrings, PyPlot, GibbsSeaWater
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

θ_dict = Dict(); z_dict = Dict(); θσ_dict = Dict()
vars = ["iter0_bulkformula", "iter129_bulkformula", "only_buoyancy", 
"only_init", "only_kappa", "only_wind"]

for expname in vars
    Pσ = jldopen(datadir(expname * region * "_AVG_P_sigma2.jld2"))["P"]
    θσ = jldopen(datadir(expname * region * "_AVG_THETA_sigma2.jld2"))["θ"]

    # avg_Pσ = mean(Pσ, dims = 2)[:]; wherenan = isfinite.(avg_Pσ .* 1)
    avg_Pσ = Pσ[:, 1]; wherenan = isfinite.(avg_Pσ .* 1)
    avg_Pσ = avg_Pσ[wherenan]

    avg_zσ = gsw_z_from_p.(avg_Pσ, 45, 0, 0)

    θ_dict[expname] = θσ[wherenan, :]
    θ_dict[expname][isnan.(θ_dict[expname])] .= 0.0
    z_dict[expname] = avg_zσ
end

fig, ax = plt.subplots()
ax.scatter(1:53, -z_dict["iter0_bulkformula"])
fig
vars = ["iter0_bulkformula", "iter129_bulkformula", "only_buoyancy", 
"only_init", "only_kappa", "only_wind"]
zs = z[10:end-3]
for expname in vars
    print(expname)
    new_θ = zeros(length(zs), 312)
    sorted_iz = sortperm(-z_dict[expname])
    for it in 1:312 
        # θ_interp = Spline1D(-z_dict[expname], θ_dict[expname][:, it];k=3, bc = "error")
        # new_θ[:, it] .= θ_interp.(zs[:])[:]
        θ_interp = interpolate((-Float32.(z_dict[expname])[sorted_iz],), θ_dict[expname][sorted_iz, it], Gridded(Linear()))
        new_θ[:, it] = θ_interp(zs[:])[:]
    end
    θ_dict[expname] = new_θ
    z_dict[expname] = zs
end

# unique(Float32.(z_dict["iter129_bulkformula"]))
# z_dict["iter0_bulkformula"]
jldsave(datadir(region * "_temp_sigma2_to_z_Nov.jld2"), θ_dict= θ_dict, z_dict = z_dict)


E, F = trend_matrices(Float32.(tecco))

fig, axs = plt.subplots( sharey = true, figsize = (15, 10))
trends_dict = Dict()
expname = "iter129_bulkformula"
trends  = θ_dict[expname][:, end] .- θ_dict[expname][:, 1]

for (i, expname) in enumerate(vars)
    trends  = θ_dict[expname][:, end] .- θ_dict[expname][:, 1]
    trends = trends[:]
    cm = axs.plot(trends .* 100, zs, label = expname)
    # axs.set_title(expname)
    axs.set_xlim(0, 5); 
end
axs.set_ylim(-5000, -1500); 
axs.invert_yaxis()
fig

jldsave(datadir(region * "_temp_diffs_sigma2_to_z.jld2"), diffs= trends_dict)

trends_dict["iter129_bulkformula"][lvls]

trends_dict["iter0_bulkformula"]
fig