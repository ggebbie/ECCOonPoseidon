using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools,
PyPlot, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions
using PyPlot   # important!
using PyCall
@pyimport seaborn as sns
sns.set(); pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;
function ref_month_1(data, ref_key)
    baseline = data[ref_key][1]
    temp_dict = Dict(key => (data[key] .- baseline) for key in keys(data))
    return temp_dict
end    
function remove_seasonal_keep_mean(x, Ecycle, Fcycle )
    temp_mean = mean(x)
    x_noszn = remove_seasonal(x, Ecycle, Fcycle)
    return x_noszn .+ temp_mean
end

θbar = Dict()
θbarLHS, θbarRHS = Dict(),Dict()
reconstruct(x0, dx) = cumsum(vcat(Float32(x0), dx .* Float32(Δt₂)))
fcycle = 1 # units: yr^{-1}
overtones = 2
E,F = trend_matrices(tecco)
Ecycle,Fcycle = seasonal_matrices(fcycle,tecco,overtones)

fig_deep, axs_deep = plt.subplots(2, 1, figsize = (7, 12))
fig_deep.suptitle("Deep θ̄ Trends")
axs_deep[1].set_title("θ̄ Trend (Using Least-Squares Fit)")
axs_deep[2].set_title("θ̄ Trend (Using Finite Difference)")

fig_surf, axs_surf = plt.subplots(2, 1, figsize =(7, 12))
fig_surf.suptitle("Surface θ̄ Trends")
axs_surf[1].set_title("θ̄ Trend (Using Least-Squares Fit)")
axs_surf[2].set_title("θ̄ Trend (Using Finite Difference)")
filedir = "ECCO_vars/"
surf_fname = "Surface_STHETA"
z = []
fig_scat, axs_scat = plt.subplots(1, 1, figsize = (7, 7))
axs_scat.set_xlabel("Deep θ Trend " * L"(º C / century)")
axs_scat.set_ylabel("Surface θ Trend"* L"(º C / yr)")

fig, axs = plt.subplots(1, 1, figsize = (12, 7))
fig_smooth, axs_smooth = plt.subplots(2, 1, figsize = (12, 7))
fig_smooth.suptitle("Interannual Temperature Anomalies")
axs_smooth[1].set_title("Surface θ'")
axs_smooth[2].set_title("Deep θ'")

pts = zeros(2, 4)
θsurf_noszn = Dict()
for (i, key) in enumerate(keys(shortnames))
    expname = key
    println("loading surface θ")
    @time @load datadir(filedir*surf_fname*"_"*expname* ".jld2") var_exp
    θ_surf = []
    for tt in 1:length(tecco)
        push!(θ_surf, volume_mean(var_exp[tt]; weights = cell_volumes[:, 1]))
    end

    θ_surf = remove_seasonal_keep_mean(θ_surf,Ecycle,Fcycle)
    θsurf_noszn[key] = θ_surf
    θ_deep = remove_seasonal_keep_mean(θz[key],Ecycle,Fcycle)

    axs_smooth[1].plot(tecco, θ_surf, label = key)
    axs_smooth[2].plot(tecco, θ_deep, label = key)
    βsurf = get_trend(θ_surf,tecco,F)
    dθsurf=(θ_surf[end]-θ_surf[1])/(tecco[end]-tecco[1])
    axs_surf[1].bar(i, βsurf, label = key)
    axs_surf[2].bar(i, dθsurf, label = key)

    βdeep = 1f2*get_trend(θ_deep,tecco,F)
    dθdeep=1f2*(θ_deep[end]-θ_deep[1])/(tecco[end]-tecco[1])
    axs_deep[1].bar(i, βdeep, label = key)
    axs_deep[2].bar(i, dθdeep, label = key)

    #need to do this plot with Qnet
    pts[:, i] = [βdeep; βsurf]
    axs_scat.scatter(βdeep, βsurf, label = key)
end
for ax in [axs_surf, axs_deep, axs_smooth], i in 1:2
    ax[i].legend()
end

c, β = lin_reg(pts[1, :], pts[2, :])
min_x, max_x = extrema(pts[1, :])
x_range = collect(min_x:0.001:max_x)
axs_scat.plot(x_range, c .+ (x_range .* β), c = "black", linestyle = "--",
label = "slope = " *string(round(β, digits = 3)))

axs_scat.legend()

fig_deep.tight_layout()
fig_deep.savefig(plotsdir() * "/OHC_Divergence/" * "DeepθTrends_" * region * suffix * "_2km3km.png")
fig_surf.tight_layout()
fig_surf.savefig(plotsdir() * "/OHC_Divergence/" * "SurfθTrends_" * region * suffix * ".png")
fig_scat.tight_layout()
fig_scat.savefig(plotsdir() * "/OHC_Divergence/" * "ScatθTrends_" * region * suffix * "_2km3km.png")
fig_smooth.tight_layout()
fig_smooth.savefig(plotsdir() * "/OHC_Divergence/" * "θAnomalies_" * region * suffix * ".png")

#Now see if there is any correlation 
#(plot surface warming as a trends as a function of deep cooling)