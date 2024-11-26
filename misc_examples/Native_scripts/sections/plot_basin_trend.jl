using DataFrames, Distributions, Statistics
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= -z[:].<= uplvl)

region = "NPAC"
fnames(section, expname) = datadir(region * "_" * expname * "_" * section * "_THETA_levels" * ".jld2")

extract_θ_ΔV(section, expname, lvls) = (load(fnames(section, expname))["θ"][lvls, :], 
                                  load(fnames(section, expname))["ΔV"][lvls])
volume_mean(x, ΔV) =  vec(sum(x .* ΔV, dims =1) ./ sum(ΔV))

mid_depth_temperatures(section, expname, lvls) = volume_mean(extract_θ_ΔV(section, expname, lvls)...)

expname = "iter129_bulkformula"
WOCE_PLines = [key for key in keys(sections)]
temperatures = Dict()
[temperatures[key] = mid_depth_temperatures(key, expname, lvls) for key in keys(sections)]
fig,axs=plt.subplots(1,1, sharey = true, figsize = (10, 7.5))

for key in keys(temperatures)
    axs.plot(tecco,temperatures[key], label = key)
end
axs.legend()
fig

nkeys = length(keys(sections))
θz = zeros(nkeys, 312)
[θz[i, :] .= temperatures[key] for (i, key) in enumerate(keys(sections))]
θz_mean = mean(θz, dims = 1)[:]
est_err = std(θz, dims = 1)[:]
est_err2 = sqrt.(mean((θz .- θz_mean').^2,  dims = 1))[:]

θz = mean(θz, dims = 1)[:]

E,F = trend_matrices(tecco)
a, b = F*θz
fit_line(x, a, b) = b .* (x .- mean(x)) .+ a
y0b = fit_line(tecco, a, b)
σ = std(θz .- y0b) .+ est_err2

fig,axs=plt.subplots(1,1, sharey = true, figsize = (10, 7.5))
axs.plot(tecco, θz, label = "ECCO Output");
axs.plot(tecco, y0b, label = "Trend = " * string(round(1000 * b, digits = 2)) * " m°C per year");
axs.fill_between(tecco, θz .- (1 .*σ), θz .+ (1 .*σ), alpha=0.3, color="k", label = "2σ estimate")
axs.legend(frameon=false)
axs.set_xlabel("time"); axs.set_ylabel("°C")
# axs.set_title("Line " * section * " Mid-Depth Temperature")
fig
fig.savefig(plotsdir("native/sections/PAC_MidDepth_Temperatures.png"), bbox_inches = "tight")


n_iterations = Int(1e5)
trends = zeros(n_iterations)
for k = 1:n_iterations
    t_rand = zeros(Int, 3)
    for (i, year) in enumerate([1994, 2004, 2013])
        random_day(year) = rand(tecco[findall( year .<= collect(tecco).<= year+1)])
        t_rand[i] = Base.findmin(abs.(collect(tecco) .- random_day(year)))[2]
    end
    tecco_sample = tecco[t_rand]
    E,F = trend_matrices(tecco_sample)

    θ_sample = θz[t_rand] .+ rand(Normal(0, mean(σ)), 3)

    a, b = F*θ_sample
    trends[k] = 1000 * b
end

fig, ax = plt.subplots()
ax.set_title("Distribution of Mid-Depth Pacific Temperature Trends")
sns.histplot(trends, ax = ax, stat = "probability")
tμ = mean(trends)
tσ = std(trends)
ax.axvline(tμ + (2 * tσ), c = "k", linestyle = "--")
ax.axvline(tμ - (2 * tσ), c = "k", linestyle = "--", label = "2σ")
pos_trend_est = sum(trends .> 0.0) / sum(trends .> -Inf)
pos_trend_est = round(pos_trend_est, digits = 2)
ax.scatter(0, nothing, label = "P(β > 0) = " * string(pos_trend_est), alpha = 0.0) 
ax.legend(frameon=false)
ax.set_xlabel("mK per year")
fig
fig.savefig(plotsdir("native/sections/" * section * "_MidDepth_Temp_Trend_Dist.png"), bbox_inches = "tight")
