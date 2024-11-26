using GibbsSeaWater

uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)

fig,axs=plt.subplots(1,1, sharey = true, figsize = (10, 7.5))
θ = load(fnames(section, expname))["θ"][lvls, :]
# SP = load(fnames(section, expname))["S"][lvls, :]
# P = zero(SP); [P[:, i] = abs.(pstdz[lvls]) for i = 1:312]
# SA = gsw_sa_from_sp.(SP, P, zero(SP), zero(SP))
# θ = gsw_ct_from_pt.(SA, θP)

ΔV = load(fnames(section, expname))["ΔV"][lvls]
volume_weight_column(x, ΔV) =  sum(x .* ΔV, dims =1) ./ sum(ΔV)
θz = volume_weight_column(θ, ΔV)[:]
E,F = trend_matrices(tecco)
a, b = F*θz
fit_line(x, a, b) = b .* (x .- mean(x)) .+ a
y0b = fit_line(tecco, a, b)
σ = std(θz .- y0b)
axs.plot(tecco, θz, label = "ECCO Output");
axs.plot(tecco, y0b, label = "Trend = " * string(round(1000 * b, digits = 2)) * " m°C per year");
axs.fill_between(tecco, θz .- (2 .*σ), θz .+ (2 .*σ), alpha=0.3, color="k", label = "2σ estimate")
axs.legend(frameon=false)
axs.set_xlabel("time"); axs.set_ylabel("°C")
axs.set_title("Line " * section * " Mid-Depth Temperature")
fig
fig.savefig(plotsdir("native/sections/" * section * "_MidDepth_Temperatures.png"), bbox_inches = "tight")

n_iterations = Int(1e6)
trends = zeros(n_iterations)
for k = 1:n_iterations
    t_rand = zeros(Int, 3)
    for (i, year) in enumerate([1994, 2004, 2013])
        random_day(year) = rand(tecco[findall( year .<= collect(tecco).<= year+1)])
        t_rand[i] = findmin(abs.(collect(tecco) .- random_day(year)))[2]
    end
    tecco_sample = tecco[t_rand]
    E,F = trend_matrices(tecco_sample)

    θ_sample = θz[t_rand] .+ rand(Normal(0, σ), 3)

    a, b = F*θ_sample
    trends[k] = 1000 * b
end

fig, ax = plt.subplots()
ax.set_title("Distribution of Mid-Depth " *section *"Temperature Trends")
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
