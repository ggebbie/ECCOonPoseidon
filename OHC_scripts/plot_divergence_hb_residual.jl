using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf,
RollingFunctions
using PyCall
cm = pyimport("cmocean.cm");colorway = cm.balance;
marks = expsymbols()
nexp = length(shortnames) # number of experiments

dθz_LHS,dθz_RHS, resids, means = Dict(), Dict(), Dict(), Dict()
dθz_LHS_Annual, dθz_RHS_Annual, resids_Annual = Dict(), Dict(), Dict()
Δt₁ = mean(((tecco[3:end] .- tecco[1:end-2]) .* (12 * 30.5 * 86400)))
Δt₂ = mean(((tecco[2:end] .- tecco[1:end-1]) .* (12 * 30.5 * 86400)))

fcycle = 1 # units: yr^{-1}
# for removing seasonal cycle from monthly averaged timeseries
overtones= 2; # semi-annual, quad-annual, etc.: helps capture an asymmetric seasonal cycle
Ecycle,Fcycle = seasonal_matrices(fcycle,tecco[2:end],overtones)

metric_name = "PercDiffPruned"
# metric_units = L"Δ ^\circ C / s"
metric_units = "% Difference"
remove_seasonal_with_mean(x, Ecycle, Fcycle) = remove_seasonal(x,Ecycle,Fcycle) .+ mean(x)
window_size=25
for (keys,values) in shortnames
    dθz_LHS[keys] =  fwd_diff(θz[keys], Δt₂)
    dθz_RHS[keys] =  fwd_mean(θC[keys])
    dθz_LHS_Annual[keys] = rollmean(dθz_LHS[keys], window_size)
    dθz_RHS_Annual[keys] = rollmean(dθz_RHS[keys], window_size)

    error = resid.(dθz_LHS[keys],dθz_RHS[keys])
    resids[keys] = prune(perc_diff.(dθz_RHS[keys],dθz_LHS[keys]))
    resids_Annual[keys] = prune(perc_diff.(dθz_RHS_Annual[keys],dθz_LHS_Annual[keys]))
    error_annual = resid.(dθz_LHS_Annual[keys],dθz_RHS_Annual[keys])

    σ=std(prune(error))/std(dθz_LHS[keys])
    println(values," ratio of STD ",string(σ))
    
    println(values," avg error ", mean(prune(error)))
    σ_annual=std(prune(error_annual))/std(dθz_LHS_Annual[keys])
    println(values," ratio of STD Annual ",string(σ_annual))
    println(values," avg error annual", mean(prune(error_annual)))
end

tecco_off = tecco[1:end-1] .+  (diff(tecco) ./2)
# tecco_off = tecco[2:end-1]
fig, ax = plt.subplots(1, figsize = (12, 7))
ax.set_title( region * " Heat Flux " * metric_name *"," * suffix)
plot_ts!(resids, tecco_off, shortnames, ignore_list, ax; 
ylabel =  metric_units, baseline =  0, colors = colors)
sns.move_legend(ax, "lower center", bbox_to_anchor=(.5, 0), ncol=4)
fig.tight_layout();fig.savefig(plotsdir() * "/OHC_Divergence/Stats/" * "TotθFlx" * metric_name * "_" 
* region * "_"  * suffix * ".png")

fig, ax = plt.subplots(1, figsize = (12, 7))
ax.set_title(region * suffix * " Heat Flux Budget " * metric_name * ", " * suffix)
ax.set_ylabel(metric_units)
sns.violinplot(data = pd.DataFrame(resids), ax = ax)
fig.tight_layout();fig.savefig(plotsdir() * "/OHC_Divergence/Stats/" * "TotθFlx" * metric_name * "Violin_" 
* region * "_" * suffix *  ".png")

fig = plt.figure( figsize = (12, 7))
gs = fig.add_gridspec(4,5)
element(i, j) = element_gs(i, j, gs)
a = fig.add_subplot(element(slice(0,4), slice(2,5)))
a.boxplot(collect(values(resids)), showmeans=true); a.set_xticklabels(collect(keys(resids)))
a.set_title(region * suffix * " Heat Flux Budget " * metric_name *", " * suffix)
a.set_ylabel(metric_units)

axs = [fig.add_subplot(element(0,slice(0,2))), fig.add_subplot(element(1,slice(0,2))), 
fig.add_subplot(element(2,slice(0,2))), fig.add_subplot(element(3,slice(0,2)))]
for (i, key) in enumerate(keys(shortnames))
    axs[i].hist(resids[key], bins = 25, alpha = 0.3, color = colors[i])
    axs[i].set_ylabel(L"Frequency")
    axs[i].set_xlabel(metric_units)
    axs[i].set_title(key)
end
fig.tight_layout();fig.savefig(plotsdir() * "/OHC_Divergence/Stats/" * "TotθFlx" * metric_name * "Stats_" 
* region * "_" * suffix *  ".png")

fig = plt.figure( figsize = (12, 7))
gs = fig.add_gridspec(4,5)
element(i, j) = element_gs(i, j, gs)
a = fig.add_subplot(element(slice(0,4), slice(2,5)))
a.boxplot(collect(values(resids_Annual)), showmeans=true); a.set_xticklabels(keys(resids))
a.set_title(region * "Annual Avg. Heat Flux Budget " * metric_name *", " * suffix)
a.set_ylabel(metric_units)
axs = [fig.add_subplot(element(0,slice(0,2))), fig.add_subplot(element(1,slice(0,2))), 
fig.add_subplot(element(2,slice(0,2))), fig.add_subplot(element(3,slice(0,2)))]
for (i, key) in enumerate(keys(shortnames))
    axs[i].hist(resids_Annual[key], bins = 25, alpha = 0.3)
    axs[i].set_ylabel(L"Frequency")
    axs[i].set_xlabel(metric_units)
    axs[i].set_title(key)
end

fig.tight_layout();fig.savefig(plotsdir() * "/OHC_Divergence/Stats/" * "TotθFlx" * metric_name * 
"StatsAnnual_" * "_" * region *  "_" * suffix * ".png")
