include("../src/intro.jl")
include("../src/OHC_helper.jl")
using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools,
PyPlot, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf
using .OHC_helper
using PyPlot   # important!
using PyCall
@pyimport seaborn as sns
sns.set(); pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;
marks = expsymbols()
nexp = length(shortnames) # number of experiments

dθz_LHS,dθz_RHS, resids, means = Dict(), Dict(), Dict(), Dict()
Δt₁ = Float32.((tecco[3:end] .- tecco[1:end-2]) .* (12 * 30.5 * 86400))
central_diff(x, dt) = (x[3:end] .- x[1:end-2]) ./ (dt)
Δt₂ = mean(Float32.((tecco[2:end] .- tecco[1:end-1]) .* (12 * 30.5 * 86400)))
fwd_diff(x, dt) = (x[2:end] .- x[1:end-1]) ./ (dt)
fwd_mean(x)=(x[2:end].+x[1:end-1])./2
perc_diff(x, y) = 100 * (x - y)/y
resid(x, y) = x - y
fcycle = 1 # units: yr^{-1}
# for removing seasonal cycle from monthly averaged timeseries
overtones= 2; # semi-annual, quad-annual, etc.: helps capture an asymmetric seasonal cycle
Ecycle,Fcycle = seasonal_matrices(fcycle,tecco[2:end],overtones)

metric_name = "PercDiffPruned"
# metric_units = L"Δ ^\circ C / s"
metric_units = "% Difference"
remove_seasonal_with_mean(x, Ecycle, Fcycle) = remove_seasonal(x,Ecycle,Fcycle) .+ mean(x)
for (keys,values) in shortnames
    dθz_LHS[keys] =  Float32.(fwd_diff(θz[keys], Δt₂))
    dθz_RHS[keys] =  Float32.(fwd_mean(θC[keys]))
    dθz_LHS[keys] = remove_seasonal_with_mean(dθz_LHS[keys],Ecycle,Fcycle)
    dθz_RHS[keys] = remove_seasonal_with_mean(dθz_RHS[keys],Ecycle,Fcycle)
    # resids[keys] = resid.(dθz_LHS[keys],dθz_RHS[keys])
    resids[keys] = prune(perc_diff.(dθz_RHS[keys],dθz_LHS[keys]))
    σ=std(prune(resid.(dθz_LHS[keys],dθz_RHS[keys])))/std(dθz_LHS[keys])
    println(values," ratio of STD ",string(σ))
end

tecco_off = tecco[1:end-1] .+  (diff(tecco) ./2)

fig, ax = plt.subplots(1, figsize = (12, 7))
ax.set_title( region * " Heat Flux " * metric_name *", z=2-3km")
plot_ts!(resids, tecco_off, shortnames, ignore_list, ax; ylabel =  metric_units, baseline =  0)
tight_layout(); savefig(plotsdir() * "/OHC_Divergence/" * "TotθFlx" * metric_name * "_" * region * "_2km3km.png")

slice(i) =pycall(pybuiltin("slice"), PyObject, i)
slice(i, j) =pycall(pybuiltin("slice"), PyObject, i, j)
fig = plt.figure( figsize = (12, 7))
gs = fig.add_gridspec(4,5)
element(i,j) = get(gs, (i,j))
a = fig.add_subplot(element(slice(0,4), slice(2,5)))
a.boxplot(collect(values(resids)), showmeans=true)
a.set_xticklabels(keys(resids))
a.set_title(region * " Heat Flux Budget " * metric_name * ", z=2-3km")
a.set_ylabel(metric_units)

axs = [fig.add_subplot(element(0,slice(0,2))), fig.add_subplot(element(1,slice(0,2))), 
fig.add_subplot(element(2,slice(0,2))), fig.add_subplot(element(3,slice(0,2)))]
i = 0 
for (keys,values) in shortnames
    i+=1
    axs[i].hist(resids[keys], bins = 25, alpha = 0.3)
    axs[i].set_ylabel(L"Frequency")
    axs[i].set_xlabel(metric_units)
    axs[i].set_title(keys)
end

