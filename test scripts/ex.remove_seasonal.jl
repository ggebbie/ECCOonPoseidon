include("../src/intro.jl")
using Revise 
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools,
    PyPlot, JLD2, DrWatson, Statistics, JLD2
import CairoMakie as Mkie
import GeoMakie
pygui(true)

fcycle = 1 # units: yr^{-1}
# for removing seasonal cycle from monthly averaged timeseries
# overtones= 6; # semi-annual, quad-annual, etc.: helps capture an asymmetric seasonal cycle
overtones= [4];

fig, axs = subplots(3, 1, figsize=(5,5))
ax = axs[1]
x = collect(0:0.05:30)
seasoned = sin.(2 * π .* (x)) 
for ov in overtones
    Ecycle,Fcycle = seasonal_matrices(fcycle,x,ov)
    deseasoned = remove_seasonal(seasoned,Ecycle,Fcycle)
    ax.plot(x,deseasoned, linewidth=1, markersize = 2, label = string(ov))
end
ax.plot(x,seasoned, linewidth=1, markersize = 2, label = "original")
# ax.hlines(mean(seasoned), xmin = 0, xmax = 30, c = "black")
ax.legend()
ax.set_xlabel("Years")
ax.set_ylabel("Value")
ax.set_title("Overtones Removed from " *  L"sin(2 \pi x) ")

ax = axs[2]
x = collect(0:0.05:30)
seasoned = sin.(2 * π .* (x)) .+ x  
for ov in overtones
    Ecycle,Fcycle = seasonal_matrices(fcycle,x,ov)
    deseasoned = remove_seasonal(seasoned,Ecycle,Fcycle)
    ax.plot(x,deseasoned, linewidth=1, markersize = 2, label = string(ov))
end

ax.plot(x,seasoned, linewidth=1, markersize = 2, label = "original")
# ax.hlines(mean(seasoned), xmin = 0, xmax = 30, c = "black")
ax.legend()
ax.set_xlabel("Years")
ax.set_ylabel("Value")
ax.set_title("Overtones Removed from " *  L"sin(2 \pi x) + x")
fig.tight_layout()


ax = axs[3]
x = collect(0:0.05:30)
seasoned = sin.(2 * π .* (x)) .+ sin.(5 * π .* (x)) .+ x  
for ov in overtones
    Ecycle,Fcycle = seasonal_matrices(fcycle,x,ov)
    deseasoned = remove_seasonal(seasoned,Ecycle,Fcycle)
    ax.plot(x,deseasoned, linewidth=1, markersize = 2, label = string(ov))
end
ax.plot(x,seasoned, linewidth=1, markersize = 2, label = "original")
# ax.hlines(mean(seasoned), xmin = 0, xmax = 30, c = "black")
ax.legend()
ax.set_xlabel("Years")
ax.set_ylabel("Value")
ax.set_title("Overtones Removed from " *  L"sin(2 \pi x) + sin(5 \pi x) + x")
fig.tight_layout()

