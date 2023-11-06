#this script uses T-S data to create T-S diagrams at a single point 

include("../../src/intro.jl")

using Revise,ECCOonPoseidon, ECCOtour, MeshArrays, 
MITgcmTools, JLD2, DrWatson, Statistics, JLD2, 
Printf, PyCall, LaTeXStrings
import NaNMath as nm, PyPlot as plt

include(srcdir("config_exp.jl"))'

(ϕ,λ) = latlonC(γ);
area = readarea(γ);
meshgrid(x, y) = (x' .* ones(length(y)), ones(length(x))' .* y);

region = "Sargasso"
tecco= 1992+1/24:1/12:2018 # ecco years
nt = length(tecco)

#load data 
savename = datadir("native/iter129_bulkformula" * region * "_WaterProp_profile_budget_z.jld2")
θS = load(savename)["θS"]
#specify maximum depth 
kmax = 30; z[kmax]

#functions that takes averages in 2-year chunks
yearly_avg(x) = hcat([mean(x[:, (i):(i+ 23)], dims = 2) for i in 1:24:312]...)

#extract data
θ = θS["θ"][1:kmax, :]; S = θS["S"][1:kmax, :]

#compute density
σ0 = densityJMD95.(θ,S, NaN, 0) .- 1000. #EOS from MITGCM 

θ = yearly_avg(θ); S = yearly_avg(S); σ0 = yearly_avg(σ0)
tecco_avg = yearly_avg(tecco')[:]; nt_avg = length(tecco_avg)

S_min, S_max = extrema(S); θ_min, θ_max = extrema(θ)
S_levels = S_min:0.005:S_max; θ_levels = θ_min:0.1:θ_max;
S_grid, θ_grid = meshgrid(S_levels, θ_levels);
σgrid = densityJMD95.(θ_grid,S_grid, 0.0, p₀) .- 1000.; #EOS from MITGCM does not require ref pressure
σ_min, σ_max = extrema(σgrid); levels = σ_min:0.3:σ_max
ns = floor(Int, size(S_grid, 2) /2)

#plot_ts diagram
fig, ax = plt.subplots(1,1, figsize = (10, 10));
cmap = "Blues"
norm = plt.matplotlib.colors.Normalize(vmin=1992, vmax=2017)

for (i, tt) in enumerate(1:nt_avg)
    ax.scatter(S[:, tt], θ[:, tt], c = repeat([tecco_avg[i]], 30), 
    cmap=cmap, vmin = 1992, vmax = 2017, edgecolors="black");
end

CS = ax.contour(S_grid[:, ns:end], θ_grid[:, ns:end], σgrid[:, ns:end], 
levels = levels, colors = "black", linewidths = 0.75);
CS = ax.contour(S_grid[:, 1:ns], θ_grid[:, 1:ns], σgrid[:, 1:ns], 
levels = levels, colors = "black", linewidths = 0.75);
ax.clabel(CS, fontsize=15, inline=true, fmt = "%.2f");
ax.set_xlabel("Practical Salinity"); ax.set_ylabel("Potential Temperature");

ax.set_title("65W, 30N; Upper 1000 meters")
fig.colorbar(plt.matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=ax, orientation="horizontal", label="time", fraction = 0.04, extend = "both")
fig.savefig(plotsdir("native/subtropical_gyre/TS_Upper1000.png"))
fig
