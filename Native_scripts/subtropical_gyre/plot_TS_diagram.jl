#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour, MeshArrays, 
MITgcmTools, JLD2, DrWatson, Statistics, JLD2, 
Printf, PyCall, LaTeXStrings
import NaNMath as nm, PyPlot as plt
@pyimport cmocean.cm as cmo
@pyimport seaborn as sns;
sns.set_theme(context = "poster", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
include(srcdir("config_exp.jl"))'

(ϕ,λ) = latlonC(γ);
area = readarea(γ);
meshgrid(x, y) = (x' .* ones(length(y)), ones(length(x))' .* y);

region = "Sargasso"
tecco= 1992+1/24:1/12:2018 # ecco years
nt = length(tecco)

# vm = maximum(abs.(θ["iter129_bulkformula"] - mean(θ["iter129_bulkformula"], dims = 2)))
savename = datadir("native/iter129_bulkformula" * region * "_WaterProp_budget_z.jld2")
θS = load(savename)["θS"]
kmax = 40; z[kmax]
yearly_avg(x) = hcat([mean(x[:, (i):(i+ 23)], dims = 2) for i in 1:24:312]...)

θ = θS["θ"][1:kmax, :]; S = θS["S"][1:kmax, :]
σ0 = OHC_helper.densityJMD95.(θ,S, NaN, 0) .- 1000. #EOS from MITGCM 
θ = yearly_avg(θ)
S = yearly_avg(S)
σ0 = yearly_avg(σ0)
tecco_avg = yearly_avg(tecco')'; nt_avg = length(tecco_avg)

fig, ax = plt.subplots(1,1, figsize = (10, 10));

for (i, tt) in enumerate(1:nt_avg)
    ax.scatter(S[:, tt], θ[:, tt], c = "blue", alpha = i* (1/nt_avg), edgecolors="black");
end
fig
S_min, S_max = extrema(S)
θ_min, θ_max = extrema(θ)


S_levels = S_min:0.005:S_max; θ_levels = θ_min:0.1:θ_max;
S_grid, θ_grid = meshgrid(S_levels, θ_levels);
σgrid = OHC_helper.densityJMD95.(θ_grid,S_grid, 0.0, p₀) .- 1000.; #EOS from MITGCM does not require ref pressure
σ_min, σ_max = extrema(σgrid); levels = σ_min:0.3:σ_max

CS = ax.contour(S_grid, θ_grid, σgrid, levels = levels, colors = "black", linewidths = 0.75);
ax.clabel(CS, fontsize=15, inline=true, fmt = "%.2f");
ax.set_xlabel("Practical Salinity"); ax.set_ylabel("Potential Temperature");
fig