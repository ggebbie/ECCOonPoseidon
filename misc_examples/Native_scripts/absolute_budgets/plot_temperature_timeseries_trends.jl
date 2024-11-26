#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, FLoops
using .OHC_helper
import PyPlot as plt 
using PyCall
include(srcdir("config_exp.jl"))
cmo = pyimport("cmocean.cm");

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
@pyimport seaborn as sns;

sns.set_theme(context = "notebook", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
              
runpath,diagpath = listexperiments(exprootdir());
diagpath["seasonalclimatology_iter0_multarr0"] = "/batou/ECCOv4r4/exps/seasonalclimatology_iter0/run_multarr0/diags/"

# abbreviations for each experiment for labels, etc.
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
nexp = length(shortnames) # number of experiments
 
ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));
sum(cell_volumes)


H = OHC_helper.smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
tecco = 1992+1/24:1/12:2018
nz = 50

ΔV = zeros(Float32, nz)
[ΔV[k] = Float32(sum(cell_volumes[:, k])) for k=1:nz]
lvls = findall( -3000 .<= z[:].<= -2000)

masked_volume = cell_volumes[:, lvls]
sum_masked_volume = Float32(sum(masked_volume))

fname = datadir("native/" * region * "_THETA_levels" * ".jld2")
θ = load(fname)["θ"]
titles = ["ECCO V4r4 with Full Forcing", "ECCO V4r4 with Climatological Forcing"]
E,F = trend_matrices(tecco)

#seconds in a year = 3.154e+7
#Area of earth = 5.1 * 10^14
fig, ax = plt.subplots(figsize = (5, 7.5))
sumsquares(x, y, dims) = sum((x .- y).^2, dims = dims)
# sqrt(sum((tecco .- mean(tecco)).^2, dims = 1))
# inv(sumsquares(tecco, mean(tecco), 1))
for (i, expname) in enumerate(["iter129_bulkformula"])
    fit = F * θ[expname]'
    θhat = (F * θ[expname]' )'* E' 
    
    se = sqrt.(sumsquares(θ[expname], θhat, 2) .* inv(sumsquares(tecco, mean(tecco), 1)[1]) * inv(310)) 
    se = vec(2 .*se)

    ax.plot(fit[2, :] * 1e17 * 3.154e-7* 5e-14 * 4e6 * 100, abs.(z))
    err1 = (vec(fit[2, :]).+se)  .* (1e17 * 3.154e-7 * 5e-14 * 4e6 * 100 )
    err2 = (vec(fit[2, :]).-se)  .* (1e17 * 3.154e-7* 5e-14 * 4e6 * 100)

    ax.fill_betweenx(y=abs.(z), x1 =err2, x2 =err1, color = "orange")
    ax.set_title(titles[i])
    ax.set_ylim(1000, 4000)

    ax.set_xlim(-2000, 100)
    ax.set_xlabel("W per m²"); ax.set_ylabel("Depth [m]")
end
ax.invert_yaxis()

fig.suptitle("North Pacific Temperature Anomaly [cK]")
fig
# fig.savefig(plotsdir("native/TempAnomalies" * region * ".png"), dpi = 500)
86400 * 365
sum(ΔV[lvls]) 
1.58e17/ sum(cell_volumes[:, lvls])