using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions
using PyCall
using DataFrames, GLM
using Plots
# import Plots: Animation, frame, gif, cgrad
# import Plots: contourf as jcontourf
# import Plots: plot as jplot
using ColorSchemes
include(srcdir("config_exp.jl"))


a = Animation()
Y = reverse(z[lvls])
X = ma_zonal_sum(ϕ .* area) ./ ma_zonal_sum(area)
ϕ_d = similar(crop_vols)
expname1 = "iter0_bulkformula"
expname2 = "iter129_bulkformula"
exps = [expname1, expname2]
c1 = minimum([nm.minimum(θ_zonal[ex]) for ex ∈ exps])
c2 = maximum([nm.maximum(θ_zonal[ex]) for ex ∈ exps])
tt = 1
for tt in 1:length(tecco)
    ps = []
    for ex in exps
        horiz_avg = @views reverse(θ_zonal[ex][:, :, tt], dims = 1)
        jcf = Plots.contourf(X, Y, horiz_avg,
        xlabel = "latitude [º]",
        xticks = -70:20:70,
        ylabel = "depth [m]", 
        linewidth = 0.1,
        clim = (c1, c2),
        c = :balance, 
        colorbar_title= "°C",
        title =ex, 
        titlefontsize = 12)
        push!(ps, jcf)
    end
    l = @layout [a;b]
    p = Plots.plot(ps[1], ps[2], size=(500,700), 
    margin = 4mm, link = :both,
    plot_title =  "Zonally averaged θ"* region * "\n Time: " * string(round(tecco[tt], digits = 6)),
    plot_titlefontsize	= 13, plot_titlevspan = 0.075, layout = l)
    frame(a, p)
end
gif(a, plotsdir() * "/OHC_Divergence/" * "θ_zonal_avg_" * region *"_.gif", fps = 10)

# fig.savefig(plotsdir() * "/OHC_Divergence/" * "θDepthTS_" * region * suffix * ".png")

