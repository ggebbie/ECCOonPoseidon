using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools,
PyPlot, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions
using PyPlot   # important!
using PyCall
using DataFrames, GLM
import Plots: Animation, frame, gif, cgrad
import Plots: contourf as jcontourf

a = Animation()
Y = reverse(z[lvls])
X = ma_zonal_sum(ϕ .* area) ./ ma_zonal_sum(area)
ϕ_d = similar(crop_vols)
# for ff in eachindex(ϕ_d)
#     ϕ_d[ff] .= ϕ[ff[1]]
# end
# X = ma_zonal_avg(ϕ_d, crop_vols)[1, :]
c1 = nanminimum(θ_zonal["iter129_bulkformula"])
c2 = nanmaximum(θ_zonal["iter129_bulkformula"])
for tt in 1:length(tecco)
    zonal_view = @views reverse(θ_zonal["iter129_bulkformula"][:, :, tt], dims = 1)
    jcf = jcontourf(X, Y, zonal_view,
    xlabel = "latitude [º]",
    ylabel = "depth [m]", 
    linewidth = 0.1,
    clim = (c1, c2),
    title = "Zonally averaged θ "* region * 
    "\n Time: " * string(round(tecco[tt], digits = 6)))
    # c = cgrad(:thermal))
    frame(a, jcf)
end
gif(a, plotsdir() * "/OHC_Divergence/" * "first_example.gif", fps = 10)

# fig.savefig(plotsdir() * "/OHC_Divergence/" * "θDepthTS_" * region * suffix * ".png")

