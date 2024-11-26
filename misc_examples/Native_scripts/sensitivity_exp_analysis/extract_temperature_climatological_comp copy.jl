include("../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
import PyPlot as plt

include(srcdir("plot_and_dir_config.jl"))
@pyimport matplotlib.patches as patches

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)

cell_depths = get_cell_thickness(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
ΔV = lateral_sum(cell_volumes)

lvls = findall( -3000 .<= -z[:].<= -2000)
mid_depths(x) = vec(sum(Float32, x[lvls, :] .* ΔV[lvls], dims = 1) / sum(Float32, ΔV[lvls]))

tecco = 1992+1/24:1/12:2018;
get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )
function get_temperature(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid, cell_volumes)

    datafilelist_θ  = get_datafiles(expname, "state_3d_set1")[end-36:end-24]
    nt = length(datafilelist_θ);
    println(nt, " months available")
    ma_template = MeshArray(γ,Float32,50)
    θ_mean = MeshArray(γ,Float32,50); fill!(θ_mean, 0.0)

    @time for tt = 1:nt
        println(tt)
        fnameθ = datafilelist_θ[tt]
        @time θ = γ.read(diagpath[expname]*fnameθ,ma_template)

        θ_mean .+= θ ./ nt
    end

    # fnameθ = datafilelist_θ[nt]
    # θ_mean = γ.read(diagpath[expname]*fnameθ,ma_template)

    return θ_mean
end

mid_depths(x, cell_thickness, lvls) = vertical_sum(x[:, lvls] .* cell_thickness[:, lvls]) ./ vertical_sum(cell_thickness[:, lvls])
filter_zeros!(x) = [x.f[a][x.f[a] .== 0.0] .= NaN for a in eachindex(x)]
filter_Inf!(x) = [x.f[a][(!isfinite).(x.f[a])] .= NaN for a in eachindex(x)]
lvls
adjust_exps = Dict()

adjust_exps["iter0_bulkformula"] = get_temperature(diagpath, "iter0_bulkformula", γ, cell_volumes)
adjust_exps["iter129_bulkformula"] = get_temperature(diagpath, "iter129_bulkformula", γ, cell_volumes)

include(srcdir("MeshArraysPlots.jl"))

reg_λ = LLCcropC360(λ, γ; modify_λ = true); reg_ϕ = LLCcropC360(ϕ, γ)

fig, axs = plt.subplots(1, 3, figsize=(20,10), subplot_kw=Dict("projection"=> proj0))
vmin = 1.4; vmax = 1.7; levels = vmin:0.02:vmax
[a.set_extent((120, 285, -40, 57),crs=projPC) for a in axs]
[a.coastlines() for a in axs]
data = mid_depths(adjust_exps["iter0_bulkformula"], cell_depths, lvls)
# data = 1 .* adjust_exps["iter0_bulkformula"][:,  42]
filter_zeros!(data); filter_Inf!(data)
axs[1].contourf(reg_λ, reg_ϕ, LLCcropC360(data, γ), vmin = vmin, vmax = vmax, 
cmap = cmo.thermal, transform=projPC, 
levels = levels, extend = "both")
axs[1].set_title("CTRL (Iteration 0) \n 2015")
fig
data = mid_depths(adjust_exps["iter129_bulkformula"], cell_depths, lvls)
# data = adjust_exps["iter129_bulkformula"][:,  42] .*1
filter_zeros!(data); filter_Inf!(data)
cb = axs[2].contourf(reg_λ, reg_ϕ, LLCcropC360(data, γ), vmin = vmin, vmax = vmax, 
cmap = cmo.thermal, transform=projPC, 
levels = levels, extend = "both")
axs[2].set_title("FULL (Iteration 129) \n 2015")

data = mid_depths(adjust_exps["iter129_bulkformula"] .- adjust_exps["iter0_bulkformula"], cell_depths, lvls)
vmin = -0.025; vmax = 0.025; levels = vmin:0.005:vmax
filter_zeros!(data); filter_Inf!(data)
cbd = axs[3].contourf(reg_λ, reg_ϕ, LLCcropC360(data, γ), vmin = vmin, vmax = vmax, 
cmap = cmo.balance, transform=projPC, 
levels = levels, extend = "both")
axs[3].set_title("Difference")

fig.colorbar(cb, ax = axs[1:2], orientation = "horizontal", fraction = 0.05, label = "[°C]"); 

fig.colorbar(cbd, ax = axs[3], orientation = "horizontal", fraction = 0.05, label = "[K]"); 

fig.savefig(plotsdir("native/sensitivity_exps/FinalTemps.png"))
fig