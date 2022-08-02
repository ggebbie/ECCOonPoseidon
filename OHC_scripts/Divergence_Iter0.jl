
function smush(ma)
    temp = 0.0 .* similar(ma[:, 1])
    for ff in eachindex(ma)
        temp[ff[1]] .+= ma[ff]
    end
    return temp 
end
function levels_2_string(uplvl, botlvl)
    return string(abs(uplvl)) * " to " * string(abs(botlvl))
end
function level_mean(ma, depths)
    weight_ma = 
    var_sum_ma = 0.0 .* similar(ma[:, 1])
    depth_sum_ma = 0.0 .* similar(ma[:, 1])
    for ff in eachindex(weight_ma)
        for ij in eachindex(weight_ma[ff[1]])
            var_sum_ma[ff[1]][ij] = ma[ff][ij] .* depths[ff][ij]
            depth_sum_ma[ff[1]][ij] = depths[ff][ij]
        end
    end
    var_mean = var_sum_ma./depth_sum_ma
    return var_mean
end
function plot_level_mean!(var, λ, ϕ, mask, ax, color)
    projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
    ax.coastlines()
    ax.gridlines(crs=projPC, draw_labels=true,
                        linewidth=2, color="gray", alpha=0, linestyle="--")
    cf = Vector{Any}(undef ,1)
    # min_val = minimum(MeshArrays.mask(var ,Inf))
    max_val = maximum(MeshArrays.mask(var,-Inf))
    min_val = -max_val
    for ff in 1:length(var)
            cf[1] = ax.pcolormesh(λ[ff], ϕ[ff], var[ff],shading="nearest", 
            transform=projPC, rasterized = true, vmin = min_val, vmax = max_val,
            cmap = color)
    end
    return cf[1]
end
# depthlbl = string(abs(round(z[lvl_idx], digits = 2)))


# this script will compile and save a basinwide average of OHC 
# for a particular basin (basin_name) 

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
sns.set()
pygui(false)
cm = pyimport("cmocean.cm")
colorway = cm.balance;
"""Plotting Function """


include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
area = readarea(γ)
tecco = collect(Float64, 1992+1/24:1/12:2018)

runpath,diagpath = listexperiments(exprootdir())
shortnames = expnames()
marks = expsymbols()
nexp = length(shortnames) # number of experiments

ocean_mask = wet_pts(Γ)
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ)

uplvl = -2000; botlvl = -3000
tt3km = findall( botlvl .<= z[:].<= uplvl)
ignore_list = ["noIA", "129ff"]
msk = PAC_msk; region = "PAC"
cell_depths = get_cell_depths(msk, ΔzF, Γ.hFacC)
cell_volumes = get_cell_volumes(area, cell_depths);
""" comparing exp sfc temperatures """

println("Surface Temps")
filedir = "ECCO_vars/"
filename = "THETAs"

nz = 50
θ_sfc= Dict()
OHC_sfc = Dict()
lvls = 1
ρ = 1029; cₚ = 3850.0
for (keys,values) in shortnames
    if values ∉ ignore_list
        expname = keys
        @load datadir(filedir*filename*"_"*expname*".jld2") var_exp
        θ_sfc[keys] = [mean(var_exp[tt][:, lvls], cell_volumes[:, lvls]) for tt in 1:length(tecco)]
        # OHC_deep[keys] = θ_deep[keys] .* (sum(cell_volumes[:, lvls]) * ρ * cₚ * 1e-21)
    end
end

fig, ax = plt.subplots(1, 2, figsize = (12, 5))
ax[1].set_title("Divergence From 0bf (" * region * ", SST)")
ax[2].set_title(" Mean " *region * " SST")
plot_div_0bf!(θ_sfc, tecco, shortnames, ignore_list, ax[1], ylabel = "Δθ [ᵒK]")
plot_ts!(θ_sfc, tecco, shortnames, ignore_list, ax[2]; ylabel = "θ [ᵒK]")
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/" * "θ_" * region * "_sfc.png")

""" comparing exp deep temperatures """

println("2km -3km  Temps")
filedir = "ECCO_vars/"
filename = "THETAs"
nz = 50
θ_deep = Dict()
OHC_deep = Dict()
lvls = tt3km
ρ = 1029; cₚ = 3850.0ˇ
vol_sum = sum(cell_volumes[:, lvls])
for (keys,values) in shortnames
    if values ∉ ignore_list
        expname = keys
        @load datadir(filedir*filename*"_"*expname*".jld2") var_exp
        θ_deep[keys] = [mean(var_exp[tt][:, lvls], cell_volumes[:, lvls]) for tt in 1:length(tecco)]
        OHC_deep[keys] = θ_deep[keys] .* (vol_sum * ρ * cₚ * 1e-21)
    end
end
fig, ax = plt.subplots(1, 2, figsize = (12, 5))
ax[1].set_title("Divergence From 0bf (" * region * ", z=2-3km)")
ax[2].set_title(" Mean " * region * " θ, z=2-3km")
plot_div_0bf!(θ_deep, tecco, shortnames, ignore_list, ax[1]; ylabel = "Δθ [ᵒK]")
plot_ts!(θ_deep, tecco, shortnames, ignore_list, ax[2]; ylabel = "Δθ [ᵒK]")
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/" * "θ_" * region * "_2km3km.png")

println("2km -3km  OHC")
fig, ax = plt.subplots(1, 2, figsize = (12, 5))
ax[1].set_title("Divergence From 0bf (" * region * ", z=2-3km)")
ax[2].set_title(" Mean " * region * " H, z=2-3km")
plot_div_0bf!(OHC_deep, tecco, shortnames, ignore_list, ax[1]; ylabel = "ΔH [ZJ]")
plot_ts!(OHC_deep, tecco, shortnames, ignore_list, ax[2]; ylabel = "H [ZJ]")
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/" * "OHC_" *region * "_2km3km.png")

""" comparing exp surface air temp """

println("Air (2m) Temps")
filedir = "ECCO_vars/"
filename = "ATEMP"
nz = 1
θ_air = Dict()
for (keys,values) in shortnames
    if values ∉ ignore_list
        expname = keys
        @load datadir(filedir*filename*"_"*expname*".jld2") var_exp
        θ_air[keys] = [mean(var_exp[tt], area .* msk) for tt in 1:length(tecco)]
    end
end
fig, ax = plt.subplots(1, 2, figsize = (12, 5))
ax[1].set_title("Divergence From 0bf, " * region)
ax[2].set_title(" Mean " * region * " θ_air")
plot_div_0bf!(θ_air, tecco, shortnames, ignore_list, ax[1], ylabel =  L"Δθ_{air} [ᵒK]")
plot_ts!(θ_air, tecco, shortnames, ignore_list, ax[2], ylabel = "θ [ᵒK]")
fig.savefig(plotsdir() * "/OHC_Divergence/" * "θair_" * region * ".png")

""" comparing swdn radiation"""

println("SWDN")
filedir = "ECCO_vars/"
filename = "SWDN"
nz = 1
SW = Dict()
for (keys,values) in shortnames
    if values ∉ ignore_list
        expname = keys
        @load datadir(filedir*filename*"_"*expname*".jld2") var_exp
        SW[keys] = [sum(var_exp[tt] .* msk)/sum(msk) for tt in 1:length(tecco)]
    end
end
fig, ax = plt.subplots(1, 2, figsize = (12, 5))
ax[1].set_title("Divergence From 0bf, " * region)
ax[2].set_title(" Mean " * region * " SW Down")
plot_div_0bf!(SW, tecco, shortnames, ignore_list, ax[1], ylabel =  L"ΔSW [W/m^{2}]")
plot_ts!(SW, tecco, shortnames, ignore_list, ax[2], ylabel =  L"SW [W/m^{2}]")
fig.savefig(plotsdir() * "/OHC_Divergence/" * "SW" * region * ".png")

""" comparing lwdn radiation"""

println("LWDN")
filedir = "ECCO_vars/"
filename = "LWDN";
nz = 1
LW = Dict()
for (keys,values) in shortnames
    if values ∉ ignore_list
        expname = keys
        @load datadir(filedir*filename*"_"*expname*".jld2") var_exp
        LW[keys] = [sum(var_exp[tt] .* msk)/sum(msk) for tt in 1:length(tecco)]
    end
end
fig, ax = plt.subplots(1, 2, figsize = (12, 5))
ax[1].set_title("Divergence From 0bf, " * region)
ax[2].set_title(" Mean " * region * " LW Down")
plot_div_0bf!(LW, tecco, shortnames, ignore_list, ax[1], ylabel =  L"ΔLW [W/m^{2}]")
plot_ts!(LW, tecco, shortnames, ignore_list, ax[2], ylabel =  L"LW [W/m^{2}]")
fig.savefig(plotsdir() * "/OHC_Divergence/" * "LW" * region * ".png")
close()

"""Compare background diffusivities init. cond."""
diffs, keys = get_diff_arrs(γ)
avg_diff = Dict()
diff_tot = Dict()
for key in keys
    diff_tot[key] = diffs[key * "_adj"] .- diffs[key * "_unadj"]
    avg_diff[key] = mean(diff_tot[key][:, lvls], cell_volumes[:, lvls])
end
vert_diff = diff_tot["diffkr"]
proj0 = ECCOonPoseidon.cartopy.crs.Robinson(central_longitude=-150)
fig, ax = plt.subplots(1, 1, figsize=(6,6), subplot_kw=Dict("projection"=> proj0))
ax.set_extent((110, -70, -55, 60))
vert_diff_mean = level_mean(vert_diff[:, lvls], cell_depths[:, lvls])
lvls_str = levels_2_string(uplvl, botlvl)
ax.set_title("BG Vert. Diff. BG. Adjustment \n at " * lvls_str * " m ")
c1 = plot_level_mean!(vert_diff_mean, λ, ϕ, msk, ax, colorway)
cbar = fig.colorbar(c1, label = L"m^2 / s")
fig.savefig(plotsdir() * "/OHC Climatologies/" * "diffkr_" * region * ".png")


"""Compare vertical diffusivities of pot. temperature"""

println("Pot Flux")
filedir = "ECCO_vars/"
filenameA = "THETA_AFLX";
filenameE = "THETA_DFLXExplicit";
filenameI = "THETA_DFLXImplicit";
filenames = [filenameA, filenameI, filenameE]
θF = Dict(f => Dict() for f in filenames)
lvls = tt3km
where_lvl = (cell_volumes[:, lvls] .!= 0.0)
temp = 0.0.*similar(msk)
for (keys,values) in shortnames
    for filename in filenames
        if values ∉ ignore_list
            θF[filename][keys] = []
            expname = keys
            @load datadir(filedir*filename*"_"*expname*".jld2") var_exp
            for tt in 1:length(tecco)
                temp .= (var_exp[tt][:, lvls[1]] .-  var_exp[tt][:, lvls[end]+1]) .* msk #extract msk points
                temp .= temp ./ area
                push!(θF[filename][keys], sum(temp)) #total vertical heat transport in msk
            end
        end
    end
end
filedir = "ECCO_vars/"
filename = "THETAs"
nz = 50
θ_z = Dict()
lvls = tt3km
temp = 0.0.*similar(cell_depths[:, lvls])

for (keys,values) in shortnames
    if values ∉ ignore_list
        expname = keys
        @load datadir(filedir*filename*"_"*expname*".jld2") var_exp
        θ_z[keys] = [sum(smush(var_exp[tt][:, lvls] .* cell_depths[:, lvls])) for tt in 1:length(tecco)]
    end
end

tot_flux = Dict()
for (keys,values) in shortnames
    if values ∉ ignore_list
        tot_flux[keys] = zeros(length(tecco))
        for filename in filenames
                tot_flux[keys] .-= θF[filename][keys] 
        end
    end     
end



dθ_z = Dict()
dθ_z_approx = Dict()
resid = Dict()
means = Dict()
for (keys,values) in shortnames
    if values ∉ ignore_list
        dθ_z[keys] =  (θ_z[keys][2:end] - θ_z[keys][1:end-1]) / 2.62e6
        dθ_z_approx[keys] = (tot_flux[keys][1:end-1] .+ tot_flux[keys][2:end])./2
        resid[keys] = dθ_z[keys] .- dθ_z_approx[keys]
        means[keys] = mean(resid[keys])
    end
end

fig, ax = plt.subplots(2, 1, figsize = (12, 7))
ax[1].set_title("Divergence From 0bf (" * region * ", z=2-3km)")
ax[2].set_title(" Total " * region * " Vert. Heat Flux (LHS), z=2-3km")
plot_div_0bf!(dθ_z, tecco[2:end], shortnames, ignore_list, ax[1]; ylabel =  L"Δ ^\circ C m / s")
plot_ts!(dθ_z, tecco[2:end], shortnames, ignore_list, ax[2]; ylabel =  L" ^\circ C m / s", baselines = 0)
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/" * "TotθVFlxLHS_" * region * "_2km3km.png")

fig, ax = plt.subplots(2, 1, figsize = (12, 7))
ax[1].set_title("Divergence From 0bf (" * region * ", z=2-3km)")
ax[2].set_title(" Total " * region * " Vert. Heat Flux (RHS), z=2-3km")
plot_div_0bf!(dθ_z_approx, tecco[2:end], shortnames, ignore_list, ax[1]; ylabel =  L"Δ ^\circ C m / s")
plot_ts!(dθ_z_approx, tecco[2:end], shortnames, ignore_list, ax[2]; ylabel =  L" ^\circ C m / s", baselines = 0)
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/" * "TotθVFlxRHS_" * region * "_2km3km.png")

fig, ax = plt.subplots(1, figsize = (12, 7))
ax.set_title( region * " Vert. Heat Flux Residual (LHS - RHS), z=2-3km")
plot_ts!(resid, tecco[2:end], shortnames, ignore_list, ax; ylabel =  L" ^\circ C m / s", baselines = 0)
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/" * "TotθVFlxResid_" * region * "_2km3km.png")




