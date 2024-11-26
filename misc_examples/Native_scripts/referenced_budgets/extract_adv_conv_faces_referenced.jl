#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall

    import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");
@pyimport seaborn as sns;

colors =  sns.color_palette("deep")[1:4]

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
ignore_list= ["noIA", "129ff"]
shortnames = reduce_dict(expnames(), ignore_list)
tecco = 1992+1/24:1/12:2018

ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region, extent = "not", include_bering = false)
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= z[:].<= uplvl)

abs_dist(x, r) = abs(x) < r

#only works for NPAC mask 
ϕ_mask_min = Float32(OHC_helper.get_min_lat(ϕ, PAC_msk)); ϕ_min_mask = ϕ .> -Inf
region2 = "PAC56"
PAC56_mask = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region2, extent = "not", include_bering = false)
[ϕ_min_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- 23.8, 0.1) .* PAC56_mask.f[ff] for ff in 1:2]
[ϕ_min_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- 22.9, 0.1) .* PAC56_mask.f[ff] for ff in 4:5]

ϕ_mask_max = Float32(OHC_helper.get_max_lat(ϕ, PAC_msk)); ϕ_max_mask = ϕ .> -Inf; 
[ϕ_max_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- 51.0, 0.1) .* PAC56_mask.f[ff] for ff in 1:2]
[ϕ_max_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- 50.3, 0.1) .* PAC56_mask.f[ff] for ff in 4:5]

cs, sn = OHC_helper.get_cs_and_sn(γ)

rotate_uv
transports_dict = Dict()
runpath,diagpath = listexperiments(exprootdir());
vars = ["iter129_bulkformula"]

cell_volume = sum(cell_volumes[:, lvls])
function filter_advection_budget(diagpath::Dict{String, String}, expname::String, γ::gcmgrid, 
    cs::MeshArrays.gcmarray{T,1,Matrix{T}}, sn::MeshArrays.gcmarray{T,1,Matrix{T}}) where T<:AbstractFloat
    filelist = searchdir(diagpath[expname],"trsp_3d_set2")
    datafilelist_H  = filter(x -> occursin("data",x),filelist) 
    filelist = searchdir(diagpath[expname],"trsp_3d_set3")
    datafilelist_R  = filter(x -> occursin("data",x),filelist)

    Wtops = Float32[]; Wbots = Float32[]; Vsouths = Float32[]; Vnorths = Float32[]
    Vsouthstrue = Float32[]
    nt = 36;
    uv_conv = MeshArray(γ,Float32,length(lvls))
    @time for tt = 1:nt
        println(tt)
        κUθ, κVθ, Uθ, Vθ = OHC_helper.extract_heatbudgetH(diagpath, expname , datafilelist_H[tt], γ)
        κzθ, wθ = OHC_helper.extract_heatbudgetR(diagpath, expname , datafilelist_R[tt], γ)

        OHC_helper.calc_UV_conv3D!(Uθ[:, lvls], Vθ[:, lvls], uv_conv)
        Etrsp, Ntrsp = rotate_UV_native(Uθ[:, lvls], Vθ[:, lvls], cs, sn) 
        Wtop = wθ[:, lvls[1]]
        Wbot = wθ[:, lvls[end] + 1]

        push!(Wtops,  sum(Wtop .* PAC_msk ) / cell_volume)
        push!(Wbots,  sum(Wbot .* PAC_msk ) / cell_volume)
        push!(Vsouths,  sum(Ntrsp .* ϕ_mask ) / cell_volume) #according to the rotations, 
        push!(Vnorths,  sum(Ntrsp .* ϕ_max_mask ) / cell_volume) #according to the rotations, 
        push!(Vsouthstrue,  sum(uv_conv .* PAC_msk ) / cell_volume)
    end
    return (Wtop = deepcopy(Wtops), Wbot = deepcopy(Wbots), 
            Vin = deepcopy(Vsouths), Vout = deepcopy(Vnorths), 
            Vinout = deepcopy(Vsouthstrue))
end

for expname in (vars)
    transports_dict[expname] = filter_volume_budget(diagpath, expname, γ, cs, sn)
end

i = 1 ;var = vars[i]
withoutBering = 1
Vout = withoutBering .* transports_dict[var].Vout
conv = (transports_dict[var].Wbot .- transports_dict[var].Wtop).+ transports_dict[var].Vin .- Vout
conv_true = (transports_dict[var].Wbot .- transports_dict[var].Wtop).+ transports_dict[var].Vinout 
nt = length(conv_true)
fig, axes = plt.subplots(1, 4, figsize = (15, 4), sharey = true)
axes[1].plot(tecco[1:nt], (transports_dict[var].Vin .-Vout    ), label =var , c = colors[i]); 
axes[1].plot(tecco[1:nt], transports_dict[var].Vinout, c = colors[i], alpha = 0.5); 

axes[2].plot(tecco[1:nt], transports_dict[var].Wtop, label = var, c = colors[i]); 
axes[3].plot(tecco[1:nt], transports_dict[var].Wbot, label = var, c = colors[i]); 
axes[4].plot(tecco[1:nt], conv, label = var, c = colors[i]); 
axes[4].plot(tecco[1:nt], conv_true, c = colors[i], alpha = 0.5, linewidth = 2); 
axes[1].set_title( L"Vθ_{in} - Vθ_{out}"); axes[2].set_title( L"Wθ_{out}"); axes[3].set_title( L"Wθ_{in}"); axes[4].set_title( "Flux Convergence")
fig.suptitle(region * " Advective Fluxes [z = 2 - 3 km]")
# axes[1].set_ylim(-10, 10)
fig

fig, axes = plt.subplots(1, 1, figsize = (15, 4), sharey = true)
axes.plot(tecco[1:nt], cumsum(conv .* 2.628e+6), label = var, c = colors[i]); 
axes.plot(tecco[1:nt], cumsum(conv_true.* 2.628e+6), c = colors[i], alpha = 0.5, linewidth = 2); 
fig

fig, axes = plt.subplots(2, 1, figsize = (10, 15), sharey = true)
axes[1].plot(tecco[1:nt], transports_dict[var].Vin, label =L"Vθ_{in}", c = "red"); 
axes[1].plot(tecco[1:nt], -Vout, label =L"Vθ_{out}", c = "red", alpha = 0.5); 
axes[1].plot(tecco[1:nt], -transports_dict[var].Wtop, label = L"-Wθ_{out}", c = "blue"); 
axes[1].plot(tecco[1:nt], transports_dict[var].Wbot, label = L"Wθ_{in}", c = "blue", alpha = 0.5); 
axes[2].plot(tecco[1:nt], conv, label = "Convergence", c = colors[i]); 
fig
# axes[1].set_title( L"V_{in} - V_{out}"); axes[2].set_title( L"W_{out}"); axes[3].set_title( L"W_{in}"); axes[4].set_title( "Flux Convergence")
fig.suptitle(region * " Advective Fluxes [z = 2 - 3 km]")
# axes[1].set_ylim(-20, 20)
axes[1].legend()
fig
fig.savefig(plotsdir("native/" * region * "_HeatBudgetFaces_" * shortnames[var]* "AdvectionFluxes_withoutBering" * "_" * suffix * ".png"), 
                     dpi = 400, bbox_inches = "tight")

# fig.savefig(plotsdir("native/" * region * var ), dpi = 400, bbox_inches = "tight")