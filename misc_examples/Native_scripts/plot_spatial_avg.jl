#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
using .OHC_helper
import NaNMath as nm
import PyPlot as plt
include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");

(ϕ,λ) = latlonC(γ)
tecco = 1992+1/24:1/12:2018
runpath,diagpath = listexperiments(exprootdir());

ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

#load in latitude mask 
ocean_mask = wet_pts(Γ)
region = "PAC56"
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
suffix = "2pt5km"
uplvl = -2.4e3; botlvl = -2.6e3
lvls = findall( botlvl .<= z[:].<= uplvl)
#create volume mask
area = readarea(γ)
cell_depths = OHC_helper.get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));

#need these to do a depth average
# H = OHC_helper.sum_vertical(cell_depths[:, lvls], γ); 
H = deepcopy(cell_depths[:, lvls])
[H[ijk][iszero.(H[ijk])] .= Inf for ijk in eachindex(H)]

runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)
nz = length(z)
nt = length(tecco)
P = MeshArray(γ,Float32, nz)
for ijk in eachindex(P)
    P[ijk] .= -pstdz[ijk[2]]
end

colors = [cmo.delta, cmo.dense, cmo.balance];
titles = ["Salinity", "Density", "Temperature"]
p₀ = 2000
ρθSavg = Dict()
for expname in ["iter129_bulkformula", "iter0_bulkformula"]
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist_UVW  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    spatial_avg = Dict()

    start = 36
    ntt = length(tecco) - start + 1

    θz_mean = MeshArray(γ,Float32, nz); fill!(θz_mean, 0.0)
    σz_mean = MeshArray(γ,Float32, nz); fill!(σz_mean, 0.0)
    Sz_mean = MeshArray(γ,Float32, nz); fill!(Sz_mean, 0.0)
    U_mean = MeshArray(γ,Float32, nz); fill!(U_mean, 0.0)
    V_mean = MeshArray(γ,Float32, nz); fill!(V_mean, 0.0)

    for tt in 36:nt
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        θname = datafilelist_θ[tt]
        UVWname = datafilelist_UVW[tt]

        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        @time θSz = γ.read(diagpath[expname]*θname,MeshArray(γ,Float32,2*nz))
        θz = θSz[:, 1:nz]; Sz = θSz[:, nz+1:end]

        @time UVW = γ.read(diagpath[expname]*UVWname,MeshArray(γ,Float32,2*nz))
        U = UVW[:, 1:nz]; V = UVW[:, nz+1:end]
        U, V = rotate_uv(U, V, Γ)
        for ijk in eachindex(P)
            σ = OHC_helper.densityJMD95.(θz.f[ijk],Sz.f[ijk], P[ijk], p₀) #EOS from MITGCM 
            σ = σ .- 1000
            θz_mean[ijk] .+= θz.f[ijk] ./ ntt
            Sz_mean[ijk] .+= Sz.f[ijk] ./ ntt
            σz_mean[ijk] .+= σ ./ ntt
            U_mean[ijk] .+= U.f[ijk] ./ntt
            V_mean[ijk] .+= V.f[ijk] ./ntt
        end 

    end

    # θ_zonal = OHC_helper.depth_average(θz_mean[:, lvls], cell_depths[:, lvls], H, γ)
    # sigma_zonal = OHC_helper.depth_average(σz_mean[:, lvls], cell_depths[:, lvls], H, γ)
    # S_zonal = OHC_helper.depth_average(Sz_mean[:, lvls], cell_depths[:, lvls], H, γ)
    # U_zonal = OHC_helper.depth_average(U_mean[:, lvls], cell_depths[:, lvls], H, γ)
    # V_zonal = OHC_helper.depth_average(V_mean[:, lvls], cell_depths[:, lvls], H, γ)
    
    θ_zonal = deepcopy(θz_mean[:, lvls])
    sigma_zonal = deepcopy(σz_mean[:, lvls])
    S_zonal = deepcopy(Sz_mean[:, lvls])
    U_zonal = deepcopy(U_mean[:, lvls])
    V_zonal = deepcopy(V_mean[:, lvls])

    spatial_avg["θ"]  =  θ_zonal .* ocean_mask
    spatial_avg["σ2"] =  sigma_zonal .* ocean_mask
    spatial_avg["S"]  =  S_zonal .* ocean_mask
    spatial_avg["U"]  =  U_zonal .* ocean_mask
    spatial_avg["V"]  =  V_zonal .* ocean_mask

    ρθSavg[expname] = spatial_avg
end

svename = datadir("native/native_sigma2_spatialavg_" * suffix *"_1995_2017.jld2")
jldsave(svename, ρθSavg = ρθSavg)

iter129 = load(svename)["ρθSavg"]["iter129_bulkformula"]
iter0 = load(svename)["ρθSavg"]["iter0_bulkformula"]

sigma_plot = Dict()
sigma_plot["Iteration 129"] = iter129["σ2"]
sigma_plot["Iteration 0"] = iter0["σ2"]
sigma_plot["Iteration 129 minus Iteration 0"] = sigma_plot["Iteration 129"] .- sigma_plot["Iteration 0"]

velocity_plot = Dict()
velocity_plot["Iteration 129"] = [deepcopy(iter129["U"]), deepcopy(iter129["V"])]
velocity_plot["Iteration 0"] = (iter0["U"], iter0["V"])
velocity_plot["Iteration 129 minus Iteration 0"] = velocity_plot["Iteration 129"] .- velocity_plot["Iteration 0"]
ocn_reg = LLCcropC(ocean_mask,γ); PAC_reg =  LLCcropC(PAC_msk,γ)
reg_λ = LLCcropC(λ,γ); reg_ϕ = LLCcropC(ϕ,γ)
levels = 36.8:0.02:37.; bounds = extrema(levels) 

proj0 = ECCOonPoseidon.cartopy.crs.Robinson(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

fig, axs = plt.subplots(1, 3, figsize=(30,10), subplot_kw=Dict("projection"=> proj0))
CF = Any[]

function plot_regular_field(data, udata, vdata, ax, cmap, bounds, i)
    # data = data .* PAC_reg
    udata = udata .* PAC_reg
    vdata = vdata .* PAC_reg
    data[data .== 0] .= NaN
    udata[udata .== 0] .= NaN
    vdata[vdata .== 0] .= NaN

    # #normalize the speeds for easy plotting
    speed = maximum(sqrt.((udata.^2) .+ (vdata.^2)))
    lw = 5*speed / nm.maximum(speed)
    # udata = udata ./ speed
    # vdata = vdata ./ speed
    println(nm.extrema(data))
    cf = ax.pcolormesh(reg_λ, reg_ϕ,  data, transform=projPC, 
    cmap = cmap, vmin = bounds[1], vmax = bounds[2]) 
    spacing = 2
    ax.streamplot(x = reg_λ[1:spacing:end, 1:spacing:end], y = reg_ϕ[1:spacing:end, 1:spacing:end],  
    u = udata[1:spacing:end, 1:spacing:end], v = vdata[1:spacing:end, 1:spacing:end], 
    transform=projPC, linewidth=lw)
    ax.set_aspect("equal")
    ax.coastlines(resolution="110m")
    ax.set_extent((-160, 120, -70, 56),crs=projPC)
    gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.top_labels = false
    gl.bottom_labels = true

    gl.right_labels = false
    if i != 1
        gl.left_labels = false
    end
    return cf
end
bounds_list = [bounds, bounds, (-0.03, 0.03)]
cmaps = [cmo.deep, cmo.deep, cmo.balance];
levels_list = [levels, levels, -0.3:0.05:0.3]
for (i, key) in enumerate(["Iteration 129", "Iteration 0", "Iteration 129 minus Iteration 0"])
    ax = axs[i]
    data = LLCcropC(sigma_plot[key], γ)
    u_data = LLCcropC(velocity_plot[key][1], γ)
    v_data = LLCcropC(velocity_plot[key][2], γ)

    cf = plot_regular_field(data, u_data, v_data, ax, cmaps[i], bounds_list[i], i)
    push!(CF, cf)
    ax.set_title(key)
end
fig
fig.subplots_adjust(wspace=0.05)
fig.colorbar(CF[1], ax = axs[1:2], orientation = "horizontal", fraction = 0.04, 
label = L"[kg/m³]", pad=0.1, extend = "both")
fig.colorbar(CF[3], ax = axs[3], orientation = "horizontal", fraction = 0.04, 
label = L"[kg/m³]", pad=0.1, extend = "both")
fig.suptitle("Time-Mean Density in ECCO (colors) \n Normalized Velocity Field (arrows) \n [z=2-3km, t = 1995-2017]", y=0.7)
fig
fig.savefig(plotsdir("DensityandVelocity" * suffix * ".png"), bbox_inches = "tight")