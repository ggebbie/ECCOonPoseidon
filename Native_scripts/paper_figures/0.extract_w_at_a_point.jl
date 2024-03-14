include("../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
import DataFrames as DF
import PyPlot as plt

include(srcdir("plot_and_dir_config.jl"))
@pyimport matplotlib.patches as patches

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)

cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
ΔV = lateral_sum(cell_volumes)

tecco = 1992+1/24:1/12:2018; nz = 50
get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )

γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
Γ=GridLoad(γ; option="full")

# lon=[i for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
# lat=[j for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
lon=[-144.9]
lat=[50.1]
(f,i,j,w,j_f,j_x,j_y)=InterpolationFactors(Γ,vec(lon),vec(lat))
function get_temperature(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid, cell_volumes)


    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    datafilelist_τ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist_w  = filter(x -> occursin("data",x),filelist) # second filter for "data"


    nt = length(datafilelist_τ);
    nϕ = length(ϕ_avg)
    
    ma_template = MeshArray(γ,Float32,50)
    ts = Any[]
    @time for tt = 1:nt
        println(tt)
        fnameuvw = datafilelist_w[tt]
        Tname = datafilelist_τ[tt]
        τx, τy = extract_ocnTAU(diagpath, expname , Tname, γ)
        τcurl = MeshArrays.curl(τx, τy, Γ)
        
        u, v, w = extract_eulerian_velocities(diagpath, expname, fnameuvw, γ)

        DD=Interpolate(w[:, 40],f,i,j,w)
        push!(ts, DD * 1)
    end

    return ts
end

adjust_exps = Dict()
adjust_exps["iter0_bulkformula"] = get_temperature(diagpath, "iter0_bulkformula", γ, cell_volumes)
# adjust_exps["only_kappa"] = get_temperature(diagpath, "only_kappa", γ, cell_volumes)
# adjust_exps["only_init"] = get_temperature(diagpath, "only_init", γ, cell_volumes)
# adjust_exps["only_sfc"] = get_temperature(diagpath, "only_sfc", γ, cell_volumes)
# adjust_exps["only_wind"] = get_temperature(diagpath, "only_wind", γ, cell_volumes)
# adjust_exps["only_buoyancy"] = get_temperature(diagpath, "only_buoyancy", γ, cell_volumes)
adjust_exps["iter129_bulkformula"] = get_temperature(diagpath, "iter129_bulkformula", γ, cell_volumes)

fig, ax = plt.subplots(figsize = (10, 10))
iter129 = [a[1] for a in adjust_exps["iter129_bulkformula"]]
iter0 = [a[1] for a in adjust_exps["iter0_bulkformula"]]

using Pandas, CSV
df = DF.DataFrame(Dict("Iteration 129" => iter129, "Iteration 0" => iter0))
CSV.write(datadir("iter0iter129_papa.csv"), df) 
ax.plot(tecco, iter129, label = "iter129")
ax.plot(tecco, iter0, label = "iter0")
ax.legend()
ax.set_xlim(2013, 2018)
fig