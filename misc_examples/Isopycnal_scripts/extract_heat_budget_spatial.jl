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
suffix = "2to3"
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)

#create volume mask
area = readarea(γ)
cell_depths = get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));
cell_volumes_flat = sum_vertical(cell_volumes[:, lvls], γ)
cell_volumes_nan = deepcopy(cell_volumes_flat)
for ijk in eachindex(cell_volumes_nan)
    cell_volumes_nan[ijk][cell_volumes_nan[ijk] .== 0] .= NaN 
end
H = sum_vertical(cell_depths[:, lvls], γ); 
[H[ijk][iszero.(H[ijk])] .= Inf for ijk in eachindex(H)]
distance(λ1, λ2, ϕ1, ϕ2) = (λ1 - λ2)^2 + (ϕ1 - ϕ2)^2

tecco = 1992+1/24:1/12:2018
nz = 50; tstart = 12*3; tstop = 312; nt = tstop - tstart + 1;
E,F = trend_matrices(tecco[tstart:tstop])

# #pre-allocate
# GTF = OHC_helper.get_geothermalheating(γ, Γ)
zero_ma() = (β = MeshArray(γ,Float32); fill!(β, 0.0))

function filter_heat_budget(savename, 
                            diagpath::Dict{String, String}, 
                            expname::String, 
                            γ::gcmgrid)
    # if isfile(savename)
    #     print("File already saved")
    #     dθ = load(datadir(expname * region * "_THETA_budget_z.jld2"), "dθ")
    #     return dθ
    # end
    filelist = searchdir(diagpath[expname],"trsp_3d_set2")
    datafilelist_H  = filter(x -> occursin("data",x),filelist) 
    filelist = searchdir(diagpath[expname],"trsp_3d_set3")
    datafilelist_R  = filter(x -> occursin("data",x),filelist)
    filelist = searchdir(diagpath[expname],"state_3d_set1") 
    datafilelist_θ  = filter(x -> occursin("data",x),filelist)

    var_names = ["dθ", "κxyθ", "uvθ", "wθ", "κzθ"]
    vars = Dict(varname => zero_ma() for varname in var_names)
    κθ_conv3D = MeshArray(γ,Float32,50); uθ_conv3D = MeshArray(γ,Float32,50);
    i = 0 
    @time for tt = tstart:tstop
        i+=1 
        println(tt)
        fnameθ = datafilelist_θ[tt]
        fnameH = datafilelist_H[tt]
        fnameR = datafilelist_R[tt]
        #horizontal convergences
        @time θ_gcm = γ.read(diagpath[expname]*fnameθ,MeshArray(γ,Float32,50))

        dθλ = γ.read(diagpath[expname]*fnameH,MeshArray(γ,Float32,200))
        κθx = dθλ[ :, 1:50]; κθy = dθλ[:, 51:100]
        uθx = dθλ[:, 101:150]; uθy = dθλ[:, 151:200]; dθλ = nothing 
        OHC_helper.calc_UV_conv3D!(κθx, κθy, κθ_conv3D); 
        OHC_helper.calc_UV_conv3D!(uθx, uθy, uθ_conv3D);

        #vertical convergences
        dθr = γ.read(diagpath[expname]*fnameR,MeshArray(γ,Float32,150))
        wθ_conv = dθr[:, 1:50];
        κzθ_conv = dθr[ :, 51:100] .+ dθr[ :, 101:150]; dθr = nothing 

        for ff=1:5, k=1:nz-1
            wθ_conv.f[ff, k] .= wθ_conv.f[ff, k+1] .- wθ_conv.f[ff, k] #in - out 
            κzθ_conv.f[ff, k] .= κzθ_conv.f[ff, k+1] .- κzθ_conv.f[ff, k] #does not do correct heatflux for z = 50 
        end

        κzθ_conv_flat =  sum_vertical(κzθ_conv[:, lvls], γ)
        uθ_conv3D_flat =  sum_vertical(uθ_conv3D[:, lvls], γ)
        wθ_conv_flat =  sum_vertical(wθ_conv[:, lvls], γ)
        κθ_conv3D_flat =  sum_vertical(κθ_conv3D[:, lvls], γ)

        θ_gcm_flat = depth_average(θ_gcm[:, lvls], cell_depths[:, lvls], H, γ)

        #normalize fluxes
        for ff in 1:5
            vars["dθ"][ff]   .+= F[2,i] .* θ_gcm_flat[ff] ./ 3.156e+7 #C / year -> C / s
            vars["κxyθ"][ff] .+= (κθ_conv3D_flat[ff] ./ cell_volumes_nan[ff]) ./ nt 
            vars["uvθ"][ff]  .+= (uθ_conv3D_flat[ff] ./ cell_volumes_nan[ff]) ./ nt 
            vars["wθ"][ff]   .+= (wθ_conv_flat[ff] ./ cell_volumes_nan[ff]) ./ nt 
            vars["κzθ"][ff]  .+= (κzθ_conv_flat[ff] ./ cell_volumes_nan[ff]) ./ nt 
        end

    end

    return vars
end
sum_fluxes(d) = sum(d[varn] for varn in ["κxyθ", "uvθ", "wθ", "κzθ"])
 
expname = "iter0_bulkformula"
svname = datadir("THETA_budget_spatial" * "_" * expname * "_" * suffix *".jld2")
dθ = filter_heat_budget(svname, diagpath, expname, γ)
jldsave(svname, dθ = dθ)

