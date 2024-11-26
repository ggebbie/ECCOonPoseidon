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

tecco = 1992+1/24:1/12:2018
E,F = trend_matrices(tecco)
nt = 36; nz = 50
zidx = 34; ptidx = 100

#create volume mask
ocean_mask = wet_pts(Γ)
area = readarea(γ)
cell_depths = get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));
H = OHC_helper.smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H

cell_volume = cell_volumes[4, zidx][ptidx]

function filter_heat_budget(diagpath::Dict{String, String}, 
                            expname::String, 
                            γ::gcmgrid)

    filelist = searchdir(diagpath[expname],"trsp_3d_set2")
    datafilelist_H  = filter(x -> occursin("data",x),filelist) 
    filelist = searchdir(diagpath[expname],"trsp_3d_set3")
    datafilelist_R  = filter(x -> occursin("data",x),filelist)
    filelist = searchdir(diagpath[expname],"state_3d_set1") 
    datafilelist_θ  = filter(x -> occursin("data",x),filelist)
    filelist = searchdir(diagpath[expname], "state_2d_set1") # first filter for state_3d_set1
    datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    var_names = ["θ", "κxyθ", "uvθ", "wθ", "κzθ"]
    vars = Dict(varname => zeros(nt) for varname in var_names)
    κθ_conv3D = MeshArray(γ,Float32,50); uθ_conv3D = MeshArray(γ,Float32,50);


    @time for tt = 1:nt
        println(tt)
        fnameθ = datafilelist_θ[tt]
        fnameH = datafilelist_H[tt]
        fnameR = datafilelist_R[tt]
        fnameS = datafilelist_S[tt]

        #horizontal convergences
        sθ = OHC_helper.extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_H)
        dθλ = γ.read(diagpath[expname]*fnameH,MeshArray(γ,Float32,200))
        κθx = dθλ[ :, 1:50]; κθy = dθλ[:, 51:100]
        uθx = dθλ[:, 101:150]; uθy = dθλ[:, 151:200]; dθλ = nothing 
        OHC_helper.calc_UV_conv3D!(κθx, κθy, κθ_conv3D); OHC_helper.calc_UV_conv3D!(uθx, uθy, uθ_conv3D);
        # κθ_conv3D = κθ_conv3D .* PAC_msk; 
        # uθ_conv3D = uθ_conv3D .* PAC_msk

        #vertical convergences
        dθr = γ.read(diagpath[expname]*fnameR,MeshArray(γ,Float32,150))
        wθ_conv = dθr[:, 1:50];
        κzθ_conv = dθr[ :, 51:100] .+ dθr[ :, 101:150]; dθr = nothing 

        for ff=1:5, k=1:nz-1
            wθ_conv.f[ff, k] .= wθ_conv.f[ff, k+1] .- wθ_conv.f[ff, k] #in - out 
            κzθ_conv.f[ff, k] .= κzθ_conv.f[ff, k+1] .- κzθ_conv.f[ff, k] #does not do correct heatflux for z = 50 
        end

        vars["θ"][tt]   =  sθ[4, zidx][ptidx] #C / year
        vars["κxyθ"][tt] = κθ_conv3D[4, zidx][ptidx] ./ cell_volume #C/ s -> C/ year
        vars["uvθ"][tt]  = uθ_conv3D[4, zidx][ptidx] ./ cell_volume
        vars["wθ"][tt]   = wθ_conv[4, zidx][ptidx] ./ cell_volume
        vars["κzθ"][tt]  = κzθ_conv[4, zidx][ptidx] ./ cell_volume

    end

    return vars
end
sum_fluxes(d) = sum(d[varn][:] for varn in ["κxyθ", "uvθ", "wθ", "κzθ"])

expname = "iter129_bulkformula"
dθ = filter_heat_budget(diagpath, expname, γ)

dθ["dθ"] = (dθ["θ"][2:end] .- dθ["θ"][1:end-1]) ./ 2.628e+6
dθ["sum"] = Float32.(sum_fluxes(dθ))
tecco_int = (tecco[2:nt] .+ tecco[1:nt-1]) ./ 2

fig, axes = plt.subplots(figsize=(10,5))
axes.plot(tecco[1:nt], dθ["sum"], label = "reconstruction")
axes.plot(tecco_int, dθ["dθ"], label = "true")
axes.legend(); axes.set_title("dθ/dt")
fig
println(mean(dθ["dθ"])) #off by an order of magnitude again...? 
println(mean(dθ["sum"]))