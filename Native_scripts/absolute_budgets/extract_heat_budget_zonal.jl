#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall, .OHC_helper
import NaNMath as nm
import PyPlot as plt
include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");

#using . is faster 
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());

# abbreviations for each experiment for labels, etc.
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
nexp = length(shortnames) # number of experiments
 
ocean_mask = OHC_helper.wet_pts(Γ)
region = "PAC56"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));

tecco = 1992+1/24:1/12:2018
nz = 50; tstart = 12*3; tstop = 312; nt = tstop - tstart + 1;
E,F = trend_matrices(tecco[tstart:tstop])
#pre-allocate
# GTF = OHC_helper.get_geothermalheating(γ, Γ)
# GTF_zonal = ma_zonal_avg(GTF, cell_volumes)
ΔV = OHC_helper.ma_zonal_sum(cell_volumes)
ΔV[ΔV .== 0] .= NaN

function filter_heat_budget(diagpath::Dict{String, String}, 
                            expname::String, 
                            γ::gcmgrid)
    filelist = searchdir(diagpath[expname],"trsp_3d_set2")
    datafilelist_H  = filter(x -> occursin("data",x),filelist) 
    filelist = searchdir(diagpath[expname],"trsp_3d_set3")
    datafilelist_R  = filter(x -> occursin("data",x),filelist)
    filelist = searchdir(diagpath[expname],"state_3d_set1") 
    datafilelist_θ  = filter(x -> occursin("data",x),filelist)

    var_names = ["dθ", "κxyθ", "uvθ", "wθ", "κzθ"]
    vars = Dict(varname => zeros(size(ΔV)) for varname in var_names)
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

        θ_gcm_zonal = OHC_helper.ma_zonal_avg(θ_gcm, cell_volumes)
        wθ_conv = wθ_conv .* PAC_msk
        κzθ_conv = κzθ_conv .* PAC_msk
        κθ_conv3D = κθ_conv3D .* PAC_msk; 
        uθ_conv3D = uθ_conv3D .* PAC_msk
        
        #
        vars["dθ"]   .+= F[2,i] .* θ_gcm_zonal ./ 3.156e+7 #C / year -> C / s
        vars["κxyθ"] .+= (OHC_helper.ma_zonal_sum(κθ_conv3D) ./ ΔV) ./nt
        vars["uvθ"]  .+= (OHC_helper.ma_zonal_sum(uθ_conv3D) ./ ΔV) ./nt
        vars["wθ"]   .+= (OHC_helper.ma_zonal_sum(wθ_conv) ./ ΔV) ./nt
        vars["κzθ"]  .+= (OHC_helper.ma_zonal_sum(κzθ_conv) ./ ΔV) ./nt

    end

    return vars
end

sum_fluxes(d) = sum(d[varn] for varn in ["κxyθ", "uvθ", "wθ", "κzθ"])
dθs = Dict()
for expname in ["iter129_bulkformula", "iter0_bulkformula"]
    dθ = filter_heat_budget(diagpath, expname, γ)
    dθs[expname] = deepcopy(dθ)
end

svename = datadir("native/THETA_budget_zonal_" * region *"_1995_2017.jld2")
jldsave(svename, dθs = dθs)
