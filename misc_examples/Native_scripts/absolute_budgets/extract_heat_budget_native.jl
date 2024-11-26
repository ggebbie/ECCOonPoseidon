#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, Plots, LaTeXStrings, FLoops
using .OHC_helper
using Plots
using ColorSchemes
using Plots.PlotMeasures

include(srcdir("config_exp.jl"))

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
region = "NPAC"; 
NPAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
region = "PAC56"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")

region = "PAC56noNPAC"
NoNPAC_msk = PAC_msk .- NPAC_msk
msk = NoNPAC_msk

cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));
sum(cell_volumes)
H = smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H
tecco = 1992+1/24:1/12:2018
#pre-allocate
GTF = OHC_helper.get_geothermalheating(cell_depths, Γ, γ)
GTF = GTF .* msk
nz = 50
GTF_z=  zeros(Float32, nz)

ΔV = zeros(Float32, nz)
for k=1:nz
    ΔV[k] = Float32(sum(cell_volumes[:, k]))
    # GTF_z[k] = sum(GTF[:, k])
end
GTF_z = repeat(GTF_z, 1, 313)
lvls = findall( -3000 .<= z[:].<= -2000)
masked_volume = cell_volumes[:, lvls]
sum_masked_volume = Float32(sum(masked_volume))

function heat_flux_profile(ds::MeshArray, ΔV)
    nz = size(ds, 2)
    vol_avg = zeros(Float32, nz)

    for ff=1:5, k=1:nz
        vol_avg[k] += Float32(sum(ds[ff, k])) ./ ΔV[k]
    end
    return vol_avg
end


function filter_heat_budget(savename, diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)
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
    filelist = searchdir(diagpath[expname], "state_2d_set1") # first filter for state_3d_set1
    datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
 
    var_names = ["θ", "κxyθ", "uvθ", "wθ", "κzθ"]
    nt = length(datafilelist_R); nz = 50
    vars = Dict(varname => zeros(Float32, nz, nt) for varname in var_names)
    κθ_conv3D = MeshArray(γ,Float32,50); uθ_conv3D = MeshArray(γ,Float32,50);

    @time for tt = 1:nt
        println(tt)
        fnameθ = datafilelist_θ[tt]
        fnameH = datafilelist_H[tt]
        fnameR = datafilelist_R[tt]
        fnameS = datafilelist_S[tt]

        #horizontal convergences
        sθ = OHC_helper.extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_H)
        #  sθ = OHC_helper.extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_depths)
        vars["θ"][:, tt]  = heat_flux_profile(sθ .* cell_volumes, ΔV); 
        dθλ = γ.read(diagpath[expname]*fnameH,MeshArray(γ,Float32,200))
        κθx = dθλ[ :, 1:50]; κθy = dθλ[:, 51:100]
        uθx = dθλ[:, 101:150]; uθy = dθλ[:, 151:200]; dθλ = nothing 
        OHC_helper.calc_UV_conv3D!(κθx, κθy, κθ_conv3D); OHC_helper.calc_UV_conv3D!(uθx, uθy, uθ_conv3D);
        κθ_conv3D = κθ_conv3D .* msk; 
        uθ_conv3D = uθ_conv3D .* msk

        vars["κxyθ"][:, tt] .= heat_flux_profile(κθ_conv3D, ΔV); 
        vars["uvθ"][:, tt] .= heat_flux_profile(uθ_conv3D, ΔV); 

        #vertical convergences
        dθr = γ.read(diagpath[expname]*fnameR,MeshArray(γ,Float32,150))
        wθ_conv = dθr[:, 1:50];
        κzθ_conv = dθr[ :, 51:100] .+ dθr[ :, 101:150]; dθr = nothing 

        for ff=1:5, k=1:nz-1
            wθ_conv.f[ff, k] .= wθ_conv.f[ff, k+1] .- wθ_conv.f[ff, k] #in - out 
            κzθ_conv.f[ff, k] .= κzθ_conv.f[ff, k+1] .- κzθ_conv.f[ff, k] #does not do correct heatflux for z = 50 
        end

        wθ_conv = wθ_conv .* msk
        κzθ_conv = κzθ_conv .* msk
        vars["wθ"][:, tt] .= heat_flux_profile(wθ_conv, ΔV); 
        vars["κzθ"][:, tt] .= heat_flux_profile(κzθ_conv, ΔV)
    end

    return vars
end

for expname in keys(shortnames)
    savename = datadir(expname * region * "_THETA_budget_z.jld2")
    dθ = filter_heat_budget(savename, diagpath, expname, γ)
    jldsave(savename, dθ = dθ)
end

flux_names = ["κxyθ", "uvθ", "wθ", "κzθ"]
d_integrate = Dict()
z_col = zeros(Float32, 50, 1)
[d_integrate[var] = cumsum(hcat(z_col, dθ[var] .* 2.628e+6), dims=2) for var in flux_names]
d_integrate["θ"]= dθ["θ"]
d_integrate["GTF"] = GTF_z
flux_names = ["κxyθ", "uvθ", "wθ", "κzθ",  "GTF"]

d_integrate["sum"] = sum(d_integrate[varn] for varn in flux_names)

flux_names = ["κxyθ", "uvθ", "wθ", "κzθ", "sum"]

anomaly(x, dim) = x .- mean(x, dims = dim)

cmax = maximum([maximum(100 *abs.(anomaly(v[10:end, :], 2))) for v in values(d_integrate)] )
cmax = 1701 .+ 0.1
levels = -4:0.5:4
levels = vcat([-cmax], levels, [cmax]) 

plts = []
labels = [L"\theta_{\kappa XY}", L"\theta_{UV}", L"\theta_{W}", L"\theta_{\kappa Z}", 
"sum"]
labels[end] = L"\theta_{\kappa XY} + \theta_{UV} + \theta_{W} + \theta_{\kappa Z}"
for (i, var) in enumerate(flux_names)    
    var_anom = anomaly(d_integrate[var][:, 1:end-1], 2)
    p = plot([NaN], [NaN], label = nothing, xlims = extrema(tecco), 
    title = labels[i])
    contourf!(p, tecco, abs.(z[20:end]), 100 .* var_anom[20:end, :], 
    color = :balance, levels = levels,   
    yflip = true, ylims = (500, 3500), 
    lw=0.5, clabels=true, cbar=false)
    push!(plts, p)
end

var = "θ"; i = 6
var_anom = var_anom = anomaly(d_integrate[var], 2)
p = plot([NaN], [NaN], label = nothing, xlims = extrema(tecco), 
title = "Actual θ ")
contourf!(p, tecco, abs.(z[20:end]), 100 .* var_anom[20:end, :], 
color = :balance, levels = levels,   
yflip = true, ylims = (500, 3500), 
lw=0.5, clabels=true, cbar=false)
push!(plts, p)

label_ys!(p) = ylabel!(p, "Depth [m]")
label_xs!(p) = xlabel!(p, "Year")
label_ys!.([plts[1], plts[4]])
label_xs!.(plts[4:6])

p = plot(title = region * " Temperature Anomaly [cK] in ECCO",
framestyle=nothing,showaxis=false,xticks=false,yticks=false,margin=0mm)
plts2 = cat(p, plts, dims = 1)

l = @layout [a{0.01h}; grid(2,3);]
p2 = plot(plts2..., layout=l, size = (1000, 650), link = :all)
p2
savefig(p2, plotsdir(region * "_HeatBudgetAnoms.png"))
