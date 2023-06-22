#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings,
Plots
gr()
import NaNMath as nm
using .OHC_helper
import PyPlot as plt
using PyCall


pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;
include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
area = readarea(γ)
runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
labels_L = [L"\theta^{129}", L"\theta^{\Delta F}",L"\theta^{\Delta T}", L"\theta^{0}"]
tecco = 1992+1/24:1/12:2018

ocean_mask = wet_pts(Γ)
region = "NPAC30"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "full")
suffix = "2to3"
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)

ϕ_mask_Inf = ϕ .* PAC_msk; ϕ_mask_Inf[findall(ϕ_mask_Inf .== 0 )] .= Inf
ϕ_mask_min = minimum(ϕ_mask_Inf)

area = readarea(γ)
area_masked = area .* PAC_msk
area_sum = sum(area_masked)

ϕ_S = (ϕ .≈ ϕ_mask_min)
ϕ_S = ϕ_S .* PAC_msk

cell_depths = get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
smush_depths = smush(cell_depths); 
cell_volumes = get_cell_volumes(area, cell_depths);

smush_depths[findall(smush_depths .== 0)] = Inf; 
inv_depths = 1 ./ smush_depths

Advective_Heat_Budgets = Dict()
masked_volume = cell_volumes[:, lvls]
sum_masked_volume = Float32(sum(masked_volume))
@time for expname in keys(shortnames)
    Advective_Heat_Budgets[expname] = Dict()
    filelist = searchdir(diagpath[expname],"trsp_3d_set2") # first filter for state_3d_set1
    datafilelist_H  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"trsp_3d_set3") # first filter for state_3d_set1
    datafilelist_R  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(datafilelist_R)
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    bot_lev = lvls[end]; top_lev = lvls[1]    
    sθ̄ = []; AdvH = Float32[]; AdvR = Float32[];DiffH = Float32[]; DiffR = Float32[]
    @time for tt=1:nt
        fnameH = datafilelist_H[tt]
        fnameR = datafilelist_R[tt]
        fnameθ = datafilelist_θ[tt]
        fnameS = datafilelist_S[tt]

        sθ = extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_depths)
        push!(sθ̄, Float32(sum(sθ[:, lvls].*masked_volume) / sum_masked_volume ))
        #northward diffustion and advection

        dθλ = γ.read(diagpath[expname]*fnameH,MeshArray(γ,Float32,200))
        κuθ = dθλ[:, lvls]; κvθ = dθλ[:, (50 .+lvls)]
        uθ = dθλ[:, (100 .+lvls)]; vθ = dθλ[:, (150 .+lvls)]

        κEθ, κNθ = rotate_uv(κuθ, κvθ, Γ)
        Eθ, Nθ = rotate_uv(uθ, vθ, Γ)

        push!(DiffH, Float32(sum(κNθ .* ϕ_S) / sum_masked_volume))
        push!(AdvH, Float32(sum(Nθ .* ϕ_S) /sum_masked_volume))
        
        #vertical advection and diffustion 
        wθ = γ.read(diagpath[expname]*fnameR,MeshArray(γ,Float32,150))
        wθ_conv = wθ[:, bot_lev+1] .- wθ[:, top_lev] #should already be masked
        κIwθ_conv = wθ[:, 50+bot_lev+1] .- wθ[:, 50+top_lev] #should already be masked
        κEwθ_conv = wθ[:, 100+bot_lev+1] .- wθ[:, 100+top_lev] #should already be masked
        κwθ_conv = κIwθ_conv + κEwθ_conv

        push!(DiffR, Float32(sum(κwθ_conv .* PAC_msk)  / sum_masked_volume))
        push!(AdvR, Float32(sum(wθ_conv .* PAC_msk)  / sum_masked_volume))

    end
    DiffH_sum = cumsum(cat(0.0, DiffH .* 2.628e+6, dims = 1))
    DiffR_sum = cumsum(cat(0.0, DiffR .* 2.628e+6, dims = 1))

    AdvH_sum = cumsum(cat(0.0, AdvH .* 2.628e+6, dims = 1))
    AdvR_sum = cumsum(cat(0.0, AdvR .* 2.628e+6, dims = 1))
    Advective_Heat_Budgets[expname]["DiffH"] = DiffH_sum
    Advective_Heat_Budgets[expname]["DiffR"] = DiffR_sum
    Advective_Heat_Budgets[expname]["AdvH"] = AdvH_sum
    Advective_Heat_Budgets[expname]["AdvR"] = AdvR_sum
    Advective_Heat_Budgets[expname]["θ"] = sθ̄
end
Heat_Budgets = Advective_Heat_Budgets
jldsave(datadir("HeatBudgetTrue_"*region*"_full_2to3.jld2"); Heat_Budgets)

# expname = "iter129_bulkformula"
# DiffH = Advective_Heat_Budgets["iter129_bulkformula"]["DiffH"]
# DiffR = Advective_Heat_Budgets["iter129_bulkformula"]["DiffR"]
# AdvH = Advective_Heat_Budgets["iter129_bulkformula"]["AdvH"]
# AdvR = Advective_Heat_Budgets["iter129_bulkformula"]["AdvR"]
# θ = Advective_Heat_Budgets[expname]["θ"]
# sum_forcing = θ[1] .+ DiffH .+ DiffR .+ AdvH .+ AdvR
# diff = (θ .- sum_forcing[1:end-1])
# m = (diff[end] - diff[1]) / length(tecco)
# geothermal(x) = m*(x - 1)
# geo = geothermal.(1:312)

# jldopen(datadir("Geothermal_"*region*"_full_2to3.jld2"); geo)

# plot(sum_forcing[1:end-1] .+ geo)
# plot!(θ)

# # 
