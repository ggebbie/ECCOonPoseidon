include("../src/intro.jl")
include("../src/OHC_helper.jl")
using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings
using .OHC_helper



filedir = "ECCO_vars/"
include(srcdir("config_exp.jl"))

const tecco = collect(Float64, 1992+1/24:1/12:2018)
runpath,diagpath = listexperiments(exprootdir())
shortnames = expnames()
ignore_list = ["noIA", "129ff"]
shortnames = reduce_dict(shortnames, ignore_list)
fname2 = "THETA_BUDG"
suffix1 = "_1to4_"
suffix2 = "_2tobot_"
suffix3 = "_surfto2_"


for (key,values) in shortnames
    expname = key
    for tt in 1:length(tecco)
        fdir = filedir * expname * "/" * fname2 * "_full_"*string(tt) *".jld2"
        jldopen(datadir(fdir), "r") do HBUDG
            uplvl = -1e3; botlvl = -4e3
            lvls = findall( botlvl .<= z[:].<= uplvl)
            savdir1 = filedir * expname * "/" * fname2 * suffix1 * string(tt) * ".jld2"
            jldsave(datadir(savdir1),
            AdvR = HBUDG["AdvR"][:, lvls], 
            AdvH = HBUDG["AdvH"][:, lvls], 
            DiffR = HBUDG["DiffR"][:, lvls], 
            DiffH = HBUDG["DiffH"][:, lvls])

            # uplvl = -2000; botlvl = -Inf
            # lvls = findall( botlvl .<= z[:].<= uplvl)
            # savdir2 = filedir * expname * "/" * fname2 * suffix2 * string(tt) * ".jld2"
            # jldsave(datadir(savdir2),
            # AdvR = HBUDG["AdvR"][:, lvls], 
            # AdvH = HBUDG["AdvH"][:, lvls], 
            # DiffR = HBUDG["DiffR"][:, lvls], 
            # DiffH = HBUDG["DiffH"][:, lvls])

        end
    end
end



# depths_array = zeros(length(lvls), length(tecco))
# for tt in 1:length(tecco)
    #@time h5read("datadir(outdir * savename *"_full_"*string(tt)*".jld2"", 
    #"AdvH")
    #@time load("datadir(outdir * savename *"_full_"*string(tt)*".jld2"", 
    #"AdvH")
    #for lvl in lvls
    # depths_array[AdvH][lvls, tt] = @views volume_average(Adv[:, lvls], cell_volumes(:, lvls)
    #end
 #end
# reconstruct()
#