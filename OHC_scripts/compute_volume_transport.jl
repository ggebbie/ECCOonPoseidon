include("../src/intro.jl")
include("../src/OHC_helper.jl")
using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings
using .OHC_helper

suffix = "_full"
fname = "UVELMASS"
savename = "UTRSP"
filedir = "ECCO_vars/"
include(srcdir("config_exp.jl"))

const tecco = collect(Float64, 1992+1/24:1/12:2018)
runpath,diagpath = listexperiments(exprootdir())
shortnames = expnames()
ignore_list = ["noIA", "129ff"]
shortnames = reduce_dict(shortnames, ignore_list)

#Mass transport in dir = (dir)velmass *  (dir)Area
fname1 = "VTRSP"
fname2 = "UTRSP"
savename = "NTRSP"
ma_zeros = MeshArray(γ, Float64, 50)
ma_zeros.=0

for ijk in eachindex(xArea)
    xArea.f[ijk] .= Γ.DXG.f[ijk[1]] .* ΔzF[ijk[2]]
    yArea.f[ijk] .= Γ.DYG.f[ijk[1]] .* ΔzF[ijk[2]]
end
#Mass transport in dir = (dir)velmass *  (opp dir)Area

for (key,values) in shortnames
    expname = key
    outdir = filedir * expname * "/"
    mkpath(datadir(outdir))
    # @time UTRSP = load_object(datadir(filedir*fname1*"_"*expname *".jld2"))
    @time UTRSP = load_object(datadir(filedir*fname*"_"*expname *".jld2"))
    for tt in 1:length(tecco)
        jldsave(datadir(outdir * savename *"_full_"*string(tt)*".jld2"),
        UTRSP = @views UTRSP[tt] .* yArea)
    end
end


for (key,values) in shortnames
    expname = key
    println(expname)
    # @time UTRSP = load_object(datadir(filedir*fname1*"_"*expname *".jld2"))
    @time for tt in 1:length(tecco)
        vfname = datadir(filedir * expname * "/" * fname1 *"_full_"*string(tt)*".jld2")
        ufname = datadir(filedir * expname * "/" * fname2 *"_full_"*string(tt)*".jld2")
        UTRSP = load_object(ufname) #weighted UVELMASS
        VTRSP = load_object(vfname) #weighted VVELMASS

        conv_UTRSP = calc_UV_conv(UTRSP, ma_zeros)
        conv_VTRSP = calc_UV_conv(ma_zeros, VTRSP)
        ETRSP, NTRSP = ECCOtour.rotate_uv(conv_UTRSP, conv_VTRSP, Γ)
        jldsave(datadir(filedir * expname * "/" *
        savename *"_full_"*string(tt)*".jld2"),
        NTRSP = NTRSP)
    end
end



# for (key,values) in shortnames
#     expname = key
#     for tt in 1:length(tecco)
#         fdir = filedir * expname * "/" * fname2 * "_full_"*string(tt) *".jld2"
#         jldopen(datadir(fdir), "r") do HBUDG
#             savdir = filedir * expname * "/" * fname2 * suffix *string(tt) * ".jld2"
#             jldsave(datadir(savdir),
#             AdvR = HBUDG["AdvR"][:, lvls], 
#             AdvH = HBUDG["AdvH"][:, lvls], 
#             DiffR = HBUDG["DiffR"][:, lvls], 
#             DiffH = HBUDG["DiffH"][:, lvls])
#         end
#     end
# end



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