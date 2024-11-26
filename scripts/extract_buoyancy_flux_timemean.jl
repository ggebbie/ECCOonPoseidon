#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl
using Pkg
Pkg.activate(".")
include("../src/intro.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, LaTeXStrings, 
PyCall
import PyPlot as plt
import NaNMath as nm
              
include(srcdir("config_exp.jl"))

tecco= 1992+1/24:1/12:2018 # ecco years
nt = length(tecco)

runpath,diagpath = listexperiments(exprootdir());
bflux_dict = Dict();

for expname in ["iter0_bulkformula"]

    bflux_dict[expname] = Dict()
    bflux_dict[expname]["TFLUX"] = MeshArray(γ,Float64); fill!(bflux_dict[expname]["TFLUX"], 0.0)
    bflux_dict[expname]["oceFWflx"] = MeshArray(γ,Float64); fill!(bflux_dict[expname]["oceFWflx"], 0.0)
    bflux_dict[expname]["oceQsw"] = MeshArray(γ,Float64); fill!(bflux_dict[expname]["oceQsw"], 0.0)
    bflux_dict[expname]["oceSflux"] = MeshArray(γ,Float64); fill!(bflux_dict[expname]["oceSflux"], 0.0)

    filelist = searchdir(diagpath[expname],"budg2d_zflux_set1") # first filter for state_3d_set1
    τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    filelist2 = searchdir(diagpath[expname],"budg2d_zflux_set2") # first filter for state_3d_set1
    τdatafilelist2 = filter(x -> occursin("data",x),filelist2) # second filter for "data"

    for tt = 1:nt
        println(tt)
        x = γ.read(diagpath[expname] * τdatafilelist[tt],MeshArray(γ,Float64, 7));
        y = γ.read(diagpath[expname] * τdatafilelist[tt],MeshArray(γ,Float64, 5));

        oceFWflx = x[:, 1]
        TFLUX = x[:, 3]
        oceQsw = x[:, 6]        
        oceSflux = y[:, 5]
        for ff = 1:5
            bflux_dict[expname]["TFLUX"].f[ff] .+= TFLUX.f[ff] ./ nt
            bflux_dict[expname]["oceFWflx"].f[ff] .+= oceFWflx.f[ff] ./ nt
            bflux_dict[expname]["oceQsw"].f[ff] .+= oceQsw.f[ff] ./ nt
            bflux_dict[expname]["oceSflux"].f[ff] .+= oceSflux.f[ff] ./ nt

        end
    end
end

i0bflux_dict = bflux_dict["iter0_bulkformula"]

for key in keys(i0bflux_dict)
    write(datadir(key * "_i0_mean.data"), Float32.(i0bflux_dict[key]))
end