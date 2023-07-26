#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, LaTeXStrings, 
PyCall
import PyPlot as plt
import NaNMath as nm
              
include(srcdir("config_exp.jl"))

tecco= 1992+1/24:1/12:2018 # ecco years
nt = length(tecco)

runpath,diagpath = listexperiments(exprootdir());
τx_dict = Dict(); τy_dict = Dict();
vars = ["iter129_bulkformula", "iter0_bulkformula"]

for expname in vars
    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    τx_dict[expname] = MeshArray(γ,Float32); fill!(τx_dict[expname], 0.0)
    τy_dict[expname] = MeshArray(γ,Float32); fill!(τy_dict[expname], 0.0)

    for tt = 1:nt
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        Tname = τdatafilelist[tt]
        
        τx, τy = OHC_helper.extract_ocnTAU(diagpath, expname , τdatafilelist, tt, γ)

        for ff = 1:5
            τx_dict[expname].f[ff] .+= τx.f[ff] ./ nt
            τy_dict[expname].f[ff] .+= τy.f[ff] ./ nt
        end
    end
end

Δdict = Dict()
Δdict["oceTAUX"] = τx_dict["iter129_bulkformula"] .- τx_dict["iter0_bulkformula"]
Δdict["oceTAUY"] = τy_dict["iter129_bulkformula"] .- τy_dict["iter0_bulkformula"]


write(datadir("oceTAUX_i129_i0_diff.data"),Δdict["oceTAUX"])
write(datadir("oceTAUY_i129_i0_diff.data"),Δdict["oceTAUY"])