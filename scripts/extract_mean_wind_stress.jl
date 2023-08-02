#compute the time-mean wind stress from 

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

τx_mean = MeshArray(γ,Float32); fill!(τx_mean, 0.0)
τy_mean = MeshArray(γ,Float32); fill!(τy_mean, 0.0)


expname = "iter129_bulkformula"

filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"


for tt = 1:nt
    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    Tname = τdatafilelist[tt]
    
    τx, τy = extract_ocnTAU(diagpath, expname, Tname, γ)

    for ff = 1:5
        τx_mean.f[ff] .+= τx.f[ff] ./ nt
        τy_mean.f[ff] .+= τy.f[ff] ./ nt
    end
end
