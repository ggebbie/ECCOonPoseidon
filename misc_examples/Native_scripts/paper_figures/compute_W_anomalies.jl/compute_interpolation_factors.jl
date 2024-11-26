include("../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
using NaNMath
import PyPlot as plt

include(srcdir("plot_and_dir_config.jl"))

lon=[i for i=-179.:0.5:179., j=-89.:0.5:89.]
lat=[j for i=-179.:0.5:179., j=-89.:0.5:89.]
(f,i,j,w)=InterpolationFactors(Γ,vec(lon),vec(lat))

jldsave(datadir("0.5deg_interpolation_factors.jld2"), lon=lon, lat=lat, f=f, i=i, j=j, w=w)
