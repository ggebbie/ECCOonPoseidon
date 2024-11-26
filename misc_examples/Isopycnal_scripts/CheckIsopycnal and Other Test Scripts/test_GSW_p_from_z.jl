using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings,
Plots
using Interpolations
import NaNMath as nm
using GibbsSeaWater

p = [   10;   50;  125;  250;  600; 1000;2000; 3000;]
# p =  [   10;   50;  125;  250;  600; 1000;]

z1 = gsw_z_from_p.(p,4, 0, 0)
z2 = gsw_z_from_p.(p,-80, 0, 0)

plot(p, z1)
plot!(p, z2)