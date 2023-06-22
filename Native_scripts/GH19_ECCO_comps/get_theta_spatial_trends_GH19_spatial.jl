using Revise, DrWatson, Statistics,
NCDatasets, Printf, 
DataFrames, LaTeXStrings, Distances, ECCOonPoseidon, ECCOtour, JLD2
include("GH19_helperfuncs.jl")

ds_OPT = NCDataset("/home/ameza/GH19.jl/data/Theta_OPT-0015.nc")

year  = reverse(ds_OPT["year"][:]); 
lon = ds_OPT["longitude"][:]
lat = ds_OPT["latitude"][:]
depth = ds_OPT["depth"][:]; nz = length(depth)
zlevs = findall(2000 .<= depth .<= 3000)

nt = length(year)
WOCE_times = findall(1990 .< year .< 2000)[1]
Challenger_times = findall(2016 .<= year )[end]

theta_OPT = reverse(ds_OPT["theta"][:, zlevs, :, :], dims =1) #reverse time for niceness
theta_OPT = theta_OPT[WOCE_times:Challenger_times, :, :, :]
year_CH_WC = year[WOCE_times:Challenger_times]
E,F = trend_matrices(year_CH_WC); nt = length(year_CH_WC)

nt = length(year_CH_WC); nz = length(zlevs)
weighted_temp = zeros(nt, 90, 180)
#volume weighted average 
for tt = 1:nt, k in 1:nz
    weighted_temp[tt, :, :] .+=  theta_OPT[tt, k, :, :] ./ nz 
end

β = zeros(90, 180)

for tt = 1:nt
    β .+= F[2,tt] .* weighted_temp[tt, :, :]
end

#coordinate meshgrid
LONS = lon' .* ones(length(lat))
LATS = lat .* ones(length(lon))'

fname = "modern_OPT-0015_θ_trends_2to3km.jld2"
jldsave(datadir(fname), β = β, λ = LONS, ϕ = LATS, years = year_CH_WC)
