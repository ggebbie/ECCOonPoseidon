include("../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools,
    PyPlot, JLD2, DrWatson

include(srcdir("config_exp.jl"))

# ϕ = lat, λ = lon
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());

# abbreviations for each experiment for labels, etc.
shortnames = expnames()
marks = expsymbols()
nexp = length(shortnames) # number of experiments
 
fileroot = "state_3d_set1"
dryval = 0.0
iswet(x) = x != dryval

# transfer to nino34 by removing seasonal climatology. Read climatology.
#inputpath = "../inputs/"
#@load inputpath*"nino34_1870-2020.jld2" nino34 nino34yr SSTclimatology
#nino34hadisst = nino34; thadisst = nino34yr; sst34hadisst = SSTclimatology;
nino34hadisst,thadisst,sst34hadisst = ECCOonPoseidon.historicalNino34((1960,2020))


# pre-define
nino34 = Dict{String,Array{Any,1}}() # don't forget trailing parentheses
nino34native = Dict{String,Array{Any,1}}() 
# to do: remove seasonal climatology from model or from reality

tecco = 1992+1/24:1/12:2018
monthsperyear = 12
fcycle = 1 # units: yr^{-1}
# for removing seasonal cycle from monthly averaged timeseries
overtones= 4; # semi-annual, quad-annual, etc.: helps capture an asymmetric seasonal cycle
Ecycle,Fcycle = seasonal_matrices(fcycle,tecco,overtones)

for (keys,values) in shortnames

    expname = keys
    println(expname)
    sst34 = extract_sst34(expname,diagpath,fileroot,γ,area,ϕ,λ,iswet)

    # remove climatology
    nino34[expname] = remove_climatology(sst34,sst34hadisst)
    nino34native[expname] = remove_seasonal(sst34,Ecycle,Fcycle)

end

# make figure with labels, legend
# get nino34 with mean of MeshArray, a patch function for this rectangle, area weighting, could involve updating "mean" function in Reemergence/MeshArrays 
ylbl  = "NINO3.4 "*L"[\degree C]"
figure(100); clf();
plot(thadisst,nino34hadisst,label="HADISST 1")
for (keys,values) in nino34
    plot(tecco,values,"-"*marks[keys],label = shortnames[keys])
end
grid(true)
legend()
ylabel(ylbl)
xlabel("calendar years")
axorig = axis()
axis((1992,2018,axorig[3],axorig[4]))
outputfile = plotsdir("nino34comparison_sameSSTscale.eps")

isdir(plotsdir()) ? nothing : mkdir(plotsdir())

savefig(outputfile)

## make second figure with different SST baseline
ylbl  = "NINO3.4 relative to modeled SST "*L"[\degree C]"
figure(101); clf();
plot(thadisst,nino34hadisst,label="HADISST 1")
for (keys,values) in nino34native
    plot(tecco,values,"-"*marks[keys],label = shortnames[keys])
end
grid(true)
legend()
ylabel(ylbl)
xlabel("calendar years")
axorig = axis()
axis((1992,2018,axorig[3],axorig[4]))
savefig(plotsdir("nino34comparison_nativeSSTscale.eps"))
