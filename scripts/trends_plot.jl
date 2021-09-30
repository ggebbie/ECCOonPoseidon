#  First steps: 1. go into Reemergence project directory. 2. go into julia REPL package mode with `]`. 3. `activate .` 4. Backspace to return to command mode in REPL.

include("../src/intro.jl")

# make spatial plots of trends computed from trends.jl.

using Revise
using MeshArrays, MITgcmTools, ECCOonPoseidon, ECCOtour
using Statistics, PyPlot, Distributions
using LinearAlgebra, StatsBase

include(srcdir("config_exp.jl"))

include(srcdir("config_regularpoles.jl"))
        
# list of experiments on poseidon
runpath,diagpath = listexperiments(exprootdir())

# abbreviations for each experiment for labels, etc.
shortnames = expnames()

## SELECT EXPERIMENTS TO ANALYZE #################################
# manually choose from available experiments listed above.
#exps = ("iter129_bulkformula","nointerannual")

# to do all experiments:
exps = keys(shortnames)
#################################################################

nexps = length(exps) # number of experiments

# assume monthly averages, ECCOv4r4
tstart = 1992 + 1/24
tend = 2018
tecco = range(tstart,step=1/12,stop=2018)
nt = length(tecco)
nx = length(λC); ny = length(ϕC)
        
cmap_seismic = get_cmap("seismic")

for exp in exps

    # read trends from file. 
    outdir = datadir("trends",exp)
    Tname = datadir(outdir,"DthetaDt.data")
    @time β = γ.read(Tname,MeshArray(γ,Float32,nz))

    figure(101)
    for zz = 1:nz
        βz = 1.0f4*β[:,zz] # units K/yr -> cK/century

        # problem here?
        βzreg = var2regularpoles(βz,γ,nx,ny,nyarc,farc,iarc,jarc,warc,nyantarc,fantarc,iantarc,jantarc,wantarc)
        #βzreg = reginterp(βz,nx,ny,f,i,j,w)
        mx = maximum(βz,NaN32) # filter out NaN32
        mn = minimum(βz,NaN32) 
        extrm = max(mx,-mn)
        clf()
        lims = range(-extrm,step=2extrm/30,stop=extrm)
        contourf(longrid,latgrid,βzreg,lims,cmap=cmap_seismic)
        colorbar(label="cK/century",orientation="vertical",ticks=lims)
        contour(longrid,latgrid,βzreg,lims,colors="k")

        depthlbl = string(round(Int,-z[zz]))*" m"
        depthlbl2 = string(round(Int,-z[zz]))*"m"
        titlelbl = exp*", "*" "*depthlbl
        outfname = plotsdir(exp,"DthetaDt_"*shortnames[exp]*"_"*depthlbl2*".eps")
        xlbl = "longitude "*L"[\degree E]"
        ylbl = "latitude "*L"[\degree N]"
        title(titlelbl)
        xlabel(xlbl)
        ylabel(ylbl)
        savefig(outfname)
    end
end
