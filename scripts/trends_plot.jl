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

# assume monthly averages, ECCOv4r4
tstart = 1992 + 1/24
tend = 2018
tecco = range(tstart,step=1/12,stop=2018)
nt = length(tecco)
nx = length(λC); ny = length(ϕC)
        
cmap_seismic = get_cmap("seismic")

for expt in exps

    # read trends from file. 
    outdir = datadir("trends",expt)
    Tname = datadir(outdir,"DthetaDt.data")
    @time β = γ.read(Tname,MeshArray(γ,Float32,nz))

    fig = figure(101)
    for zz = 1:nz
        βz = 1.0f4*β[:,zz] # units K/yr -> cK/century

        #βzreg = var2regularpoles(βz,γ,nx,ny,nyarc,farc,iarc,jarc,warc,nyantarc,fantarc,iantarc,jantarc,wantarc)
        βzreg = var2regularpoles(βz,γ,nx,ny,nyarc,λarc,nyantarc,λantarc)
        
        mx = maximum(MeshArrays.mask(βz,-Inf))
        mn = minimum(MeshArrays.mask(βz,Inf))

        extrm = max(mx,-mn)

        clf()
        lims = range(-extrm,step=2extrm/30,stop=extrm)

        cenlon = -160.0
        proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree()
        proj = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=cenlon)
        ax = fig.add_subplot(projection = proj)
        ax.set_global()
        ax.coastlines()
        
        contourf(λC,ϕC,βzreg',cmap=cmap_seismic, transform = proj0)
        # contourf(λC,ϕC,βzreg',lims,cmap=cmap_seismic)
        colorbar(label="cK/century",orientation="vertical",ticks=lims)
        contour(λC,ϕC,βzreg',colors="k",transform = proj0)
        # contour(λC,ϕC,βzreg',lims,colors="k")

        depthlbl = string(round(Int,z[zz]))*" m"
        depthlbl2 = string(round(Int,z[zz]))*"m"
        titlelbl = expt*", "*" "*depthlbl
        outdir = plotsdir(expt)
        !isdir(outdir) ? mkpath(outdir) : nothing
        outfname = plotsdir(expt,"DthetaDt_"*shortnames[expt]*"_"*depthlbl2*".pdf")
        xlbl = "longitude "*L"[\degree E]"
        ylbl = "latitude "*L"[\degree N]"
        ax.set_title(titlelbl)
        ax.set_xlabel(xlbl)
        ax.set_ylabel(ylbl)
        savefig(outfname)
    end
end
