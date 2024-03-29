# Work in progress.
include("../../src/intro.jl")

using Revise
using MITgcmTools
using MeshArrays
using ECCOonPoseidon
using ECCOtour
using GeoPythonPlot
using LaTeXStrings
using JLD2
import GeoPythonPlot as GPP

# recompute or read from file?
readfromfile = true

# list of experiments on poseidon
runpath,diagpath = listexperiments(exprootdir())

# abbreviations for each experiment for labels, etc.
shortnames = expnames()

## SELECT EXPERIMENTS TO COMPARE #################################
# baseline experiment
expbase = "iter129_bulkformula"

# comparison experiment(s)
#expcompare = "iter129_fluxforced"
#expcompare = "nosfcadjust"
expcompare = "noinitadjust"

if readfromfile
    outfile = datadir("divergence","faststats_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*".jld2")
    @load outfile xbar σx xmax xmin absxbar z 
else
    # or re-run diagnostics
    include("experiment_divergence.jl")
end

# Look at both min and max to get the most extreme value.
xextreme = max.(-xmin,xmax);

!ispath(plotsdir("divergence")) && mkpath(plotsdir("divergence"))

# Choose the depth level and variable of interest.
zlev = 20  # there are nz=50 levels for ECCOv4r4.
zlbl  = "depth [m]"
depthlbl = string(round(Int,z[zlev]))*"m"
for cval = 1:2 # variable 1 = theta, variable 2 = practical salinity

    ################ Linear scaling of divergence ##################
    xlbl = "month starting in Jan 1992"
    if cval == 1
        str1 = raw"$\theta_{"*shortnames[expbase]*raw"}$"
        str2 = raw"$\theta_{"*shortnames[expcompare]*raw"}$"
        titlelbl = str1*raw"$-$"*str2*","*depthlbl
        linfname = plotsdir("divergence","dtheta_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_1992-2017_"*depthlbl*".pdf")
        logfname = plotsdir("divergence","logdtheta_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_1992-2017_"*depthlbl*".pdf")
        ylbl = L"\theta"*" "*L"[^{\degree}C]"
    elseif cval == 2
        str1 = raw"$S_{"*shortnames[expbase]*raw"}$"
        str2 = raw"$S_{"*shortnames[expcompare]*raw"}$"
        titlelbl = str1*raw"$-$"*str2*","*depthlbl
        linfname = plotsdir("divergence","dsalt_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_1992-2017_"*depthlbl*".pdf")
        logfname = plotsdir("divergence","logdsalt_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_1992-2017_"*depthlbl*".pdf")
        ylbl = "salinity  [PSS-1978]"
    end

    GPP.figure(10+cval)
    GPP.clf()
    GPP.title(titlelbl)
    lev = zlev + (cval-1)*50 # 2 properties in 1 variable
    GPP.plot(σx[:,lev],label=L"std()")
    GPP.plot(absxbar[:,lev],label=L"mean(||)")
    GPP.plot(xextreme[:,lev],label=L"max(||)")
    GPP.xlabel(xlbl)
    GPP.ylabel(ylbl)
    GPP.grid(true)
    GPP.legend()
    GPP.savefig(linfname)

    GPP.figure(20+cval)
    GPP.clf()
    GPP.title(titlelbl)
    lev = zlev + (cval-1)*50
    GPP.semilogy(σx[:,lev],label=L"std()")
    GPP.semilogy(absxbar[:,lev],label=L"mean(||)")
    GPP.semilogy(xextreme[:,lev],label=L"max(||)")
    GPP.xlabel(xlbl)
    GPP.ylabel(ylbl)
    GPP.grid(true)
    GPP.legend()
    GPP.savefig(logfname)

    ## make divergence at a snapshot as a function of depth
    #  pick a time index: here, beginning and end
    local nt = size(xbar,1) # where can be 
    tlist = (1,nt)
    for tt ∈ tlist
        tlbl = time_label(tt-1) # subtract one, months since Jan 1992
        if cval == 1
            timefname = plotsdir("divergence","dtheta_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_"*tlbl[1:3]*tlbl[5:8]*".pdf")
        elseif cval ==2
            
            timefname = plotsdir("divergence","dsalt_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_"*tlbl[1:3]*tlbl[5:8]*".pdf")
        end
        
        levs = 1+(cval-1)*50:cval*50
            
        GPP.figure(12+tt+cval)
        GPP.clf()
        GPP.semilogx(σx[tt,levs],z,label=L"std()")
        GPP.semilogx(absxbar[tt,levs],z,label=L"mean(||)")
        GPP.semilogx(xextreme[tt,levs],z,label=L"max(||)")
        GPP.title(str1*raw"$-$"*str2*", "*tlbl)
        GPP.xlabel(ylbl)
        GPP.ylabel(zlbl)
        GPP.grid(true)
        GPP.legend()
        GPP.gca().invert_yaxis()
        GPP.savefig(timefname)
    end
end

# Ratio of divergence 
# zlev = 20
# #figure()
# #plot(σ[:,20,1])
# title(L"(\theta_{129} - \theta_{ff})/(\theta_{129} - \theta_0) "*string(floor(Int,-Z[zlev]))*"m")
# plot(σstdff[:,zlev,1]./σstd0[:,zlev,1],label=L"std")
# plot(σmedff[:,zlev,1]./σmed0[:,zlev,1],label=L"median")
# plot(σmeanff[:,zlev,1]./σmean0[:,zlev,1],label=L"mean")
# #plot(σmax[:,zlev,1],label=L"max(|\theta_{129}-\theta_{0}|)")
# xlabel("month starting in Jan. 1992")
# ylabel(L"\theta  [^{\circ}C]")
# grid(true)
# legend()
# savefig("dtheta_129v0vff_1992-2017_299m_nomax_12feb2021.pdf")


# plot a timeseries of T vs T and S vs S would be interesting. 

# make a timeseries of the statistics. would be cool.

# make one plan view.
