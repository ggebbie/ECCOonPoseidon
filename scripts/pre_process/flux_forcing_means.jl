include("../../src/intro.jl")

using Revise
using ECCOtour
using ECCOonPoseidon
using Statistics
using PythonPlot
using Distributions
using FFTW
using LinearAlgebra
using StatsBase
using MeshArrays
using MITgcmTools
using MAT

expt == "nointerannual" ? keepregion = false : keepregion = true
println("Experiment: ",expt)
include(srcdir("config_exp.jl"))
include(srcdir("config_regularpoles.jl"))

inputdir_exp = fluxdir("nointerannual_Feb5")
inputdir = fluxdir()
outputdir_exp = fluxdir("nointerannual_meantest")

# read lat, lon at center of grid cell
(ϕC,λC) = latlonC(γ)
# on the vector (Staggered) grid
(ϕG,λG) = latlonG(γ)

# set 160 West as the center of the regpoles grid
lonmid =  -160
centerlon!(λC,lonmid)
centerlon!(λG,lonmid)
midname = "_6hourlyavg_"
varnames = ("atmPload","oceFWflx","oceQsw","oceSflux","oceSPflx","oceTAUE","oceTAUN","oceTAUX",
            "oceTAUY","sIceLoad","sIceLoadPatmPload","sIceLoadPatmPload_nopabar","TFLUX")
years = 1992:2017

#calculate mean from original flux files
vname = varnames[2]
#for vname ∈ varnames
    filein = joinpath(inputdir,vname*midname)
    println(filein)
    #initialize matrix to hold added fields
    #force mean fluxes to match original
    filein_exp = joinpath(inputdir_exp,vname*midname)
    #initialize matrix to hold added fields
    meanout = MeshArray(γ,Float32,1) 
    meanout_exp = MeshArray(γ,Float32,1)
    fill!(meanout,0.0)
    fill!(meanout_exp,0.0)

    #add each year of flux data together, then add to meanout
    @time for tt in 1:length(years)
        fname = filein*string(years[tt]) #get filename for this variable and year (original fluxes)
        field = read_bin(fname,Float32,γ); #read in data
        fname_exp = filein_exp*string(years[tt]) #get filename for this variable and year (modified fluxes)
        field_exp = read_bin(fname_exp,Float32,γ); #read in data
        for fv in 1:5 #for each of the five grid faces
            tmplat = ϕC[fv]; #get lat grid
            tmplon = λC[fv]; #get lon grid
            l1 = size(meanout[fv],1);
            l2 = size(meanout[fv],2);
            tmptot = fill(0.0, l1, l2); #array of zeros to hold added fluxes for this grid face in this year
            tmptot_exp = fill(0.0, l1, l2);
            for ss in 1:size(field,2) #for each 6 hour timestep in this year
                field_face = field[fv,ss]; #lon/lat grid of flux for this timestep
                tmptot .+= field_face; 
                field_face_exp = field_exp[fv,ss];
                tmptot_exp .+= field_face_exp;
            end
            meanout.f[fv] .+= tmptot;
            meanout_exp.f[fv] .+= tmptot_exp;
        end
    end
    meanout = meanout./37987; #divide by total number of 6 hour timesteps in 26 year timeseries
    meanout_exp = meanout_exp./37987; #divide by total number of 6 hour timesteps in 26 year timeseries
    mean_diff = meanout_exp .- meanout;

    mean_diff_regpoles = regularpoles(mean_diff,γ,rp_params)
    meanout_regpoles = regularpoles(meanout,γ,rp_params)
    meanout_exp_regpoles = regularpoles(meanout_exp,γ,rp_params)

    figure()
    clf()
    cmap_seismic =get_cmap("seismic")
    lims = range(-1.0,step=0.05,stop=1.0)
    contourf(rp_params.λC,
    rp_params.ϕC,
    mean_diff_regpoles',
    lims,
    cmap=cmap_seismic)
    colorbar(label="weight",orientation="vertical",ticks=lims)
    !ispath(plotsdir()) && mkpath(plotsdir())
    outfname = plotsdir("meandiff.pdf")
    xlbl = "longitude "*L"[\degree E]"
    ylbl = "latitude "*L"[\degree N]"
    xlabel(xlbl)
    ylabel(ylbl)
    savefig(outfname)

    figure()
    clf()
    cmap_seismic =get_cmap("seismic")
    lims = range(-200.0,step=10.0,stop=200.0)
    contourf(rp_params.λC,
    rp_params.ϕC,
    meanout_regpoles',
    lims,
    cmap=cmap_seismic)
    colorbar(label="weight",orientation="vertical",ticks=lims)
    !ispath(plotsdir()) && mkpath(plotsdir())
    outfname = plotsdir("meanout.pdf")
    xlbl = "longitude "*L"[\degree E]"
    ylbl = "latitude "*L"[\degree N]"
    xlabel(xlbl)
    ylabel(ylbl)
    savefig(outfname)

    figure()
    clf()
    cmap_seismic =get_cmap("seismic")
    lims = range(-200.0,step=10.0,stop=200.0)
    contourf(rp_params.λC,
    rp_params.ϕC,
    meanout_exp_regpoles',
    lims,
    cmap=cmap_seismic)
    colorbar(label="weight",orientation="vertical",ticks=lims)
    !ispath(plotsdir()) && mkpath(plotsdir())
    outfname = plotsdir("meanout_exp.pdf")
    xlbl = "longitude "*L"[\degree E]"
    ylbl = "latitude "*L"[\degree N]"
    xlabel(xlbl)
    ylabel(ylbl)
    savefig(outfname)

    #now load in each experiment flux file and subtract mean_diff, then re-save updated flux values
    for tt in 1:length(years)
        fname_exp = filein_exp*string(years[tt])
        field_exp = read_bin(fname_exp,Float32,γ);
        field_exp_fixmean = field_exp .- mean_diff;
        #save new files
        fileout = joinpath(outputdir_exp,vname*midname)
        fnameout = fileout*string(years[tt])
        println("saving file "*fnameout)
        write(fnameout,field_exp_fixmean)
    end
#end


#bad individual timeseries method ----------------------------
#for vname = varnames[end] #∈ varnames
#    filein = joinpath(inputdir,vname*midname)
#    println(filein)
#    frootsample = inputdir*vname*midname
#   #initialize matrix to hold added fields
#    meanout = MeshArray(γ,Float32,1)
#    fill!(meanout,NaN)
#    for fv = 1 #in 1:5
#        tmplat = ϕC[fv];
#        tmplon = λC[fv];
#        l1 = size(meanout[fv],1);
#        l2 = size(meanout[fv],2);
#        tmpmean = fill(NaN, l1, l2)
#        for yv = 1 #in 1:size(tmplat,2)
#            for xv = 1 #in size(tmplat,1)
#                @time fluxsample_point,nseries = extract_timeseries(filein,years,γ,xv,yv,fv)
#                tmpmean[xv,yv] = mean(fluxsample_point)
#            end
#        end
#        meanout[fv] = tmpmean;
#    end
#end


#force mean fluxes to match original - this would only get rid of the small difference from the filtering step -----------
#meanout = MeshArray(γ,Float32,1); #initialize array
#fill!(meanout,0.0);

#calculate mean at each point of flux_14day_lopass -- should be zero but it's not, which is the whole problem
#@time for fv in 1:5 #for each of the five grid faces
#        tmplat = ϕC[fv];
#        tmplon = λC[fv];
#        l1 = size(flux_14day_lopass[fv],1);
#        l2 = size(flux_14day_lopass[fv],2);
#        tmptot = fill(0.0, l1, l2); #empty array to hold added fluxes for this grid face
#        for ss in 1:size(flux_14day_lopass,2) #for each 14 day timestep
#            field_face = flux_14day_lopass[fv,ss]; #lon/lat grid of flux for this timestep
#           tmptot = tmptot .+ field_face;
#        end
#       meanout[fv] = meanout[fv] .+ tmptot;
#end
#meanout = meanout./size(flux_14day_lopass,2); #divide by total number of 14 day timesteps in 26 year timeseries

#then subtract mean of flux_14day_lopass from flux_14day_lopass
#flux_14day_lopass_meanfix = flux_14day_lopass .- meanout;
