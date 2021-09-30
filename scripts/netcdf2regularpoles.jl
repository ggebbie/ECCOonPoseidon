# Special case: Take EVEL* NVEL* from nctiles_monthly and map onto regular grid. 

# 1. read monthly-average fields
# 2. interpolate to Cartesian grid
# 3. save to self-describing file.
# 4. repeat

include("../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour, MeshArrays, MITgcmTools
using NetCDF

#exppath = "/batou/ECCOv4r4/MITgcm/exps/"
#runpath,diagpath,regpolespath = listexperiments(exppath);
varnames = ("EVEL","EVELMASS","EVELSTAR","NVEL","NVELMASS","NVELSTAR","WVELMASS","WVELSTAR","oceTAUN","oceTAUE")

include(srcdir("config_exp.jl"))
include(srcdir("config_regularpoles.jl"))

# Read nctiles_monthly from ECCO Drive
# this kludge needs to be replaced
exppath = "/batou/eccodrive/nctiles_monthly/"

# varname = varnames[1] # for interactive use
for varname in varnames

    println("varname ",varname)
    exppathvar = exppath*varname*"/"

    pathoutvar = pathout*varname*"/"
    !isdir(pathoutvar) ? mkdir(pathoutvar) : nothing;

    # make directory if doesn't exist
    for tt = 1:length(tecco)

        year,month = timestamp_monthly_v4r4(tt)

        if month < 10
            filesuffix = "_"*string(year)*"_0"*string(month)*".nc"
            nctilename = exppathvar*string(year)*"/"*varname*filesuffix
        else
            filesuffix = "_"*string(year)*"_"*string(month)*".nc"
            nctilename = exppathvar*string(year)*"/"*varname*filesuffix
        end
        
        @time varsregpoles =  netcdf2regularpoles(nctilename,varname,γ,nx,ny,nyarc,farc,iarc,jarc,warc,nyantarc,fantarc,iantarc,jantarc,wantarc)

        @time writeregularpoles(varsregpoles,γ,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts,z,depthatts)
        
    end
end
