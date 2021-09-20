# print output here
path_out = sig1dir(expt)
diagpath = diagdir(expt)
pathout = regpolesdir(expt)
!isdir(pathout) ? mkdir(pathout) : nothing;

# DEFINE THE LIST OF SIGMA1 VALUES.
sig1grid = sigma1grid()

################################################################
# get MITgcm / ECCOv4r4 LLC grid and depth information. Store in γ.
γ = setupLLCgrid(datadir("grid/"))
nf = length(γ.fSize)

# get standard levels of MITgcm
z = depthlevels(γ)
pstdz = pressurelevels(z)
p₀ = 1000.0 ; # dbar
