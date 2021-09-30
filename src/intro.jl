using DrWatson
@quickactivate "ECCOonPoseidon"

println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())
"""
)

## SELECT EXPERIMENTS TO COMPARE #################################
# manually choose from available experiment
if length(ARGS) > 0
    expt = ARGS[1]
else
    # for interactive use, ARGS may be set this way:
    # use nointerannual as a default
    # following filter_interannual
    #push!(empty!(ARGS), "nointerannual")
    expt = "nointerannual"
    println("No experiment chosen")
    println("Default experiment is ",expt)
end
