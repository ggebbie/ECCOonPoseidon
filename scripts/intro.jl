using DrWatson
@quickactivate "ECCOonPoseidon"

println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())

Have fun with your new project!

You can help us improve DrWatson by opening
issues on GitHub, submitting feature requests,
or even opening your own Pull Requests!
"""
)

using ECCOtour

println("success, ECCOtour loaded")

if isempty(ARGS)
    println("ARGS empty")
else
    println(ARGS[1])
    println(typeof(ARGS[1]))
end

