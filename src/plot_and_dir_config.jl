@pyimport seaborn as sns
custom_params = Dict("lines.linewidth" => 2.4)
colors = sns.color_palette("colorblind")
colors_muted = sns.color_palette("muted")

colors = vcat(reverse(colors[1:4]), colors[5:end])
sns.set_theme(context = "talk", style = "ticks",
              palette = colors, rc = custom_params);
              
include(srcdir("config_exp.jl"))

proj0 = ECCOonPoseidon.cartopy.crs.Robinson(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

runpath,diagpath = listexperiments(exprootdir()); 

diagpath["only_init"] = vastdiagdir("nosfcadjust_exps", "run_adjust_init")
diagpath["only_kappa"] = vastdiagdir("nosfcadjust_exps", "run_adjust_kappa")
diagpath["only_sfc"] = vastdiagdir("nooceanadjust", "run_noadjusts")
diagpath["only_buoyancy"] = vastdiagdir("nooceanadjust", "run_noadjusts_nowind")
diagpath["only_wind"] = vastdiagdir("nooceanadjust", "run_noadjusts_nobuoyancy")

diagpath["test_iter0_forc"] = vastdiagdir("iter0_climatological_tau", "run_iter0_tau")
diagpath["test_iter129_forc"] = vastdiagdir("iter129_climatological_tau", "run_iter129_tau")

# diagpath["no_kappa"] = vastdiagdir("nooceanadjust", "run_nokappaadjust")
exps_list = ["iter0_bulkformula" , "only_init", "only_kappa", "only_sfc", "iter129_bulkformula", "only_wind", "only_buoyancy"]
lin_exps = ["CTRL", "Initial", "Kappa", "Forcing", "FULL", "Wind", "Buoyancy"]; nexps = length(lin_exps)
plot_labels = ["Iteration 0", "INIT", "MIXING", "FORCING", "Iteration 129", "WIND", "BUOYANCY"]
plot_labels = Dict(exps_list[i] => plot_labels[i] for i in 1:nexps)

plot_labels_effects = ["Iteration 0", "INIT Effect", "MIXING Effect", 
                        "FORCING Effect", "Iteration 129", "WIND Effect", "BUOYANCY Effect"]
plot_labels_effects = Dict(exps_list[i] => plot_labels_effects[i] for i in 1:nexps)

exp_colors = vcat(colors[1:4], ["k"], [colors[4]], [colors_muted[5]]) 
exp_colors = Dict(exps_list[i] => exp_colors[i] for i in 1:nexps)
exp_colors["iter129_bulkformula"] = deepcopy(exp_colors["only_wind"])
exp_colors["only_wind"] = "#934B00"

