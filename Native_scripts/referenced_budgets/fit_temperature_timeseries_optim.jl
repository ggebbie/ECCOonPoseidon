include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise, ECCOonPoseidon, ECCOtour, JLD2, 
    DrWatson, PyCall, MITgcmTools, 
    MeshArrays, Distributions, Optim, DataFrames, 
    Measurements

using .OHC_helper
import PyPlot as plt
@pyimport seaborn as sns;
@pyimport pandas as pd;

include(srcdir("config_exp.jl"))
include("optim_funcs.jl")

sns.set_theme(context = "poster", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));

tecco = 1992+1/24:1/12:2018
region = "NPAC"; suffix = "2to3"

fig,axs=plt.subplots(1,1, sharey = true, figsize = (10, 7.5))

fname(x) = datadir("native/" * x * region * "_THETA_budget_ref_" * suffix * ".jld2")
load_θ(expname) = load(fname(expname))["dθ"]["θ"]

θ_iter129 = load_θ("iter129_bulkformula"); θ_iter0 = load_θ("noinitadjust");

x0 = [0.1, 0.1, 0.001, 0.001, 1.5]; #initial ;guess
t = (collect(tecco) .- tecco[1]);

f(x) = loss_function(θ_iter129, θ_iter0, t, x);
sol = optimize(f, x0, NelderMead());
a1, a2, b1, b2, c = Optim.minimizer(sol);
println(Optim.minimizer(sol))

ỹiter129 = obj_func(t, [a1, b1, c]); 
ỹiter0 = obj_func(t, [a2, b2, c])

ỹiter129_sim, ỹiter0_sim, x_ens = error_bars(θ_iter129, ỹiter129, θ_iter0, ỹiter0, t, t);
x_mean, x_std = mean_std(x_ens);

#PLOTTING
fig,axs=plt.subplots(1,1, sharey = true, figsize = (10, 7.5));
axs.plot(t, θ_iter129, label = "Iter129");
axs.plot(t, θ_iter0, label = "Iter0");
plot_mean_and_error!(axs, t, ỹiter129_sim);
plot_mean_and_error!(axs, t, ỹiter0_sim);
fig

a1, a2, b1, b2, c = measurement.(x_mean, x_std)
df = DataFrame(a=zeros(Measurement{Float64}, 2), 
b=zeros(Measurement{Float64}, 2), c=zeros(Measurement{Float64}, 2), 
names = ["i129", "i0"])

df.a .= (a1, a2); df.b .= (b1, b2); df.c .= c

# #RUN FORWARD
# t_ext = 0:1000
# ỹiter129_sim, ỹiter0_sim = error_bars(θ_iter129, ỹiter129, θ_iter0, ỹiter0, t, t_ext);

# fig,axs=plt.subplots(1,1, sharey = true, figsize = (10, 7.5));
# axs.plot(t, θ_iter129, label = "Iter129");
# axs.plot(t, θ_iter0, label = "Iter0");
# plot_mean_and_error!(axs, t_ext, ỹiter129_sim);
# plot_mean_and_error!(axs, t_ext, ỹiter0_sim);
# fig
