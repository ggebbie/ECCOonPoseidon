include("../src/intro.jl")
using Revise 
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, Statistics, JLD2, DSP, 
    FFTW
import PyPlot as plt 

fc1 = 5; fc2 = 15; fs = 1000
fil = digitalfilter(Bandpass(fc1,fc2,fs=fs), Butterworth(5));

t = collect(0:1/1000:1) # 1 second
sig = sin.((2*π*10).*t) .+ sin.((2*π*20).*t)
bp_sig = filt(fil,sig)

fig, axs = plt.subplots(2, 1, figsize=(5,5), sharey = true)
axs[1].plot(t, sig)
axs[2].plot(t, bp_sig, label = "filtered signal")
axs[2].plot(t,  sin.((2*π*10).*t), label = "Signal I want to isolate")
axs[1].legend()
fig