#this script cannot be run from teh repl 
using PyCall
@pyimport matplotlib.animation as anim
using PyPlot, IJulia
# export PATH="/home/ameza/batou-home/anaconda3/pkgs/ffmpeg-5.1.2-gpl_hbd009f3_102/bin/:$PATH"
#need to export this to path^
global fig, ax 
fig, ax = plt.subplots(figsize=(4,4))
x = [0:0.01:2pi;]
# i=0,1,...,frames-1
function animate(i)
    ax.clear()
    ax.plot(x, sin.(x.+i/10.0))
end

myanim = anim.FuncAnimation(fig, animate, frames=200, interval=25, blit=false)
myanim.save("test3.mp4")

