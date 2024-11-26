Ω = (2π * 86400) / (60 * 60 *24) #cycles per day
f(ϕ) = (2 * Ω * sind(ϕ)) #f in 1/day
β(ϕ) = 2 * Ω * cosd(ϕ) / 6378 #beta in 1/(kilometers * day)

N0 = 7.3e-3 * 86400 #cps -> cpd 
b = 1.3

function trapping_depths(k, ω)
   a = (k^2 + β(45)*k*inv(ω)) * inv(f(45)^2 - ω^2)
   print(a)
   coef1 = sqrt(a) * b * N0; coef2 = sqrt(a) * b * N0 - 1
   return b* log(abs(coef1 / coef2))
end

cpd = collect(0:0.01:0.9)[2:end]
wavenumebers = collect(0.001:0.003:0.01)
CPD, WAVENUMBERS = cpd .* ones(length(wavenumebers))', wavenumebers' .* ones(length(cpd))
import PyPlot as plt
fig, ax = plt.subplots()
cf = ax.contour(WAVENUMBERS, CPD, trapping_depths.(WAVENUMBERS, CPD),colors="black", levels = [1, 5])
ax.clabel(cf, inline=true, fontsize=8)
ax.set_ylabel("Frequency [cpd]")
ax.set_xlabel("Wavenumber [1/km]")
fig