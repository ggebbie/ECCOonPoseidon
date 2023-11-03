using FFTW, RollingFunctions

nt = 312
freqs = 1 ./ reverse(FFTW.rfftfreq(nt)[2:end]) 
freqs ./= 12
Δt = 1; N = nt; T = N*Δt

window_size = 3
roll_freqs = rollmean(freqs, window_size)


spectral_density(x, N, T) = reverse((2*T * inv(N^2)) .* abs2.(rfft(x .- mean(x))[2:end] ))
rolling_spectral_density(x, N, T) = rollmean(spectral_density(x, N, T), window_size)
perc_var_spectral_dens(x, N, T) = rolling_spectral_density(x, N, T) ./ sum(rolling_spectral_density(x, N, T))