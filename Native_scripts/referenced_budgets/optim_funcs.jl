obj_func(t, x) =  (x[1] .* exp.(-t .* x[2])) .+ x[3] 
sum_squares(y, t, x) = mean((y .- obj_func(t, x)).^2)
loss_function(y1, y2, t, x) = sum_squares(y2, t,  x[[1, 3, 5]]) .+ sum_squares(y1, t,  x[[2, 4, 5]])


function error_bars(y1, ỹ1, y2, ỹ2, t, t_est)
    σ1 = std(y1 .- ỹ1)
    σ2 = std(y2 .- ỹ2)

    n1 = Normal(0, σ1); n2 = Normal(0, σ2)
    nt = length(t); nt_est = length(t_est)
    nrealizations = 10000

    ỹ1_sims = zeros(nt_est, nrealizations)
    ỹ2_sims = zeros(nt_est, nrealizations)
    x_ens = zeros(5, nrealizations)
    for i in 1:nrealizations
        x0 = [0.1, 0.1, 0.001, 0.001, 1.5] #initial guess
        y1_samp = y1 .+ rand(n1, nt)
        y2_samp = y2 .+ rand(n2, nt)

        g(x) = loss_function(y1_samp, y2_samp, t, x)
        sol = optimize(g, x0, NelderMead())
        a1, a2, b1, b2, c = Optim.minimizer(sol);

        x_ens[:, i] .= [a1, a2, b1, b2, c]
        ỹ1_sims[:, i] .= obj_func(t_est, [a1, b1, c])
        ỹ2_sims[:, i] .= obj_func(t_est, [a2, b2, c])
    end

    return ỹ1_sims, ỹ2_sims, x_ens
end

mean_std(x) = (mean(x, dims = 2)[:], std(x, dims = 2)[:])
function plot_mean_and_error!(ax, t, y_ens; c = "k" )
    ens_avg, ens_σ = mean_std(y_ens)
    ax.plot(t, ens_avg, color = c)
    ax.fill_between(t, ens_avg .- (2 .*ens_σ), ens_avg .+ (2 .*ens_σ), alpha=0.3, color=c)
end
