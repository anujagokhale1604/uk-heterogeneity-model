# Two-Agent New Keynesian Model: Heterogeneous Response to Monetary Tightening
# Gokhale (2026) — companion code to "Who Leaves When Rates Rise?"
# github.com/anujagokhale1604

using Parameters, Plots

@with_kw struct Params
    β::Float64 = 0.99
    σ::Float64 = 1.0
    κ::Float64 = 0.1
    ϕ_π::Float64 = 1.5
    ϕ_y::Float64 = 0.5
    ρ_r::Float64 = 0.7
    θ_O::Float64 = 0.6    # mortgage exposure: owners
    θ_R::Float64 = 0.0    # mortgage exposure: renters
    T::Int = 20
end

function simulate(p::Params, shock_size::Float64 = 0.01)
    @unpack β, σ, κ, ϕ_π, ϕ_y, ρ_r, θ_O, θ_R, T = p
    r = zeros(T); π = zeros(T); y = zeros(T)
    c_O = zeros(T); c_R = zeros(T)
    r[1] = shock_size
    for t in 2:T
        r[t] = ρ_r * r[t-1]
        π[t] = β * (t < T ? π[t+1] : 0.0) + κ * y[t]
        y[t] = (t < T ? y[t+1] : 0.0) - (1/σ) * (r[t] - π[t])
        c_O[t] = -(1/σ) * (r[t] + θ_O * r[t])
        c_R[t] = -(1/σ) * (r[t] + θ_R * r[t])
    end
    return r, π, y, c_O, c_R
end

p = Params()
r, π, y, c_O, c_R = simulate(p, 0.01)
t = 1:p.T

# Figure 1: Consumption differential
plt1 = plot(size=(700,400), xlabel="Quarters after shock", ylabel="Consumption response",
    title="Heterogeneous Consumption Response to Monetary Tightening",
    legend=:bottomright, grid=true, gridcolor=:lightgrey, bg=:white, framestyle=:box)
plot!(plt1, t, c_O, color=:crimson, lw=2.5, label="Mortgaged owners (θ = 0.6)")
plot!(plt1, t, c_R, color=:steelblue, lw=2.5, label="Renters (θ = 0.0)")
hline!(plt1, [0], color=:black, lw=0.8, linestyle=:dash, label="")
savefig(plt1, "irf_clean.png")

# Figure 2: Emigration intent gap
emigration_gap = abs.(c_O .- c_R) .* 8
plt2 = plot(size=(750,420), xlabel="Quarters after shock", ylabel="Response (normalised)",
    title="Consumption Contraction and Emigration Intent Gap\nFollowing Monetary Tightening Shock",
    legend=:right, grid=true, gridcolor=:lightgrey, bg=:white, framestyle=:box)
plot!(plt2, t, c_O, color=:crimson, lw=2.5, label="Owner consumption")
plot!(plt2, t, c_R, color=:steelblue, lw=2.5, label="Renter consumption")
plot!(plt2, t, emigration_gap, color=:darkorange, lw=2, linestyle=:dash, label="Emigration intent gap (scaled ×8)")
hline!(plt2, [0], color=:black, lw=0.8, linestyle=:dot, label="")
savefig(plt2, "irf_emigration.png")

println("Figures saved.")
