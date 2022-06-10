using Revise
using Random
import Pkg
Pkg.activate(".")
using Scalar2D, ReverseDiff
using ADerrors
using Statistics, Plots
using DelimitedFiles
pyplot()

import RecipesBase.plot
plot(f::Function, x::Vector, args...; kwargs...) = plot(x, f.(x), args..., kwargs...)

function print_errors(obs, ID)
    uwobs = ADerrors.uwreal(obs, ID)
    ADerrors.uwerr(uwobs)
    println("Value  =   ", uwobs)
    println("τᵢ     =   ", ADerrors.taui(uwobs,ID), " +/- ", ADerrors.dtaui(uwobs,ID))
    # iw = window(uwobs, ID)
    # r = rho(uwobs, ID)
    # dr = drho(uwobs, ID)
    # println("Window: ", iw)
    # plot(collect(1:iw), r[1:iw], yerr = dr[1:iw], seriestype = :scatter)
end

# Theory Parameters
lsize_1 = 6
lsize_2 = lsize_1
lambda  = 0.5
beta    = 0.537
# msq     = 2 * ((1 - 2*lambda) / beta - 2)
# lam     = 2*lambda/beta^2
# lam     = 1.18
# prm     = LattParmA((lsize_1, lsize_2), msq, lam)
prm     = LattParmB((lsize_1, lsize_2), beta, lambda)

# HMC parameters
tau     = 1.0
nsteps  = 10
epsilon = tau/nsteps
n_traj  = 10

# Initialize phi field with random numbers
phi = zeros(prm.iL[1], prm.iL[2])
randn!(phi)

# Precomile the tape. Only has to be done once, at the beginning of the program.
cf_tape = ReverseDiff.compile(ReverseDiff.GradientTape(p -> -action(p,prm),
                                                       phi))

# Thermalization
acc = Vector{Int64}()
corrs = Vector{Float64}()
M = Vector{Float64}()
chi = Vector{Float64}()
chi_t = Vector{Float64}()
for i in 1:200000
    @time HMC!(phi, epsilon, nsteps, acc, prm)
    # push!(corrs, correlation_function(phi, (20,0)))
    # push!(M, magnetization(phi))
    push!(chi, susceptibility(phi))
    phi_t = he_flow(phi, prm, prm.iL[1]^2/64)
    push!(chi_t, susceptibility(phi_t))
    # println("$i")
end

tag = "test60"
print_errors(chi[1000:end], tag)
print_errors(chi_t[1000:end], tag)

# Perform n_traj HMC steps
acc = Vector{Int64}()
M = Vector{Float64}()
chi = Vector{Float64}()
corrs = []
# cfgs = reshape(phi, (1,8,8))

for i in 1:100000
    @time HMC!(phi, epsilon, nsteps, acc, prm, integrator = Leapfrog());
    push!(M, magnetization(phi))
    # push!(corrs, correlation_function(phi))
    # push!(chi, susceptibility(phi))
    # if i%100 == 0
    #     cfgs = vcat(cfgs, reshape(phi, (1,8,8)))
    # end
    # println("Config: ", i)
end

avg_M = mean(M)

chi_bad = chi2(M, prm) |> mean

tag = "test44"
print_errors(M[10000:1:end], tag)

print_errors(M[1000:10:end] .|> abs, tag)
println("Acc = ", mean(acc))

plot(1:length(M[1:20000]), M[1:20000])

histogram([phi...], bins=40)

histogram(M[1:1:end], bins=50, normed=true)

writedlm("/home/david/git/dalbandea/phd/codes/3-Phi4/dav-multilevel-flow-experiments/flows/main/b0.7-cfgs-well.txt", reshape(permutedims(cfgs, [3,2,1]), length(cfgs)))

histogram([mean(cfgs, dims=(2,3))...], bins=40, normed=true)

Cs = hcat(corrs...) |> permutedims

Csmean = []
for i in 1:size(Cs)[2]
    push!(Csmean, Cs[1000:end,i] |> mean)
end

wdir = "/home/david/git/dalbandea/phd/codes/3-Phi4/Scalar2D.jl/results/1-Checks/L80_b0.677_l0.5_D2022-05-18-19-36-22/measurements/L80_b0.677_l0.5_mag.txt"

M = [readdlm(wdir)...]

print_errors(M[10:1:end], wdir) 

print_errors(M[10:1:end] .|> abs, wdir) 

corr = correlation_matrix(phi)
heatmap(circshift(corr ./ corr[1], (10,10)), c = cgrad([:darkviolet, :lightblue, :green, :yellow]))

function distance_matrix(phi)
    L = size(phi)[1]
    dist = similar(phi)
    for i in 1:L
        for j in 1:L
            dist[i,j] = sqrt((min(i-1, L-i+1))^2 + (min(j-1, L-j+1))^2)
        end
    end
    return dist
end

dists = distance_matrix(phi)

corr_rsh = reshape(corr, (:))
dists_rsh = reshape(dists, (:))
order = sortperm(dists_rsh)

plot(dists_rsh[order], corr_rsh[order], seriestype = :scatter)

## Gradient flow


M_t = Vector{Float64}()

for t in 0:0.01:10
    phi_t = he_flow(phi, prm, t)
    push!(M_t, magnetization(phi_t))
end
