using Revise
using Random
import Pkg
Pkg.activate(".")
using Scalar2D, ReverseDiff
using ADerrors
using Statistics, Plots
using DelimitedFiles
pyplot()

function print_errors(obs, ID)
    uwobs = ADerrors.uwreal(obs, ID)
    ADerrors.uwerr(uwobs)
    println("Value  =   ", uwobs)
    println("τᵢ     =   ", ADerrors.taui(uwobs,ID), " +/- ", ADerrors.dtaui(uwobs,ID))
    # iw = window(uwobs, ID)
    # r = rho(uwobs, ID)
    # dr = drho(uwobs, ID)
    # println("Window: ", iw)
    # plot(collect(1:500), r[1:500], yerr = dr[1:500], seriestype = :scatter)
end

# Theory Parameters
lsize_1 = 20
lsize_2 = lsize_1
lambda  = 0.5
beta    = 0.645
msq     = 2 * ((1 - 2*lambda) / beta - 2)
lam     = 2*lambda/beta^2
# lam     = 1.18
# prm     = LattParmA((lsize_1, lsize_2), msq, lam)
prm     = LattParmB((lsize_1, lsize_2), beta, lambda)

# HMC parameters
tau     = 1.0
nsteps  = 8
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
for i in 1:10
    @time HMC!(phi, epsilon, nsteps, acc, prm)
end

# Perform n_traj HMC steps
acc = Vector{Int64}()
M = Vector{Float64}()
chi = Vector{Float64}()
# cfgs = reshape(phi, (1,8,8))
for i in 1:200000
    @time HMC!(phi, epsilon, nsteps, acc, prm, integrator = Leapfrog());
    push!(M, magnetization(phi))
    push!(chi, susceptibility(phi))
    # if i%100 == 0
    #     cfgs = vcat(cfgs, reshape(phi, (1,8,8)))
    # end
end

avg_M = mean(M)

chi_bad = chi2(M, prm) |> mean

tag = "test6"
print_errors(M[1000:1:end], tag)
print_errors(M[1000:1:end] .|> abs, tag)

plot(1:length(M[1:20000]), M[1:20000])

histogram([phi...], bins=40)

histogram(M[1:1:end], bins=50, normed=true)

writedlm("/home/david/git/dalbandea/phd/codes/3-Phi4/dav-multilevel-flow-experiments/flows/main/b0.7-cfgs-well.txt", reshape(permutedims(cfgs, [3,2,1]), length(cfgs)))

histogram([mean(cfgs, dims=(2,3))...], bins=40, normed=true)
