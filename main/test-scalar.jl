using Revise
using Random
import Pkg
Pkg.activate(".")
using Scalar2D, ReverseDiff

# Theory Parameters
lsize_1 = 32
lsize_2 = lsize_1
lambda  = 0.02
beta    = 0.6
msq     = 2 * ((1 - 2*lambda) / beta - 2)
lam     = 2*lambda/beta^2
prm     = LattParmA((lsize_1, lsize_2), msq, lam)
# prm     = LattParmB((lsize_1, lsize_2), beta, lambda)

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
for i in 1:10
    @time HMC!(phi, epsilon, nsteps, acc, prm)
end

# Perform n_traj HMC steps
acc = Vector{Int64}()
M = Vector{Float64}()
chi = Vector{Float64}()
for i in 1:10
    @time HMC!(phi, epsilon, nsteps, acc, prm, integrator = OMF4())
    push!(M, magnetization(phi))
    push!(chi, susceptibility(phi))
end

avg_M = mean(M)

chi_bad = chi2(M, prm) |> mean
