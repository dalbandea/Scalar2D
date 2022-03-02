using Random
import Pkg
Pkg.activate(".")
using Revise, Scalar2D, ReverseDiff

# Theory Parameters
lsize_1 = 24
lsize_2 = 24
msq     = 0.1
lam     = 0.5
beta    = 2.0
prm     = LattParm((lsize_1, lsize_2), msq, lam, beta)

# HMC parameters
tau     = 1.0
nsteps  = 100
epsilon = tau/nsteps
n_traj  = 10


# Initialize phi field with random numbers
phi = zeros(prm.iL[1], prm.iL[2])
randn!(phi)

# Precomile the tape. Only has to be done once, at the beginning of the program.
cf_tape = ReverseDiff.compile(ReverseDiff.GradientTape(p -> -action(p,prm),
                                                       phi))

# Perform n_traj HMC steps
acc = Vector{Int64}()
for i in 1:n_traj
    @time HMC!(phi, epsilon, nsteps, acc, prm, AD_tape = cf_tape)
end
