module Scalar2D

import Random, ReverseDiff, FFTW

abstract type LattParm end

struct LattParmA <: LattParm
    iL::Tuple{Int64,Int64}
    msq::Float64
    lam::Float64
end

struct LattParmB <: LattParm
    iL::Tuple{Int64,Int64}
    beta::Float64
    lambda::Float64
end

export LattParmA, LattParmB

include("integrators/integrators.jl")
export Integrators, Leapfrog, OMF4

include("Scalar2DAction.jl")
export action, force!

include("Scalar2DHMC.jl")
export HMC!, Hamiltonian, OMF4!, leapfrog!, update_momenta!, update_field!

include("Scalar2DMeasurements.jl")
export magnetization, susceptibility, chi2, correlation_function, correlation_function2, G0, correlation_matrix

include("Scalar2DGradientFlow.jl")
export he_flow

end # module
