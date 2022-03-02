module Scalar2D

import Random, ReverseDiff

struct LattParm
    iL::Tuple{Int64,Int64}
    msq::Float64
    lam::Float64
    beta::Float64
end
export LattParm

include("Scalar2DAction.jl")
export action, force!

include("Scalar2DHMC.jl")
export HMC!, Hamiltonian, OMF4!, leapfrog!, update_momenta!, update_field!

end # module
