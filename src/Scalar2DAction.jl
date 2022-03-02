
function action(phi, prm) where {T}

    # For type stability initialize to correct type
    # https://www.juliabloggers.com/writing-type-stable-julia-code/
    act = zero(eltype(phi))

    Nx = size(phi,1)
    Ny = size(phi,2)
    for i in 1:Nx, j in 1:Ny
        iu = mod1(i+1, Nx)
        ju = mod1(j+1, Ny)
        act = act + (phi[i,j]-phi[iu,j])^2 + (phi[i,j]-phi[i,ju])^2 +
            prm.msq*phi[i,j]^2 + prm.lam*phi[i,j]^4
    end

    return act/2
end

# Force is computed in place for performance
function force!(frc, phi, prm) 

    Nx = size(phi,1)
    Ny = size(phi,2)
    for i in 1:Nx, j in 1:Ny
        iu = mod1(i+1, Nx)
        ju = mod1(j+1, Ny)
        id = mod1(i-1, Nx)
        jd = mod1(j-1, Ny)
        
        frc[i,j] = phi[iu,j] + phi[id,j] + phi[i,ju] + phi[i,jd] -
            (prm.msq+4)*phi[i,j] - 2*prm.lam*phi[i,j]^3
    end

    return nothing
end

# cf_tape has recorded the operations to do the derivative
# https://github.com/JuliaDiff/ReverseDiff.jl/blob/master/examples/gradient.jl
forcead!(frc, phi, cf_tape) where {T} =  ReverseDiff.gradient!(frc, cf_tape, phi)


