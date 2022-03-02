function HMC!(phi, eps, ns, acc, prm::LattParm; AD_tape = 0)

    phi_cp = similar(phi)
    phi_cp .= phi
    
    mom = Random.randn(Float64, prm.iL[1], prm.iL[2])
    hini = Hamiltonian(mom, phi, prm)

    # OMF4!(mom, phi, eps, ns, prm::LattParm)
    leapfrog!(mom, phi, eps, ns, prm::LattParm, AD_tape)
    
    hfin = Hamiltonian(mom, phi, prm)

    pacc = exp(-(hfin-hini))
    if (pacc < 1.0)
        r = rand()
        if (pacc > r) 
            push!(acc, 1)
        else
            phi .= phi_cp
            push!(acc, 0)
        end
    else
        push!(acc, 1)
    end

    if (acc[end] == 0)
        @info("    REJECT: Energy [inital: $hini; final: $hfin; difference: $(hfin-hini)]; Pacc = $(count(acc.==1)/length(acc))")
    else
        @info("    ACCEPT:  Energy [inital: $hini; final: $hfin; difference: $(hfin-hini)]; Pacc = $(count(acc.==1)/length(acc))")
    end
    
    return nothing
end


function Hamiltonian(mom, phi, prm::LattParm)

    hini = mapreduce(x -> x^2, +, mom)/2.0 + action(phi, prm)

    return hini
end


function OMF4!(mom, phi, eps, ns, prm::LattParm)

    r1::Float64 =  0.08398315262876693
    r2::Float64 =  0.2539785108410595
    r3::Float64 =  0.6822365335719091
    r4::Float64 = -0.03230286765269967
    r5::Float64 =  0.5-r1-r3
    r6::Float64 =  1.0-2.0*(r2+r4)

    frc = zeros(Float64, prm.iL[1], prm.iL[2])

    for i in 1:ns
        # STEP 1
        force!(frc, phi, prm)
        mom .= mom .+ (r1*eps).*(frc)
        update_field!(phi, mom, r2*eps) 

        # STEP 2
        force!(frc, phi, prm)
        mom .= mom .+ (r3*eps).*(frc)
        update_field!(phi, mom, r4*eps) 

        # STEP 3
        force!(frc, phi, prm)
        mom .= mom .+ (r5*eps).*(frc)
        update_field!(phi, mom, r6*eps) 

        # STEP 4
        force!(frc, phi, prm)
        mom .= mom .+ (r5*eps).*(frc)
        update_field!(phi, mom, r4*eps) 

        # STEP 5
        force!(frc, phi, prm)
        mom .= mom .+ (r3*eps).*(frc)
        update_field!(phi, mom, r2*eps) 

        # STEP 6
        force!(frc, phi, prm)
        mom .= mom .+ (r1*eps).*(frc)
    end

    return nothing
end


function leapfrog!(mom, phi, eps, ns, prm::LattParm, AD_tape)

	# First half-step for momenta
	update_momenta!(mom, phi, eps/2.0, prm, AD_tape)

	# ns-1 steps
	for i in 1:(ns-1) 
		# Update gauge links
        update_field!(phi, mom, eps) 

		#Update momenta
        update_momenta!(mom, phi, eps, prm, AD_tape)
	end
	# Last update for gauge links
    update_field!(phi, mom, eps) 

	# Last half-step for momenta
    update_momenta!(mom, phi, eps/2.0,prm, AD_tape)

	return nothing
end


function update_momenta!(mom, phi, eps, prm, AD_tape)
    frc = zeros(Float64, prm.iL[1], prm.iL[2])

    if AD_tape == 0
        force!(frc, phi, prm)
    elseif isa(AD_tape, ReverseDiff.CompiledTape)
        forcead!(frc, phi, AD_tape)
    else
        throw("Incorrect tape")
    end
    mom .= mom .+ eps .* frc

    return nothing
end

function update_field!(phi, mom, eps)

    phi .= phi + eps * mom
    
    return nothing
end
