hat_p2(i::Int64, j::Int64, L::Int64) = 4 * (sin(pi * i / L)^2 + sin(pi * j / L)^2)

"""
Returns 
"""
function he_flow_p(phi, prm, t)
    L = prm.iL[1]
    ep2t = zeros((L,L)) 

    for j in 1:L, i in 1:L
        p2 = hat_p2(i-1, j-1, L)
        ep2t[i,j] = exp(-p2*t)
    end

    res = ep2t .* phi

    return res
end

function he_flow(phi, prm, t)
    phi_p0 = FFTW.fft(phi)
    phi_pt = he_flow_p(phi_p0, prm, t)
    phi_t  = FFTW.ifft(phi_pt)

    return real.(phi_t)
end
