function he_flow_p(phi, prm, t)
    L = prm.iL[1]
    ep2t = zeros((L,L)) 

    for j in 1:L, i in 1:L
        p2 = 4 * (sin(2 * pi * (i-1) / L/ 2)^2 + sin(2 * pi * (j-1) / L / 2)^2)
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
