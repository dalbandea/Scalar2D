function magnetization(phi)
    return sum(phi) / length(phi)
end

function susceptibility(phi)
    return sum(phi[1] * phi) / length(phi)
end

function chi2(mags, prm)
    V = prm.iL[1] * prm.iL[2]
    M = sum(mags) / length(mags)
    M2 = mags.^2
    chi = V * (M2 .- M^2)
    return chi
end
