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

circindex(i::Int,N::Int) = 1 + mod(i-1,N) # for symmetric BC

function correlation_function(phi, y)
    L = size(phi)[1]
    # C = phi .* (circshift(phi, (y,0)) .+ circshift(phi, (0,y)))
    C = phi .* circshift(phi, (0,y))

    return sum(C)/L^2
end

# function correlation_function(phi, y)
#     C = 0.0
#     L = size(phi)[1]

#     for i in 1:L
#         for j in 1:L
#             C += phi[i,j] * (phi[circindex(i+y,L), j] + phi[i, circindex(j+y, L)])
#         end
#     end

#     return C/L^2/2
# end

function correlation_function(phi)
    L = size(phi)[1]
    C = Vector{Float64}(undef, L)

    for i in 1:L
        C[i] = correlation_function(phi, i-1)
    end

    return C
end
