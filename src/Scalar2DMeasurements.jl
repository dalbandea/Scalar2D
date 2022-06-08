function magnetization(phi)
    return sum(phi) / length(phi)
end

function G0(phi)
    return sum(phi.^2) / length(phi)
end

# function susceptibility(phi)
#     return sum(phi[1] * phi) / length(phi)
# end

function susceptibility(phi)
    return sum(phi.^2) / length(phi)
end

function chi2(mags, prm)
    V = prm.iL[1] * prm.iL[2]
    M = sum(mags) / length(mags)
    M2 = mags.^2
    chi = V * (M2 .- M^2)
    return chi
end

function correlation_function(phi, shift::Tuple{Int64, Int64})
    L = size(phi)[1]
    C = phi .* circshift(phi, shift)

    return sum(C)/L^2
end

circindex(i::Int,N::Int) = 1 + mod(i-1,N) # for symmetric BC

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
    C = zeros(L)

    for i in 1:L
        for j in 1:L
            C[i] += correlation_function(phi, (i-1, j-1))
        end
    end

    return C ./ L
end

function correlation_matrix(phi)
    L = size(phi)[1]
    C = zeros((L,L))

    for i in 1:L
        for j in 1:L
            C[i,j] = correlation_function(phi, (i-1, j-1))
        end
    end

    return C
end
