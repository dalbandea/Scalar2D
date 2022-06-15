using FFTW
using Plots

t0 = 0              # Start time 
fs = 44100          # Sampling rate (Hz)
tmax = 0.01          # End time       
t = t0:1/fs:tmax;   
signal = sin.(2π * 60 .* t)

F = fftshift(fft(signal))
freqs = fftshift(fftfreq(length(t), fs))

# plots 
time_domain = plot(t, signal, title = "Signal", label='f',legend=:top, seriestype = :scatter)
freq_domain = plot(freqs, abs.(F), title = "Spectrum", xlim=(-100, +100), xticks=-100:20:100, label="abs.(F)",legend=:top, seriestype = :scatter) 
plot(time_domain, freq_domain, layout = (2,1))

###

sr = 100
tmax = 2.0
ts = 1.0/sr
t = 0:ts:tmax
f = 1.5
x = sin.(2*pi*f*t)
time_domain = plot(t, x, seriestype = :scatter);
X = fft(x)
freqs = fftfreq(length(t),sr)
freq_domain = plot(freqs, abs.(X), seriestype = :scatter, xlims=(-2,2));
plot(time_domain, freq_domain, layout = (2,1));

plot(freqs)

plot(abs.(myft(x)), xlims=(0,20), seriestype = :scatter)

function myft(x)
    N = length(x)
    X = zeros(ComplexF64, N)

    for k in 1:N, n in 1:N
        X[k] += exp(-im * 2*pi/N * (k-1) * (n-1)) * x[n]
    end
    return X
end


### 2 dimensions phi4

function fill_ept(phi, t::Float64)
    L = size(phi)[1]
    ept = zeros((L,L)) 
    for j in 1:L, i in 1:L
        ept[i,j] = exp(-4*(sin(2*pi*(i-1)/L/2)^2 + sin(2*pi*(j-1)/L/2)^2)*t)
    end

    return ept
end

function he_flow_p(phi, t)
    L = size(phi)[1]
    ept = zeros((L,L)) 
    for j in 1:L, i in 1:L
        ept[i,j] = exp(-4*(sin(2*pi*(i-1)/L/2)^2 + sin(2*pi*(j-1)/L/2)^2)*t)
    end

    res = ept .* phi

    return res
end

X = randn((10,10))

function he_flow(phi, t::Float64)
    phi_p0 = fft(phi)
    phi_pt = he_flow_p(phi_p0, t)
    phi_t = ifft(phi_pt)
end

function he_flow(phi, ept)
    phi_p0 = fft(phi)
    phi_pt = phi_p0 .* ept
    phi_t = ifft(phi_pt)
    return real.(phi_t)
end

X_p0 = fft(X)
X_pt = he_flow_p(X_p0, 1)
X_t = ifft(X_pt)


### finite step flow

function force(phi)
    L = size(phi)[1]
    frc = zeros((L,L)) .- 4*phi
    for shift in [(0,1), (0,-1), (1,0), (-1,0)]
        frc .+= circshift(phi, (shift))
    end

    return frc
end

function finite_he_flow(phi, t, ns)
    phi_t = copy(phi)
    eps = t/ns
    for i in 1:ns
        phi_t .= phi_t .+ eps * force(phi_t) 
    end

    return phi_t
end



### Checks

X = randn((10,10))

finite_he_flow(X, 1, 10000) .- he_flow(X, 1) .|> abs |> maximum

he_flow(he_flow(X,1), -1) .- X .|> abs |> maximum

finite_he_flow(finite_he_flow(X,1, 100000), -1, 100000) .- X .|> abs |> maximum
