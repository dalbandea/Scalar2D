using Revise
using Random
import Pkg
import RecipesBase.plot
using ADerrors
using Statistics, Plots
using DelimitedFiles
Pkg.develop(path="/home/david/.julia/personal-utils/Utils/")
using Utils
Pkg.activate(".")
using Scalar2D, ReverseDiff
pyplot()



filterdir(dir::String, text::String) = filterdir(dir, [text])

"""
Returns vector of Strings with occurrences in the directory using regex, but must escape \$ with \\.
"""
function filterdir(dir::String, texts::Array{String})
    dirfiles = readdir(dir)
    occurrences = filter(s -> occursin(Regex("$(texts[1])"), s), dirfiles)
    for text in texts
        occurrences = filter(s -> occursin(Regex("$text"), s), occurrences)
    end
    return joinpath.(dir, occurrences)
end


function extract_symmetrize(corr_file::String, ID::String; burnout::Int64=1) #{{{
	# Load correlation data
    corr = readdlm(corr_file)[burnout:end,:]
	T = length(corr[1,:]);
	nconf = length(corr[:,1]);

	println("Size: ", size(corr))
	println("Temporal length: ", T)
	
	T2p1 = convert(Int64, T/2+1)
	xdata = range(0, length=T2p1)
	ydata = Vector{uwreal}(undef, T2p1);

	for i in eachindex(ydata)
		ydata[i] = uwreal( (corr[:,i] .+ corr[:,1+(T-i+1)%T])/2, ID)
		uwerr(ydata[i])
	end

	return xdata, ydata, T
end #}}}

function extract_symmetrize(corr_file::String, ID::String; burnout::Int64=1) #{{{
	# Load correlation data
    corr = readdlm(corr_file)[burnout:end,:]
	T = length(corr[1,:]);
	nconf = length(corr[:,1]);

	println("Size: ", size(corr))
	println("Temporal length: ", T)
	
	T2p1 = convert(Int64, T/2+1)
	xdata = range(0, length=T2p1)
	ydata = Vector{uwreal}(undef, T2p1);

	for i in eachindex(ydata)
		ydata[i] = uwreal( (corr[:,i] .+ corr[:,1+(T-i+1)%T])/2, ID)
		uwerr(ydata[i])
	end

	return xdata, ydata, T
end #}}}

function plot_correlator(corr_file::String, ID::String; burnout::Int64=1) #{{{
	# Load correlation data
    corr = readdlm(corr_file)[burnout:end,:]
	T = length(corr[1,:]);

	println("Size: ", size(corr))
	println("Temporal length: ", T)
	
	xdata = range(0, length=T)
	ydata = Vector{uwreal}(undef, T);

	for i in eachindex(ydata)
		ydata[i] = uwreal(corr[:,i], ID)
		uwerr(ydata[i])
	end

	pl = plot(xdata, value.(ydata), yerr=err.(ydata), seriestype=:scatter, yscale=:identity)
	display(pl)

    return xdata, ydata, T

end #}}}

function tmin_fit(f::UtilsFunc , xdata, ydata::Vector{uwreal}, tmin::Int64, prms0::Vector{Float64}; Tmax::Int64=length(ydata)-1, plot_tmin_fit::Bool=false, yscale::Symbol=:identity) #{{{

	if (Tmax > xdata[end][1])
		throw("Tmax value greater than greatest value in xdata")
	end

	pos_tmin = findfirst(x->x[1]==tmin, xdata)
	pos_Tmax = findfirst(x->x[1]==Tmax, xdata)
	fit_region = pos_tmin:pos_Tmax

	# Fit
	fitp, cse, cs = fit_data(f, xdata[fit_region], ydata[fit_region], prms0)
	uwerr.(fitp) # compute errors

	# Plot fit
	pyplot()
	if(plot_tmin_fit==true)
		if(size(xdata[1],1)==1)
			pl = plot_fit(xdata, ydata, f, fitp)
			plot!(pl, xdata, value.(ydata), yerr=err.(ydata), seriestype=:scatter, title="tmin = "*string(tmin), yscale=yscale)
			# plot!(ylim=(9000,10000))
			display(pl)
		else
			res = f.(xdata, [fitp for i in 1:length(xdata)])
			uwerr.(res)
			pl_x = (hcat(xdata...) |> permutedims)[:,1]
			pl = plot( pl_x, value.(res), ribbons=err.(res), reuse=false, title="tmin="*string(tmin)*", 1")
			plot!(pl, pl_x, value.(ydata), yerr=err.(ydata), seriestype=:scatter)
			display(pl)
		end
	end
	
	return fitp, cse, cs

end #}}}

function tmin_loop(f::UtilsFunc , xdata, ydata, tmin::Int64, tmax::Int64, prms0::Vector{Float64}; update_prms::Bool=true, Tmax::Int64=length(ydata)-1, plot_tmin_fit::Bool=false, plot_column::Int64=0) #{{{

	fitps = [] # tmin and parameters to be returned
	tminvalues = tmin:tmax

	prms = prms0

	for itmin in tminvalues

		# Fit
		fitp, cse, cs = tmin_fit(f, xdata, ydata, itmin, prms, Tmax=Tmax, plot_tmin_fit=plot_tmin_fit, yscale=:log10)
		display(cs)

		if(itmin == tmin)
			fitps = hcat(itmin, permutedims( fitp ) )
		else
			aux_fitps = hcat(itmin, permutedims( fitp ))
			fitps = vcat(fitps, aux_fitps)
		end

		if(update_prms)
			prms = value.(fitp)
		end
	end

	if(plot_column > 0)
		pl = plot(fitps[:,1], value.(fitps[:,plot_column]), yerr=err.(fitps[:,plot_column]), seriestype=:scatter)
		display(pl)
	end

	return fitps
end #}}}

function derivate_sym_correlator(corr_x, corr_y) #{{{

	T2p1 = length(corr_y)
	dcorr_x = corr_x[2:end-1]
	dcorr_y = Vector{uwreal}(undef, T2p1-2 )

	for i in eachindex(dcorr_y)
		dcorr_y[i] = (corr_y[i+2]-corr_y[i])/2
		uwerr(dcorr_y[i])
	end

	return dcorr_x, dcorr_y
end #}}}

function eq_44(ydata::Array{uwreal})
    xis = Vector{uwreal}()
    for i in 1:length(ydata)-2
        xi = 1 / acosh((ydata[i] + ydata[i+2])/(2*ydata[i+1]))
        uwerr(xi)
        push!(xis, xi)
    end

    return xis
end

plot(xdata, ydata::Vector{uwreal}; kwargs...) = plot(xdata, value.(ydata), yerr=err.(ydata), seriestype=:scatter, kwargs...)
    

wdir = "./results/1-Checks/L20_b0.645_l0.5_D2022-05-13-16-06-25/measurements/"

burnout = 1000

M = readdlm(filterdir(wdir, "mag")[1])[burnout:end]

strID = wdir

uwM = uwreal(vcat(M...), strID)

Mabs = uwreal(abs.(M), strID)
uwerr(Mabs)
Mabs

corrs = readdlm(filterdir(wdir, "corr")[1])[burnout:end]

xdata, ydata, T = extract_symmetrize(filterdir(wdir, "corr")[1], strID, burnout=burnout)

corrf = similar(ydata)
for i in eachindex(corrf)
    corrf[i] = ydata[i] - uwM^2
    uwerr(corrf[i])
end

xis = eq_44(corrf)

plot(2:19, xis)


xdata, ydata, T = plot_correlator(filterdir(wdir, "corr")[1], strID, burnout=burnout)

f = SymCorrelator(T, 1, 2)

fitp, cse, cs = tmin_fit(f, xdata, corrf, 8, [0.1, 0.2, 1.0, 1.0])

pl = plot_fit(xdata, corrf, f, fitp)
plot!(pl, xdata, value.(corrf), yerr=err.(corrf), seriestype=:scatter)

fitps = tmin_loop(f, xdata, corrf, 1,8, [0.05, 0.20, 1.5, 1.0], plot_column=3)

xis = []
for i in 1:length(corrf)-2
    xi = 1 / acosh((corrf[i] + corrf[i+2])/(2*corrf[i+1]))
    uwerr(xi)
    push!(xis, xi)
end

xi

dxdata, dcorrf = derivate_sym_correlator(xdata, corrf)

g = SymCorrelator(T=T, s=-2)

fitp, cse, cs = tmin_fit(g, dxdata, dcorrf, 7, [0.1, 0.2, 1.0, 1.0])

fitps = tmin_loop(g, dxdata, dcorrf, 1,7, [0.05, 0.20, 1.0, 1.0], plot_column=3)
