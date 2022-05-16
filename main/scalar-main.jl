using Revise
using Random
import Pkg
import Dates
Pkg.activate(".")
using Scalar2D, ReverseDiff
using DelimitedFiles
using ArgParse

# julia16 main/scalar-main.jl -L 10 -b 0.5 --nsteps 10 -n 1 --wdir=trash/
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "-L"
        help = "lattice size"
        required = true
        arg_type = Int

        "-l"
        help = "lambda"
        required = false
        arg_type = Float64
        default = 0.5

        "-b"
        help = "beta"
        required = true
        arg_type = Float64

        "-t"
        help = "tau"
        required = false
        arg_type = Float64
        default = 1.0

        "--nsteps"
        help = "integration steps per trajectory"
        required = true
        arg_type = Int

        "-n"
        help = "number of trajectories"
        required = true
        arg_type = Int

        "--wdir"
        help = "path to directory to save configurations and logs"
        required = true
        arg_type = String
        # default = "configs/"
    end

    return parse_args(s)
end

parsed_args = parse_commandline()

# Theory Parameters
lsize_1 = parsed_args["L"]
lsize_2 = lsize_1
lambda  = parsed_args["l"]
beta    = parsed_args["b"]
prm     = LattParmB((lsize_1, lsize_2), beta, lambda)

# HMC parameters
tau     = parsed_args["t"]
nsteps  = parsed_args["nsteps"]
epsilon = tau/nsteps
n_traj  = parsed_args["n"]

# Initialize phi field with random numbers
phi = zeros(prm.iL[1], prm.iL[2])
randn!(phi)

wdir_prefix = "L"*string(prm.iL[1])*"_b"*string(beta)*"_l"*string(lambda)
dt = Dates.now()
wdir_sufix = "_D"*Dates.format(dt, "yyyy-mm-dd-HH-MM-SS")
wdir = joinpath(parsed_args["wdir"], wdir_prefix*wdir_sufix)

logdir = joinpath(wdir, "logs")
logfile = joinpath(logdir, "log.txt")
mesdir = joinpath(wdir, "measurements")
magfile = joinpath(mesdir, wdir_prefix*"_mag.txt")
corrfile = joinpath(mesdir, wdir_prefix*"_corr.txt")

if isdir(wdir) == false
    mkpath(wdir)
    mkpath(logdir)
    mkpath(mesdir)
else
    error("Directory wdir already exists.")
end

# Perform n_traj HMC steps
acc = Vector{Int64}()

for i in 1:n_traj
    println("Config $i")

    @time HMC!(phi, epsilon, nsteps, acc, prm, integrator = Leapfrog());

    global io_stat = open(magfile, "a")
    write(io_stat, "$(magnetization(phi))\n")
    close(io_stat)

    global io_stat = open(corrfile, "a")
    for value in correlation_function(phi)
        write(io_stat, "$(value)\t")
    end
    write(io_stat, "\n")
    close(io_stat)
end
