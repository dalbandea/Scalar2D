using Revise
using Random
import Pkg
import Dates
Pkg.activate(".")
using Scalar2D, ReverseDiff
using DelimitedFiles
using ArgParse
using InteractiveUtils

# julia16 main/scalar-main.jl -L 10 -b 0.5 --nsteps 10 -n 1 --wdir=trash/
# julia16 --threads=auto  main/scalar-main.jl -L 80 -b 0.677 -l 0.5 --nmeas 400 --nsteps 18 -n 2000000 --wdir=results/1-Checks/
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

        "--nmeas"
        help = "measure every nmeas trajectories"
        required = false
        arg_type = Int
        default = 1

        "-n"
        help = "number of trajectories"
        required = true
        arg_type = Int

        "--wdir"
        help = "path to directory to save configurations and logs"
        required = true
        arg_type = String
        # default = "configs/"

        "--replica"
        help = "Replica number. If !=0, wdir must point to existing directory
        with existing replica 0"
        required = false
        default = 0
        arg_type = Int
    end

    return parse_args(s)
end


########################
# PARSE ARGUMENTS ######
########################

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
n_meas  = parsed_args["nmeas"]
replica = parsed_args["replica"]


###########################
# CREATE DIRECTORIES ######
###########################

wdir_prefix = "L"*string(prm.iL[1])*"_b"*string(beta)*"_l"*string(lambda)
dt = Dates.now()
wdir_sufix = "_D"*Dates.format(dt, "yyyy-mm-dd-HH-MM-SS")

if replica == 0
    wdir = joinpath(parsed_args["wdir"], wdir_prefix*wdir_sufix, "0-r0")
elseif replica > 0
    if isdir(joinpath(parsed_args["wdir"], "0-r0"))
        wdir = joinpath(parsed_args["wdir"], string(replica)*"-r"*string(replica))
    else
        error("Replica 0 folder not found")
    end
else
    error("Replica number not valid")
end

# trururu

logdir = joinpath(wdir, "logs")
logfile = joinpath(logdir, "log.txt")
mesdir = joinpath(wdir, "measurements")
magfile = joinpath(mesdir, wdir_prefix*"_mag.txt")
suscfile = joinpath(mesdir, wdir_prefix*"_susc.txt")
susctfile = joinpath(mesdir, wdir_prefix*"_susct.txt")
g0file = joinpath(mesdir, wdir_prefix*"_g0.txt")
corrfile = joinpath(mesdir, wdir_prefix*"_corr.txt")

if isdir(wdir) == false
    mkpath(wdir)
    mkpath(logdir)
    mkpath(mesdir)
else
    error("Directory wdir already exists.")
end


#############
# LOGS ######
#############

global io_stat = open(logfile, "a")
versioninfo(io_stat)
write(io_stat, "\n\n")
write(io_stat, "$(PROGRAM_FILE) $(join(ARGS, " "))")
write(io_stat, "\n\n")
write(io_stat, "RNG_state: ", string(Random.default_rng()))
write(io_stat, "\n")
write(io_stat, "Recover with `copy!(Random.default_rng(), state)`")
close(io_stat)

# Create patch with non-staged changed (from tracked files)
write(joinpath(logdir, "gitpatch.patch"), readchomp(`git diff`))

global io_stat = open(joinpath(logdir, "branchinfo.txt"), "w")
# Print current branch to branchinfo file
write(io_stat, "Current branch: ")
write(io_stat, readchomp(`git rev-parse --abbrev-ref HEAD`))
write(io_stat, "\n")
# Print current commit of branch to branchinfo file
write(io_stat, "Current commit: ")
write(io_stat, readchomp(`git rev-parse --short HEAD`))
write(io_stat, "\n\n")
# Print current estate of the repository to branchinfo file
write(io_stat, "Estate of repository:")
write(io_stat, "\n")
write(io_stat, readchomp(`git show-ref`))
close(io_stat)



#############
# HMC #######
#############

# Initialize phi field with random numbers
phi = zeros(prm.iL[1], prm.iL[2])
randn!(phi)

# Perform n_traj HMC steps
acc = Vector{Int64}()

for i in 1:n_traj
    # println("Config $i")

    @time HMC!(phi, epsilon, nsteps, acc, prm, integrator = Leapfrog());

    if i % n_meas == 0
        global io_stat = open(magfile, "a")
        write(io_stat, "$(magnetization(phi))\n")
        close(io_stat)

        # global io_stat = open(g0file, "a")
        # write(io_stat, "$(G0(phi))\n")
        # close(io_stat)

        global io_stat = open(suscfile, "a")
        write(io_stat, "$(susceptibility(phi))\n")
        close(io_stat)

        phi_t = he_flow(phi, prm, prm.iL[1]^2 / 64)
        global io_stat = open(susctfile, "a")
        write(io_stat, "$(susceptibility(phi_t))\n")
        close(io_stat)

        # global io_stat = open(corrfile, "a")
        # for value in correlation_function(phi)
        #     write(io_stat, "$(value)\t")
        # end
        # write(io_stat, "\n")
        # close(io_stat)
    end
end
