using Revise, DelimitedFiles, Random
import Pkg
Pkg.activate(".")
using Scalar2D
using ArgParse


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--infile"
        help = "path to file containing directories with replicas. Will also be working directory"
        required = true
        arg_type = String
    end

    return parse_args(s)
end


########################
# PARSE ARGUMENTS ######
########################

parsed_args = parse_commandline()

infile  = parsed_args["infile"]

########################
# LOGS #################
########################

wdir    = dirname(infile)
paths   = vec(readdlm(infile, String))
logdir  = joinpath(wdir, "logs")
out_file = joinpath(wdir, "ana.out")

if isdir(wdir)
    mkpath(logdir)
else
    error("Parent directory of input file does not exist")
end

create_patch(logdir)
create_branchinfo(logdir)


########################
# ANALYSIS #############
########################

obs_lst = ["mag"]
transform_lst  = [[identity, x -> abs.(x)]]

ana_observables(paths, obs_lst, transform_lst, skp=100, out_file = out_file)
