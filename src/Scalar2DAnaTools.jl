using ADerrors, DelimitedFiles, Statistics

function read_replicas(path::String, obs_string::String; burnout::Int64 = 1000,
        skp::Int64 = 1)
    isdir(path) ? nothing : error("path is not a directory")
    find_str = """.*$obs_string\\.txt\$"""
    cmd1 = `find $path -regex $find_str`
    obs_files = pipeline(cmd1, `sort -V`) |> readchomp

    obs_files = split(obs_files, '\n')

    obs = Vector{Float64}()
    lengths = Vector{Int64}()
    for file in obs_files
        data = vec(readdlm(file))[burnout:skp:end]
        append!(obs, data...)
        append!(lengths, length(data))
    end
    return obs, lengths
end

ana_observables(rep_path::String, obs_list::Vector{String}; kwargs...) = ana_observables(rep_path, obs_list, [identity for i in obs_list]; kwargs...)

# Fallback for input vector of transforms
function ana_observables(rep_path::String, obs_list::Vector{String}, transforms::Vector; kwargs...)
    vectorized_transforms = [ i isa Vector ? i : [i] for i in transforms ] |> x -> convert(Vector{Vector{Function}}, x)
    return ana_observables(rep_path, obs_list, vectorized_transforms; kwargs...)
end


function ana_observables(rep_paths::Vector{String}, args...; kwargs...)
    for path in rep_paths
        ana_observables(path, args...; kwargs...)
    end
    return nothing
end

"""
Given a path directory, analyzes observables in `obs_list` transformed with the functions in `transform`, and prints mean and autocorrelation for each transform.

- `transform`: must be a vector of functions. If for a given observables one wants to perform several transforms, then put them as a vector inside the vector. If not provided, it defaults to the identity transformation `identity`

Examples
========
    path = "my/path/"
    obs_lst = ["mag", "susc"]
    transform_lst = [identity, [x -> x.^2, x -> abs.(x)]]]
    ana_observables(path, obs_lst, transform_lst, burnout = 1000, skp=1)
"""
function ana_observables(rep_path::String, obs_list::Vector{String}, transforms::Vector{Vector{T}}; burnout::Int64 = 1000, skp::Int64 = 1, out_file::String) where T <: Function
    if length(transforms) != length(obs_list)
        error("Legnth of transforms != Length of obs_list")
    end

    for i in eachindex(obs_list)
        obs, lengths = read_replicas(rep_path, obs_list[i], burnout = burnout, skp = skp)
        for transform in transforms[i]
            ana_observable(transform(obs), rep_path*obs_list[i], lengths, skp = skp, out_file = out_file)

        end

        if i == length(obs_list)
            let
                acc     = parse_acc(rep_path)
                nsteps  = parse_log(rep_path)
                
                open(out_file, "a") do io
                    write(io, "$(sum(lengths)*skp),$(nsteps),$(acc)\n")
                end;
            end
        end
    end
    return nothing
end

function ana_observable(obs, ID::String, lengths::Vector{Int64}; skp::Int64, out_file::String)
        uwobs = uwreal(obs, ID, lengths)
        uwerr(uwobs)
        tau = taui(uwobs, ID)
        dtau = dtaui(uwobs, ID)
        open(out_file, "a") do io
            write(io, "$(value(uwobs)) +/- $(err(uwobs)),$(tau*skp) +/- $(dtau*skp),")
        end
        return nothing
end

function parse_acc(path::String)
    find_str = """.*mag\\.txt\$"""
    cmd1 = `find $path -regex $find_str`
    cmd2 = `head -n 1`
    mag_file = pipeline(cmd1, cmd2) |> readchomp
    mags =  vec(readdlm(mag_file))
    accs = mags .!= circshift(mags, -1)
    return round(mean(accs), digits=2)
end

function parse_log(path::String)
    find_str = """.*log\\.txt\$"""
    cmd1 = `find $path -regex $find_str`
    cmd2 = `head -n 1`
    cmd3 = `xargs cat`
    cmd4 = `grep """^{.*}\$"""`
    cmd5 = `sed """s/[{}]//g"""`
    cmd6 = `sed """s/, /\\n/g"""`
    cmd7 = `grep """nsteps"""`
    cmd8 = `awk """{print \$NF}"""`
    nstps = pipeline(cmd1, cmd2, cmd3, cmd4, cmd5, cmd6, cmd7, cmd8) |> readchomp
    return parse(Int64, nstps)
end
