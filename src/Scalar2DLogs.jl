"""
Create patch with non-staged changed (from tracked files)
"""
function create_patch(logdir::String)
    write(joinpath(logdir, "gitpatch.patch"), readchomp(`git diff`))
end

function create_branchinfo(logdir::String)
    io_stat = open(joinpath(logdir, "branchinfo.txt"), "w")
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
end
