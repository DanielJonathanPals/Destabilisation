# Here is a collection of functions that can be used to test wether certain objects have the correct format

module FormatTests

export check_traj

"""
    check_traj(traj)

Checks if `traj` has the format of a trajectory i.e. if `traj` is a two dimentional matrix with all entries
    given by numbers.
"""
function check_traj(traj)
    if !(typeof(traj) <: Matrix)
        error("You are currently treating an object as an trajectory which is not a matrix")
    end
    if any(x -> !(typeof(x) <: Number), traj)
        error("You are currently treating an object as an trajectory which does not only contain numbers")
    end
end

end