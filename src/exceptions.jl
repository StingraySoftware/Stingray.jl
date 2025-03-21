module Exceptions

export StingrayError

struct StingrayError <: Exception
    msg::String
end

function Base.showerror(io::IO, e::StingrayError)
    print(io, "StingrayError: ", e.msg)
end

end 
