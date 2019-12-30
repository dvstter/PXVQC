module MyVQC

using Zygote
using Zygote: @adjoint

include("tensorop.jl")
include("struct.jl")
include("gates.jl")
include("circuit.jl")
include("diff.jl")

end # end for module
