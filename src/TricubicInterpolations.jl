module TricubicInterpolations

export Tricubic
export partial_derivative_x, partial_derivative_y, partial_derivative_z

include("cube.jl")
include("finite_difference.jl")
include("interpolation.jl")

end # module TricubicInterpolations
