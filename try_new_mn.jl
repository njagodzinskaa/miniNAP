

include("Utils.jl") # Link script containing the functions

using MAT

dt = loadDataMat("/Users/jjc/Downloads/OWT220207_2A_DIV63_BASELINE.mat")

[print(x) for x in eachindex(dt[2])]

