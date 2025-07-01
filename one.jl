run(`clear`) # clear REPL

using Pkg
Pkg.activate(".")
# Load libraries
using Glob, HDF5, DSP, StatsBase, Findpeaks, LinearAlgebra, MatrixEquations, ColorSchemes, CSV, DataFrames, MAT
using Combinatorics: combinations
using LinearAlgebra
using SparseArrays
using Arpack
using Plots
using CSV

include("Utils.jl") # Link script containing the functions
theme(:dracula) # Set plotting theme

# Specify channel layout
# channels = [47.0; 48.0; 46.0; 45.0; 38.0; 37.0; 28.0; 36.0; 27.0; 17.0; 26.0; 16.0; 35.0; 25.0; 15.0; 14.0; 24.0; 34.0; 13.0; 23.0; 12.0; 22.0; 33.0; 21.0; 32.0; 31.0; 44.0; 43.0; 41.0; 42.0; 52.0; 51.0; 53.0; 54.0; 61.0; 62.0; 71.0; 63.0; 72.0; 82.0; 73.0; 83.0; 64.0; 74.0; 84.0; 85.0; 75.0; 65.0; 86.0; 76.0; 87.0; 77.0; 66.0; 78.0; 67.0; 68.0; 55.0; 56.0; 58.0; 57.0]
exclude_chans = [11, 18, 81, 88]
channels = [i*10 + j for i in 1:8 for j in 1:8 if !(i*10 + j in exclude_chans)]

# Specify path to where the .h5 files are stored
# ppath = "/Volumes/NO NAME/StimExp"
# ppath = "/Users/jjc/MEANAP/Data/StimExp"
# ppath = "/Users/jjc/MEANAP/Data/"
ppath = "/Users/nataliajagodzinska/Desktop/SAND_summer"

# processed_files = glob("*.h5", ppath)

processed_files = ["/Users/nataliajagodzinska/Desktop/SAND_summer/OWT220207_2C_DIV63_BASELINE.mat"]

# Remove the ending of each filename
# processed_trunc = map(x -> x[length(ppath)+1:end-3], processed_files)
processed_trunc = map(x -> x[length(ppath)+1:end-4], processed_files)


# Create empty dataframes
df_all = DataFrame()
res = DataFrame(recording=[], high_coord=[], high_AC=[], high_fr=[], low_coord=[], low_AC=[], low_fr=[])

@time main(processed_files[1])



