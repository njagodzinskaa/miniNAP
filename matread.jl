using MAT, HDF5, DataFrames

filename = "/Users/jjc/Downloads/OWT220207_2A_DIV63_BASELINE.mat"
vars = matread(filename)

els = split(vars["header"], "\n")
re = r"/^\d+$/"

m = match.(r"\d{2}", split(els[7], ";"))

m = [m[x].match for x in 1:length(m)]


info = DataFrame(Label = "")

info.Label