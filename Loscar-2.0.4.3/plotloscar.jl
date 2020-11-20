using DelimitedFiles, Plots, Interpolations
co2doubling = 3;

time_vals = readdlm("tmv.dat", '\t', Float32, '\n')  
d13c = readdlm("d13c.dat")

