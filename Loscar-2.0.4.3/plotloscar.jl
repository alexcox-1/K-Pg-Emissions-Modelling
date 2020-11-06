using DelimitedFiles, Plots, Interpolations
co2doubling = 3;

time_vals = readdlm("tmv.dat", '\t', Float32, '\n')  
pco2 = readdlm("pco2a.dat", '\t', Float32, '\n')
temp = (pco2./600) .- 1;
temp = 3 .*temp;
plot(time_vals,pco2,lw=2)
