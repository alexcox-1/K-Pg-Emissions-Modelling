# model time goes from 0 to 1.4 million years
# 301 bins
# create a 3 column csv
# read in the bsr temperature curve
# define doublingrate
timestep = 5000;
nbins = 301;

timevals = 0:timestep:timestep*(nbins-1);
carbonvals = zeros(nbins);
sulfatevals = zeros(nbins);

# writedlm a loscar input file (tab separated) for carbon and sulfate
# write loscar input file (call a function from another file) use interpolation
# compile loscar using run(`make loscar PALEO=1`)

#= time_vals = readdlm("tmv.dat", '\t', Float32, '\n')  
pco2 = readdlm("pco2a.dat", '\t', Float32, '\n')
temp = (pco2./600) .- 1;
temp = doublingrate .*temp; 
temp = temp - sulfate temp ( t = -11.3 * (1-e(-0.0466*pinatubos/year))) =#

#= log likelihoods time, mean_temp, and error compare to output temperature curve
do we accept this scenario? perturb co2 and so2 separately, one bin at a time 
replica exchange parallelization =# 


