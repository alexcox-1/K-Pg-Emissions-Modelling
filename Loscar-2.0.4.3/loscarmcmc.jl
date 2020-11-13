## do the loscar mcmc ac.gr@ 11/12/20
using StatGeochem
bsrtemps = importdataset("tempdatabsr.csv",',');
time = bsrtemps["time"];
temp = bsrtemps["temp"];
temperror = bsrtemps["error"];

# get the time to start at zero, and be in years
time .= (time .- minimum(time)) .* 1000000;

# in this time frame, time goes from 0 to 1500 kyr and
# the kpg boundary is at 500 kyr.

timevals1 = 1:5000:500001;
timevals2 = 500002:5000:1500002;

timevals = [timevals1;timevals2];


