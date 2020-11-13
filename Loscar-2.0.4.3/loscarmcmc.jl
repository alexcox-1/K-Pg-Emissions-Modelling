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
# we have to include a special one year bin for the 
# so-called impact event
timevals2 = 500002:5000:1500002;

timevals = [timevals1;timevals2];

# we now have a length 302 time array. let's fill it with
# co2 and so2 emissions.
# characteristic Pg/y will be 0.01 - 0.1
co2vals = zeros(302) .+ 0.01;
svals = zeros(302) .+ 0.01;

co2doublingrate = 3;

# do a loscar run!

tmv,pco2,loscartemp = runloscar(timevals,co2vals,svals,co2doublingrate);

# apply the sulfate correction
# convert svals to pinatubos a year
pinatubos = svals ./ 0.009;
sulfatecorr = -11.3 .* (1 .- exp.(-0.0466.*pinatubos));

loscartempwsulf = loscartemp .- sulfatecorr;




