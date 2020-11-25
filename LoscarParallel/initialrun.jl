# see what's up
using StatGeochem
temps2e4 = importdataset("gauss_100_temps.csv",',');
temps2e4 = temps2e4["temps"];
temps2e4 = reshape(temps2e4,300,12800);
co22e4 = importdataset("gauss_100_co2vals.csv",',');
co22e4 = co22e4["co2"];
co22e4 = reshape(co22e4,300,12800);
so22e4 = importdataset("gauss_100_so2vals.csv",',');
so22e4 = so22e4["so2"];
so22e4 = reshape(so22e4,300,12800);
ll2e4 = importdataset("gauss_100_ll.csv",',');
ll2e4 = ll2e4["ll"];
ll2e4 = reshape(ll2e4,100,128);
co2step = importdataset("gauss_100_co2step.csv",',');
co2step = co2step["co2step"];
co2step = reshape(co2step,100,128);
so2step = importdataset("gauss_100_so2step.csv",',');
so2step = so2step["so2step"];
so2step = reshape(so2step,100,128);
