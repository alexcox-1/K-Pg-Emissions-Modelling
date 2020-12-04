# see what's up
using StatGeochem
temps2e4 = importdataset("d13c_temps.csv",',');
temps2e4 = temps2e4["temps"];
temps2e4 = reshape(temps2e4,300,102400);
co22e4 = importdataset("d13c_co2dist.csv",',');
co22e4 = co22e4["co2"];
co22e4 = reshape(co22e4,300,102400);
so22e4 = importdataset("d13c_so2dist.csv",',');
so22e4 = so22e4["so2"];
so22e4 = reshape(so22e4,300,102400);
ll2e4 = importdataset("d13c_ll.csv",',');
ll2e4 = ll2e4["ll"];
ll2e4 = reshape(ll2e4,200,512);
co2step = importdataset("d13c_co2step.csv",',');
co2step = co2step["co2step"];
co2step = reshape(co2step,200,512);
so2step = importdataset("d13c_so2step.csv",',');
so2step = so2step["so2step"];
so2step = reshape(so2step,200,512);
d13c = importdataset("d13c_d13c.csv",',');
d13c = d13c["d13c"];
d13c = reshape(d13c,300,102400);
