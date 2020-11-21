# see what's up
using StatGeochem
temps2e4 = importdataset("5e4_temps.csv",',');
temps2e4 = temps2e4["temps"];
temps2e4 = reshape(temps2e4,300,51200);
co22e4 = importdataset("5e4_co2.csv",',');
co22e4 = co22e4["co2"];
co22e4 = reshape(co22e4,300,51200);
so22e4 = importdataset("5e4_so2.csv",',');
so22e4 = so22e4["so2"];
so22e4 = reshape(so22e4,300,51200);
ll2e4 = importdataset("5e4_ll.csv",',');
ll2e4 = ll2e4["ll"];
