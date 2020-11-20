# see what's up
using StatGeochem
temps1e4 = importdataset("1e5_temps.csv",',');
temps1e4 = temps1e4["temps"];
temps1e4 = reshape(temps1e4,300,1600);
co21e4 = importdataset("1e5_co2.csv",',');
co21e4 = co21e4["co2"];
co21e4 = reshape(co21e4,300,1600);
so21e4 = importdataset("1e5_so2.csv",',');
so21e4 = so21e4["so2"];
so21e4 = reshape(so21e4,300,1600);
ll1e4 = importdataset("1e5_ll.csv",',');
ll1e4 = ll1e4["ll"];
