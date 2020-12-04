using StatGeochem, DelimitedFiles
d13cdata = importdataset("d13cdat.csv",',');
d13cages = d13cdata["age"];
d13cvals = d13cdata["d13c"];
d13csigma = zeros(length(d13cvals)) .+ 0.1;
relative_age_sigma = 0.02/100/2; # 0.5% 2-sigma
Age_sigma = d13cages * relative_age_sigma;
# the K/Pg age tie point from ODP 1262 
d13cages = 66.022 .- d13cages;
(c, m, e) = bin_bsr(d13cages, d13cvals, -0.5, 1, 300; x_sigma=Age_sigma, y_sigma=d13csigma, nresamplings=10000, p=1.0)
d13cbsr = [c m e];
writedlm("d13cdatabsr.csv",d13cbsr,',')
# add error to this in quadrature 