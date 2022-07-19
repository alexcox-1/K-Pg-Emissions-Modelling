using StatGeochem, DelimitedFiles
d13cdata = importdataset("d13cdata.csv",',');
d13cages = d13cdata["age"];
d13cvals = d13cdata["d13"];
d13csigma = zeros(length(d13cvals)) .+ 0.1;
relative_age_sigma = 0.02/100/2; # 0.5% 2-sigma
Age_sigma = d13cages * relative_age_sigma;
# the K/Pg age tie point from ODP 1262 
d13cages = 66.022 .- d13cages;
(c, m, e) = bin_bsr(d13cages, d13cvals, -1, 1, 400; x_sigma=Age_sigma, y_sigma=d13csigma, nresamplings=10000, p=1.0)
d13cbsr = [c m e];
writedlm("d13cdatabsr.csv",d13cbsr,',')
# do the same for benthic values
d13cdata = importdataset("d13cdat.csv",',');
d13cbenthicages = d13cdata["age"][1:780];
d13cbenthicvals = d13cdata["d13c"][1:780];
d13cbenthicsigma = zeros(length(d13cbenthicvals)) .+ 0.1;
relative_benthic_age_sigma = 0.02/100/2; # 0.5% 2-sigma
Age_benthic_sigma = d13cbenthicages * relative_benthic_age_sigma;
# the K/Pg age tie point from ODP 1262 
d13cbenthicages = 66.022 .- d13cbenthicages;
(c, m, e) = bin_bsr(d13cbenthicages, d13cbenthicvals, -1, 1, 400; x_sigma=Age_benthic_sigma, y_sigma=d13cbenthicsigma, nresamplings=10000, p=1.0)
d13cbenthicbsr = [c m e];
writedlm("d13cdatabenthicbsr.csv",d13cbenthicbsr,',')
e = sqrt(e.^2 .+ 0.3^2);
