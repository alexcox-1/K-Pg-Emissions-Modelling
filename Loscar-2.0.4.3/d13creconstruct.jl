d13cdata = importdataset("d13cdat.csv",',');
d13cages = d13cdata["age"];
d13cvals = d13cdata["d13c"];
d13csigma = zeros(length(d13cvals)) .+ 0.05;
relative_age_sigma = 0.05/100/2; # 0.5% 2-sigma
Age_sigma = d13cages * relative_age_sigma;
d13cages = d13cages .- 66.02;
(c, m, e) = bin_bsr(d13cages, d13cvals, -0.5, 1, 300; x_sigma=Age_sigma, y_sigma=d13csigma, nresamplings=10000, p=1.0)
p3 = plot(c, m, seriestype=:scatter, xflip=true, label="", framestyle=:box,ribbon=(e),size=(900,300),xtickfont=14,ytickfont=14)
savefig(p3,"/Users/alexcox/Documents/manuscripts/loscarletter/d13c.pdf")