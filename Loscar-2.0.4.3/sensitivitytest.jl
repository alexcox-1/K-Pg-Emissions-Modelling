## do a loscar sensitivity test - ac.gr@ July 2022

using StatGeochem, Statistics, DelimitedFiles
    include("runloscar.jl")
    co2vals = readdlm("co2dist.csv",',');
    so2vals = readdlm("so2dist.csv",',');
    Expvals = readdlm("expdist.csv",',');
    carbvals = readdlm("carbdist.csv",',');
    Reminvals = readdlm("remindist.csv",',');
    doublevals = readdlm("doubledist.csv",',');

    rank = 2
    co2vals = co2vals[:,rank+1];
    svals = so2vals[:,rank+1];
    Expvals = Expvals[:,rank+1];
    Carbvals = carbvals[:,rank+1];
    Reminvals = Reminvals[:,rank+1];
    co2doublingrate = doublevals[rank+1];

    logco2vals = log.(co2vals);
    logsvals = log.(svals);
    logexpvals = log.(Expvals);
    logCarbvals = log.(Carbvals);
    bsrtemps = importdataset("tempdatabsr.csv",',');
    timev = bsrtemps["time"];
    temp = bsrtemps["temp"];
    temperror = bsrtemps["temperror"];
    timev .= (timev .- minimum(timev)) .* 1000000;

    tmv,pco2,loscartemp, d13csa, d13cba, d13csi, d13csp, d13cst = runloscar(timev,co2vals,svals,Expvals,co2doublingrate,Reminvals,Carbvals);

    plot(tmv,d13csa)