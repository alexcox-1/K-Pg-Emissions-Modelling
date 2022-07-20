## do a loscar sensitivity test - ac.gr@ July 2022

using StatGeochem, Statistics
    include("runloscar.jl")
    co2vals = zeros(400);
    co2vals[101:200] .= 0.01;
    co2vals[201:300] .= 0.04;
    co2vals[301:400] .= 0.0;
    svals = zeros(400) .+ 0.01;
    Expvals = ones(400);
    co2doublingrate = 3;
    bsrtemps = importdataset("tempdatabsr.csv",',');
    timev = bsrtemps["time"];
    temp = bsrtemps["temp"];
    temperror = bsrtemps["temperror"];
    timev .= (timev .- minimum(timev)) .* 1000000;

    tmv,pco2,loscartemp, d13csa, d13cba = runloscar(timev,co2vals,svals,Expvals,co2doublingrate);

    plot(tmv,d13csa)