## do a loscar sensitivity test - ac.gr@ July 2022

using StatGeochem, Statistics
    include("runloscar.jl")
    co2vals = zeros(400);
    svals = zeros(400);
    Expvals = ones(400);
    co2doublingrate = 3;
    bsrtemps = importdataset("tempdatabsr.csv",',');
    timev = bsrtemps["time"];
    temp = bsrtemps["temp"];
    temperror = bsrtemps["temperror"];
    timev .= (timev .- minimum(timev)) .* 1000000;
    Reminvals = ones(400);
    Reminvals[1:200] .= 0.9995
    Reminvals[201:400] .= 1.002

    tmv,pco2,loscartemp, d13csa, d13cba = runloscar(timev,co2vals,svals,Expvals,co2doublingrate,Reminvals);

    plot(tmv,d13csa)