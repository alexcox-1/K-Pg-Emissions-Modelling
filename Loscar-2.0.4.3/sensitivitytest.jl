## do a loscar sensitivity test - ac.gr@ July 2022

using StatGeochem, Statistics
    include("runloscar.jl")
    co2vals = zeros(400) .+ 0.02;
    svals = zeros(400) .+ 0.01;
    Expvals = ones(400);
    Carbvals = ones(400);
    co2doublingrate = 3;
    bsrtemps = importdataset("tempdatabsr.csv",',');
    timev = bsrtemps["time"];
    temp = bsrtemps["temp"];
    temperror = bsrtemps["temperror"];
    timev .= (timev .- minimum(timev)) .* 1000000;
    Reminvals = ones(400);

    tmv,pco2,loscartemp, d13csa, d13cba, d13csi, d13csp, d13cst = runloscar(timev,co2vals,svals,Expvals,co2doublingrate,Reminvals,Carbvals);

    plot(tmv,d13csa)