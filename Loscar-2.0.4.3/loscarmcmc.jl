## do the loscar mcmc ac.gr@ 11/12/20
using StatGeochem
    include("runloscar.jl")
    bsrtemps = importdataset("tempdatabsr.csv",',');
    timev = bsrtemps["time"];
    temp = bsrtemps["temp"];
    temperror = bsrtemps["temperror"];

    # get the time to start at zero, and be in years
    timev .= (timev .- minimum(timev)) .* 1000000;

    # in this time frame, time goes from 0 to 1500 kyr and
    # the kpg boundary is at 500 kyr.

    # we have to include a special one year bin for the 
    # so-called impact event


    # we now have a length 302 time array. let's fill it with
    # co2 and so2 emissions.
    # characteristic Pg/y will be 0.01 - 0.1
    # change these to log
    co2vals = zeros(300) .+ 0.04; 
    svals = zeros(300) .+ 0.01;
    logco2vals = log.(co2vals);
    logsvals = log.(svals);

    co2doublingrate = 3;

    # do a loscar run!

    tmv,pco2,loscartemp = runloscar(timev,co2vals,svals,co2doublingrate);

    # apply the sulfate correction
    # convert svals to pinatubos a year
    pinatubos = svals ./ 0.009;
    sulfatecorr = -11.3 .* (1 .- exp.(-0.0466.*pinatubos));
    loscartimebin = Int.(floor.(tmv ./ 5000) .+ 1);
    loscartempwsulf = similar(loscartemp);
    for i = 1:length(loscartemp)
        loscartempwsulf[i] = loscartemp[i] .+ sulfatecorr[loscartimebin[i]];
    end
    # put the output temp into bins 

    # discard the last time val which is outside the range
    loscartempwsulf = loscartempwsulf[loscartimebin .<= 300];
    mu = Array{Float64,1}(undef,length(temp));
    nanmean!(mu,vec(tmv),loscartempwsulf,first(timev),last(timev),length(timev));
    mu = fillnans(mu,50);

    ll = normpdf_ll(temp,temperror,mu);

    numiter = 5;

    ## monte carlo loop
    # perturb one of the co2 vals and one of the svals
    # work with co2 vals and svals in logspace
    # randomindex = rand(1:length(logco2vals));
    # logco2vals[randomindex] += randn()/100;
    # logsvals[randomindex] +=  randn()/100;
    # careful to not overallocate
    # svals .= exp.(logsvals)
    # runloscar(......svals.....)
    # save the log ones...? porque no los dos
    # accept or reject proposal based on ll
    # also have llprop
    logco2valsᵣ = copy(logco2vals);
    logsvalsᵣ = copy(logsvals);
    lldist = Array{Float64,1}(undef,numiter);
    co2dist = Array{Float64,2}(undef,length(logco2vals),numiter);
    sdist = Array{Float64,2}(undef,length(logsvals),numiter);
    tempwsulfarray = Array{Float64,2}(undef,length(mu),numiter);
    co2_step_sigma = 0.1;
    so2_step_sigma = 0.1;
    for i = 1:numiter
        copyto!(logco2valsᵣ,logco2vals);
        copyto!(logsvalsᵣ,logsvals);
        randsigma = rand()*length(co2vals)/100
        randsigma2 = rand()*length(co2vals)/100
        randsigma3 = rand()*length(co2vals)/100
        randmu = rand()*length(co2vals)
        randmu2 = rand()*length(co2vals)
        randmu3 = rand()*length(co2vals)
        randamplitude = randn()*co2_step_sigma*2.9
        randamplitude2 = randn()*co2_step_sigma*2.9
        randamplitude3 = randn()*co2_step_sigma*2.9
        for j=1:length(co2vals)
            logco2valsᵣ[j] += randamplitude * normpdf(randmu, randsigma, j)
            logco2valsᵣ[j] += randamplitude2 * normpdf(randmu2, randsigma2, j)
            logco2valsᵣ[j] += randamplitude3 * normpdf(randmu3, randsigma3, j)
        end
        randsigmas = rand()*length(svals)/100
        randsigma2s = rand()*length(svals)/100
        randsigma3s = rand()*length(svals)/100
        randmus = rand()*length(svals)
        randmu2s = rand()*length(svals)
        randmu3s = rand()*length(svals)
        randamplitudes = randn()*so2_step_sigma*2.9
        randamplitude2s = randn()*so2_step_sigma*2.9
        randamplitude3s = randn()*so2_step_sigma*2.9
        for j=1:length(svals)
            logsvalsᵣ[j] += randamplitudes * normpdf(randmus, randsigmas, j)
            logsvalsᵣ[j] += randamplitude2s * normpdf(randmu2s, randsigma2s, j)
            logsvalsᵣ[j] += randamplitude3s * normpdf(randmu3s, randsigma3s, j)
        end

        tmv,pco2,loscartemp = runloscar(timev,exp.(logco2valsᵣ),exp.(logsvalsᵣ),co2doublingrate);

        pinatubos = exp.(logsvalsᵣ) ./ 0.009;
        sulfatecorr = -11.3 .* (1 .- exp.(-0.0466.*pinatubos));
        loscartimebin = Int.(floor.(tmv ./ 5000) .+ 1);
        loscartempwsulf = similar(loscartemp);
        for i = 1:length(loscartemp)
            loscartempwsulf[i] = loscartemp[i] .+ sulfatecorr[loscartimebin[i]];
        end
        # put the output temp into bins 

        # discard the last time val which is outside the range
        loscartempwsulf = loscartempwsulf[loscartimebin .<= 300];
        mu = Array{Float64,1}(undef,length(temp));
        nanmean!(mu,vec(tmv),loscartempwsulf,first(timev),last(timev),length(timev));
        mu = fillnans(mu,5);

        llᵣ = normpdf_ll(temp,temperror,mu);

        # is this allowed?
        if log(rand()) < (llᵣ-ll)
            ll = llᵣ
            logco2vals .= logco2valsᵣ  
            logsvals .= logsvalsᵣ  
            co2_step_sigma = (randamplitude+randamplitude2+randamplitude3)/3
            so2_step_sigma = (randamplitudes+randamplitude2s+randamplitude3s)/3
        end
        lldist[i] = ll;
        co2dist[:,i] = logco2vals;
        sdist[:,i] = logsvals;
        tempwsulfarray[:,i] = mu;

    end
