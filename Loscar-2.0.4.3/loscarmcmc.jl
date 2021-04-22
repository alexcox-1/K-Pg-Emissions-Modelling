## do the loscar mcmc ac.gr@ 11/12/20
using StatGeochem
    include("runloscar.jl")
    include("fillnans.jl")
    bsrtemps = importdataset("tempdatabsr.csv",',');
    timev = bsrtemps["time"];
    temp = bsrtemps["temp"];
    temperror = bsrtemps["temperror"];
    bsrd13c = importdataset("d13cdatabsr.csv",',')
    d13cvals = bsrd13c["d13cval"];
    d13cerror = bsrd13c["d13cerror"]
    bsrd13cb = importdataset("d13cdatabenthicbsr.csv",',')
    d13cbvals = bsrd13cb["d13cval"];
    d13cberror = bsrd13cb["d13cerror"]
    # add an error in quadrature given LOSCAR uncertainties (assumed 0.5)
    d13cerror .= sqrt.((d13cerror.^2) .+ 0.3^2)
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
    co2vals = zeros(300) .+ 0.01; 
    svals = zeros(300);
    Expvals = ones(300) .- 0.3;
    logco2vals = log.(co2vals);
    logsvals = log.(svals);

    co2doublingrate = 3;

    # do a loscar run!

    tmv,pco2,loscartemp, d13csa, d13cba = runloscar(timev,co2vals,svals,Expvals,co2doublingrate);

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
    d13cmu = Array{Float64,1}(undef,length(temp));
    nanmean!(d13cmu,vec(tmv),d13csa,first(timev),last(timev),length(timev));
    d13cmu = fillnans(mu,50);
    ll = normpdf_ll(temp,temperror,mu) + normpdf_ll(d13cvals,d13cerror,d13cmu);

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
    d13carray = Array{Float64,2}(undef,length(mu),numiter);
    co2_step_sigma = 0.1;
    so2_step_sigma = 0.1;
    for i = 1:numiter
        copyto!(logco2valsᵣ,logco2vals);
        copyto!(logsvalsᵣ,logsvals);
        randhalfwidth = rand()*length(co2vals)/100

        randmu = rand()*length(co2vals)
        randamplitude = randn()*co2_step_sigma*2.9

        for j=1:length(co2vals)
            logco2valsᵣ[j] += randamplitude * ((randmu-randhalfwidth)<j<(randmu+randhalfwidth))

        end
        randhalfwidths = rand()*length(svals)/100

        randmus = rand()*length(svals)

        randamplitudes = randn()*so2_step_sigma*2.9

        for j=1:length(svals)
            logsvalsᵣ[j] += randamplitudes * ((randmus-randhalfwidths)<j<(randmus+randhalfwidths))

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
        d13cmu = Array{Float64,1}(undef,length(temp));
        nanmean!(d13cmu,vec(tmv),d13csa,first(timev),last(timev),length(timev));
        d13cmu = fillnans(d13cmu,50);
        llᵣ = normpdf_ll(temp,temperror,mu) + normpdf_ll(d13cvals,d13cerror,d13cmu);

        # is this allowed?
        if log(rand()) < (llᵣ-ll)
            ll = llᵣ
            logco2vals .= logco2valsᵣ  
            logsvals .= logsvalsᵣ  
            co2_step_sigma = min(abs(randamplitude),1)
            so2_step_sigma = min(abs(randamplitudes),1)
        end
        lldist[i] = ll;
        co2dist[:,i] = logco2vals;
        sdist[:,i] = logsvals;
        tempwsulfarray[:,i] = mu;
        d13carray[:,i] = d13cmu;
    end
