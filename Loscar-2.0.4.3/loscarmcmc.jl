## do the loscar mcmc ac.gr@ 11/12/20
using StatGeochem

include("runloscar.jl")
include("fillnans.jl")

    bsrtemps = importdataset("tempdatabsr.csv",',');
    timev = bsrtemps["time"];
    temp = bsrtemps["temp"];
    temperror = bsrtemps["temperror"];
    # import the bootstrap resampled surface Atl d13C.
    bsrd13c = importdataset("d13cdatabsr.csv",',')
    d13cvals = bsrd13c["d13cval"];
    d13cerror = bsrd13c["d13cerror"]
    # add an error in quadrature given LOSCAR uncertainties (assumed 0.3)
    d13cerror .= sqrt.((d13cerror.^2) .+ 0.3^2)
    # add benthic d13c values
    bsrd13cb = importdataset("d13cdatabenthicbsr.csv",',')
    d13cbvals = bsrd13cb["d13cval"];
    d13cberror = bsrd13cb["d13cerror"]
    # add an error in quadrature
    d13cberror = sqrt.((d13cberror.^2) .+ 0.3^2)
    # get the time to start at zero, and be in years
    timev .= (timev .- minimum(timev)) .* 1000000;

    # in this time frame, time goes from 0 to 1500 kyr and
    # the kpg boundary is at 500 kyr.


    # we now have a length 400 time array. let's fill it with
    # co2 and so2 emissions.
    # characteristic Pg/y will be 0.01 - 0.1
    # change these to log
    co2vals = zeros(400) .+ 0.02;
    co2vals[301:400] .= 0.000001
    svals = zeros(400) .+ 0.01;
    svals[301:400] .= 0.000001
    # add the export reduction factor solved earlier
    expvals = ones(400);
    logco2vals = log.(co2vals);
    logsvals = log.(svals);
    logexpvals = log.(expvals);
    Reminvals = ones(400);
    Reminvals[1:200] .= 0.9995;
    Reminvals[201:400] .= 1.002;
    co2doublingrate = 3.0;
    Carbvals = ones(400);
    logCarbvals = log.(Carbvals);

    # do a loscar run!

    tmv,pco2,loscartemp, d13csa, d13cba = runloscar(timev,co2vals,svals,expvals,co2doublingrate,Reminvals,Carbvals);
    if isnan(tmv[1])
        ll = NaN
    else
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
        loscartempwsulf = loscartempwsulf[loscartimebin .<= 400];
        mu = Array{Float64,1}(undef,length(temp));
        nanbinmean!(mu,vec(tmv),loscartempwsulf,first(timev),last(timev),length(timev));
        # get rid of the internal NaNs
        mu = fillnans(mu,150);
        # get rid of the NaNs at the end
        if isnan(mu[400])
            mu[first(findall(x->isnan(x),mu)):400] .= mu[first(findall(x->isnan(x),mu))-1]
        end
        d13cmu = Array{Float64,1}(undef,length(temp));
        nanbinmean!(d13cmu,vec(tmv),d13csa,first(timev),last(timev),length(timev));
        d13cmu = fillnans(d13cmu,150);
        if isnan(d13cmu[400])
            d13cmu[first(findall(x->isnan(x),d13cmu)):400] .= d13cmu[first(findall(x->isnan(x),d13cmu))-1]
        end
        d13cbmu = Array{Float64,1}(undef,length(temp));
        nanbinmean!(d13cbmu,vec(tmv),d13cba,first(timev),last(timev),length(timev));
        d13cbmu = fillnans(d13cbmu,150);
        if isnan(d13cbmu[400])
            d13cbmu[first(findall(x->isnan(x),d13cbmu)):400] .= d13cbmu[first(findall(x->isnan(x),d13cbmu))-1]
        end
        ll = normpdf_ll(temp,temperror,mu) + normpdf_ll(d13cvals,d13cerror,d13cmu) + normpdf_ll(d13cbvals,d13cberror,d13cbmu) + normpdf_ll(3,0.1,co2doublingrate) + normpdf_ll(1,0.004,Reminvals) + normpdf_ll(1,1,exp.(logCarbvals));
    end
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
    logexpvalsᵣ = copy(logexpvals);
    Reminvalsᵣ = copy(Reminvals);
    logCarbvalsᵣ = copy(logCarbvals)
    co2doublingrateᵣ = copy(co2doublingrate);
    lldist = Array{Float64,1}(undef,numiter);
    doubledist = Array{Float64,1}(undef,numiter);
    remindist = Array{Float64,2}(undef,length(Reminvals),numiter);
    carbdist = Array{Float64,2}(undef,length(Carbvals),numiter);
    co2dist = Array{Float64,2}(undef,length(logco2vals),numiter);
    sdist = Array{Float64,2}(undef,length(logsvals),numiter);
    tempwsulfarray = Array{Float64,2}(undef,length(mu),numiter);
    expdist = Array{Float64,2}(undef,length(logexpvals),numiter);
    d13carray = Array{Float64,2}(undef,length(mu),numiter);
    d13cbarray = Array{Float64,2}(undef,length(mu),numiter);
    step_sigma_co2_array = Array{Float64,1}(undef,numiter);
    step_sigma_so2_array = Array{Float64,1}(undef,numiter);
    # do the mcmc
    # set the std of the proposal amplitude distribution
    co2_step_sigma = 0.1;
    so2_step_sigma = 0.1;
    exp_step_sigma = 0.01;
    remin_step_sigma = 0.0008;
    carb_step_sigma = 0.01;
    halfwidthc = 0.5;
    halfwidths = 0.5;
    halfwidthexp = 0.25;
    halfwidthcarb = 0.25;
    halfwidthremin = 0.25;
    counter = 0;
    for i = 1:numiter
        copyto!(logco2valsᵣ,logco2vals);
        copyto!(logsvalsᵣ,logsvals);
        copyto!(logexpvalsᵣ,logexpvals);
        copyto!(Reminvalsᵣ,Reminvals)
        copyto!(logCarbvalsᵣ,logCarbvals);
        co2doublingrateᵣ = co2doublingrate;
        randhalfwidth = halfwidthc * rand()*length(co2vals)

        randmu = rand()*length(co2vals)
        randamplitude = randn()*co2_step_sigma*2.9

        for j=1:length(co2vals)
            logco2valsᵣ[j] += randamplitude * ((randmu-randhalfwidth)<j<(randmu+randhalfwidth))

        end
        logco2valsᵣ[301:400] .= -20;
        randhalfwidths = halfwidths * rand()*length(co2vals)

        randmus = rand()*length(svals)

        randamplitudes = randn()*so2_step_sigma*2.9

        for j=1:length(svals)
            logsvalsᵣ[j] += randamplitudes * ((randmus-randhalfwidths)<j<(randmus+randhalfwidths))

        end
        logsvalsᵣ[301:400] .= -20;

        randhalfwidthexp = halfwidthexp * rand()*length(expvals)

        randmuexp = rand()*length(expvals)

        randamplitudeexp = randn()*exp_step_sigma*2.9

        for j=1:length(expvals)
            logexpvalsᵣ[j] += randamplitudeexp * ((randmuexp-randhalfwidthexp)<j<(randmuexp+randhalfwidthexp))

        end

        randhalfwidthremin = halfwidthremin * rand()*length(Reminvals)

        randmuremin = rand()*length(Reminvals)

        randamplituderemin = randn()*remin_step_sigma

        for j=1:length(Reminvals)
            Reminvalsᵣ[j] += randamplituderemin * ((randmuremin-randhalfwidthremin)<j<(randmuremin+randhalfwidthremin))

        end

        randhalfwidthcarb = halfwidthcarb * rand()*length(Carbvals)

        randmucarb = rand()*length(Carbvals)

        randamplitudecarb = randn()*carb_step_sigma

        for j=1:length(Carbvals)
            logCarbvalsᵣ[j] += randamplitudecarb * ((randmucarb-randhalfwidthcarb)<j<(randmucarb+randhalfwidthcarb))

        end
        # perturb the co2doubling rate normally 
        co2doublingrateᵣ += (randn() / 5)
        # run loscar with the new values
        tmv,pco2,loscartemp, d13csa, d13cba = runloscar(timev,exp.(logco2valsᵣ),exp.(logsvalsᵣ),exp.(logexpvalsᵣ),co2doublingrateᵣ,Reminvalsᵣ,exp.(logCarbvalsᵣ));

        pinatubos = exp.(logsvalsᵣ) ./ 0.009;
        sulfatecorr = -11.3 .* (1 .- exp.(-0.0466.*pinatubos));
        loscartimebin = Int.(floor.(tmv ./ 5000) .+ 1);
        loscartempwsulf = similar(loscartemp);
        for i = 1:length(loscartemp)
            loscartempwsulf[i] = loscartemp[i] .+ sulfatecorr[loscartimebin[i]];
        end
        # put the output temp into bins 

        # discard the last time val which is outside the range
        loscartempwsulf = loscartempwsulf[loscartimebin .<= 400];
        muᵣ = Array{Float64,1}(undef,length(temp));
        nanbinmean!(muᵣ,vec(tmv),loscartempwsulf,first(timev),last(timev),length(timev));
        # get rid of the internal NaNs
        muᵣ = fillnans(muᵣ,150);
        # get rid of the NaNs at the end
        if isnan(muᵣ[400])
            muᵣ[first(findall(x->isnan(x),muᵣ)):400] .= muᵣ[first(findall(x->isnan(x),muᵣ))-1]
        end
        d13cmuᵣ = Array{Float64,1}(undef,length(temp));
        nanbinmean!(d13cmuᵣ,vec(tmv),d13csa,first(timev),last(timev),length(timev));
        d13cmuᵣ = fillnans(d13cmuᵣ,150);
        if isnan(d13cmuᵣ[400])
            d13cmuᵣ[first(findall(x->isnan(x),d13cmuᵣ)):400] .= d13cmuᵣ[first(findall(x->isnan(x),d13cmuᵣ))-1]
        end
        d13cbmuᵣ = Array{Float64,1}(undef,length(temp));
        nanbinmean!(d13cbmuᵣ,vec(tmv),d13cba,first(timev),last(timev),length(timev));
        d13cbmuᵣ = fillnans(d13cbmuᵣ,150);
        if isnan(d13cbmuᵣ[400])
            d13cbmuᵣ[first(findall(x->isnan(x),d13cbmuᵣ)):400] .= d13cbmuᵣ[first(findall(x->isnan(x),d13cbmuᵣ))-1]
        end
        llᵣ = normpdf_ll(temp,temperror,muᵣ) + normpdf_ll(d13cvals,d13cerror,d13cmuᵣ) + normpdf_ll(d13cbvals,d13cberror,d13cbmuᵣ) + normpdf_ll(3,0.1,co2doublingrateᵣ) + normpdf_ll(1,0.004,Reminvalsᵣ) + normpdf_ll(1,1,exp.(logCarbvalsᵣ));

        # is this allowed?
        if log(rand()) < (llᵣ-ll)
            counter = 0
            ll = llᵣ
            logco2vals .= logco2valsᵣ  
            logsvals .= logsvalsᵣ  
            logexpvals .= logexpvalsᵣ
            co2doublingrate = co2doublingrateᵣ
            Reminvals = Reminvalsᵣ
            logCarbvals .= logCarbvalsᵣ
            co2_step_sigma = min(abs(randamplitude),1);
            so2_step_sigma = min(abs(randamplitudes),1)
            exp_step_sigma = min(abs(randamplitudeexp),0.1)
            remin_step_sigma = min(abs(randamplituderemin),0.005)
            carb_step_sigma = min(abs(randamplitudecarb),0.1)
            mu = muᵣ
            d13cmu = d13cmuᵣ
            d13cbmu = d13cbmuᵣ
        end
        lldist[i] = ll;
        co2dist[:,i] = logco2vals;
        sdist[:,i] = logsvals;
        expdist[:,i] = logexpvals;
        carbdist[:,i] = exp.(logCarbvals)
        tempwsulfarray[:,i] = mu;
        d13carray[:,i] = d13cmu;
    end
