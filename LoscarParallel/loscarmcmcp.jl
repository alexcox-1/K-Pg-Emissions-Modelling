## do the loscar mcmc ac.gr@ 11/12/20
using StatGeochem, MPI
MPI.Init()

include("runloscarp.jl")
include("fillnans.jl")

let
# Get MPI properties
    comm = MPI.COMM_WORLD
    ntasks = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    # Test stdout
    print("Hello from $rank of $ntasks processors!\n")
    # make scratch folder for each task
    prefix = "/dartfs-hpc/scratch/alex/loscar$rank"
    loscdir = "/dartfs-hpc/scratch/alex/K-Pg-Emissions-Modelling/LoscarParallel"
    system("rm -rf $prefix")
    system("mkdir -p $prefix")
    system("cp -r $loscdir $prefix")
    cd(prefix)
    bsrtemps = importdataset("LoscarParallel/tempdatabsr.csv",',');
    timev = bsrtemps["time"];
    temp = bsrtemps["temp"];
    temperror = bsrtemps["temperror"];
    # import the bootstrap resampled surface Atl d13C.
    bsrd13c = importdataset("LoscarParallel/d13cdatabsr.csv",',')
    d13cvals = bsrd13c["d13cval"];
    d13cerror = bsrd13c["d13cerror"]
    # add an error in quadrature given LOSCAR uncertainties (assumed 0.3)
    d13cerror .= sqrt.((d13cerror.^2) .+ 0.3^2)
    # add benthic d13c values
    bsrd13cb = importdataset("LoscarParallel/d13cdatabenthicbsr.csv",',')
    d13cbvals = bsrd13cb["d13cval"];
    d13cberror = bsrd13cb["d13cerror"]
    # add an error in quadrature
    d13cberror = sqrt.((d13cberror.^2) .+ 0.3^2)
    # get the time to start at zero, and be in years
    timev .= (timev .- minimum(timev)) .* 1000000;

    # in this time frame, time goes from 0 to 1500 kyr and
    # the kpg boundary is at 500 kyr.


    # we now have a length 300 time array. let's fill it with
    # co2 and so2 emissions.
    # characteristic Pg/y will be 0.01 - 0.1
    # change these to log
    co2vals = zeros(300) .+ 0.02;
    svals = zeros(300) .+ 0.01;
    # add the export reduction factor solved earlier
    expvals = ones(300);
    logco2vals = log.(co2vals);
    logsvals = log.(svals);
    logexpvals = log.(expvals);
    co2doublingrate = 3;

    # do a loscar run!

    tmv,pco2,loscartemp, d13csa, d13cba = runloscarp(timev,co2vals,svals,expvals,co2doublingrate);
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
        loscartempwsulf = loscartempwsulf[loscartimebin .<= 300];
        mu = Array{Float64,1}(undef,length(temp));
        nanbinmean!(mu,vec(tmv),loscartempwsulf,first(timev),last(timev),length(timev));
        # get rid of the internal NaNs
        mu = fillnans(mu,150);
        # get rid of the NaNs at the end
        if isnan(mu[300])
            mu[first(findall(x->isnan(x),mu)):300] .= mu[first(findall(x->isnan(x),mu))-1]
        end
        d13cmu = Array{Float64,1}(undef,length(temp));
        nanbinmean!(d13cmu,vec(tmv),d13csa,first(timev),last(timev),length(timev));
        d13cmu = fillnans(d13cmu,150);
        if isnan(d13cmu[300])
            d13cmu[first(findall(x->isnan(x),d13cmu)):300] .= d13cmu[first(findall(x->isnan(x),d13cmu))-1]
        end
        d13cbmu = Array{Float64,1}(undef,length(temp));
        nanbinmean!(d13cbmu,vec(tmv),d13cba,first(timev),last(timev),length(timev));
        d13cbmu = fillnans(d13cbmu,150);
        if isnan(d13cbmu[300])
            d13cbmu[first(findall(x->isnan(x),d13cbmu)):300] .= d13cbmu[first(findall(x->isnan(x),d13cbmu))-1]
        end
        ll = normpdf_ll(temp,temperror,mu) + normpdf_ll(d13cvals,d13cerror,d13cmu) + normpdf_ll(d13cbvals,d13cberror,d13cbmu);
    end
    numiter = 10;
    num_per_exchange = 1;
    ## monte carlo loop
    # perturb one of the co2 vals and one of the svals
    # work with co2 vals and svals in logspace
    # for parallel, do some number of iterations, collect answers, 
    # then continue 

    logco2valsᵣ = copy(logco2vals);
    logsvalsᵣ = copy(logsvals);
    logexpvalsᵣ = copy(logexpvals);
    lldist = Array{Float64,1}(undef,numiter);
    co2dist = Array{Float64,2}(undef,length(logco2vals),numiter);
    sdist = Array{Float64,2}(undef,length(logsvals),numiter);
    tempwsulfarray = Array{Float64,2}(undef,length(mu),numiter);
    expdist = Array{Float64,2}(undef,length(logexpvals),numiter);
    d13carray = Array{Float64,2}(undef,length(mu),numiter);
    d13cbarray = Array{Float64,2}(undef,length(mu),numiter);
    # create a record of what other MPI tasks have right now
    all_log_co2 = Array{Float64}(undef, length(logco2vals), ntasks);
    all_log_s = Array{Float64}(undef, length(logsvals), ntasks);
    all_log_exp = Array{Float64}(undef, length(logexpvals), ntasks);
    all_lls = Array{Float64}(undef,ntasks);
    step_sigma_co2_array = Array{Float64,1}(undef,numiter);
    step_sigma_so2_array = Array{Float64,1}(undef,numiter);
    # do the mcmc
    # set the std of the proposal amplitude distribution
    co2_step_sigma = 0.1;
    so2_step_sigma = 0.1;
    exp_step_sigma = 0.01;
    halfwidthc = 0.5;
    halfwidths = 0.5;
    halfwidthexp = 0.25;
    counter = 0;
    @inbounds for i = 1:numiter
        (rank == 0) && @warn "Iteration $i"
        print("Iteration $i")
        # update current prediction
        copyto!(logco2valsᵣ,logco2vals);
        copyto!(logsvalsᵣ,logsvals);
        copyto!(logexpvalsᵣ,logexpvals);
        # Exchange proposals, sometimes
        if i % num_per_exchange == 0
            # Exchange current proposals across all MPI tasks
            MPI.Allgather!(logco2vals, all_log_co2, comm)
            MPI.Allgather!(logsvals, all_log_s, comm)
            MPI.Allgather!(logexpvals, all_log_exp, comm)
            MPI.Allgather!(lldist[i:i], all_lls, comm)
            # Choose which proposal we want to adopt
            all_lls .-= maximum(all_lls) # rescale
            all_lls .= exp.(all_lls) # Convert to plain (relative) likelihoods
            lsum = sum(all_lls)
            chosen = 1
            cl = all_lls[chosen]
            r = rand()
            while r > (cl / lsum) && (chosen < length(all_lls))
                chosen += 1
                cl += all_lls[chosen]
            end
            logco2valsᵣ .= view(all_log_co2, :, chosen)
            logsvalsᵣ .= view(all_log_s, :, chosen)
            logexpvalsᵣ .= view(all_log_exp, :, chosen)
        end
        # choose which indices to perturb
        randhalfwidth = halfwidthc * rand()*length(co2vals)

        randmu = rand()*length(co2vals)
        randamplitude = randn()*co2_step_sigma*2.9

        for j=1:length(co2vals)
            logco2valsᵣ[j] += randamplitude * ((randmu-randhalfwidth)<j<(randmu+randhalfwidth))

        end
        randhalfwidths = halfwidths * rand()*length(co2vals)

        randmus = rand()*length(svals)

        randamplitudes = randn()*so2_step_sigma*2.9

        for j=1:length(svals)
            logsvalsᵣ[j] += randamplitudes * ((randmus-randhalfwidths)<j<(randmus+randhalfwidths))

        end

        randhalfwidthexp = halfwidthexp * rand()*length(expvals)

        randmuexp = rand()*length(expvals)

        randamplitudeexp = randn()*exp_step_sigma*2.9

        for j=1:length(expvals)
            logexpvalsᵣ[j] += randamplitudeexp * ((randmuexp-randhalfwidthexp)<j<(randmuexp+randhalfwidthexp))

        end

        # run loscar with the new values
        tmv,pco2,loscartemp, d13csa, d13cba = runloscarp(timev,exp.(logco2valsᵣ),exp.(logsvalsᵣ),exp.(logexpvalsᵣ),co2doublingrate);
        if all(isnan.(tmv))
            llᵣ = NaN
        else
        # do the sulfate corrections
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
            muᵣ = Array{Float64,1}(undef,length(temp));
            nanbinmean!(muᵣ,vec(tmv),loscartempwsulf,first(timev),last(timev),length(timev));
            # get rid of the internal NaNs
            muᵣ = fillnans(muᵣ,150);
            # get rid of the NaNs at the end
            if isnan(muᵣ[300])
                muᵣ[first(findall(x->isnan(x),muᵣ)):300] .= muᵣ[first(findall(x->isnan(x),muᵣ))-1]
            end
            d13cmuᵣ = Array{Float64,1}(undef,length(temp));
            nanbinmean!(d13cmuᵣ,vec(tmv),d13csa,first(timev),last(timev),length(timev));
            d13cmuᵣ = fillnans(d13cmuᵣ,150);
            if isnan(d13cmuᵣ[300])
                d13cmuᵣ[first(findall(x->isnan(x),d13cmuᵣ)):300] .= d13cmuᵣ[first(findall(x->isnan(x),d13cmuᵣ))-1]
            end
            d13cbmuᵣ = Array{Float64,1}(undef,length(temp));
            nanbinmean!(d13cbmuᵣ,vec(tmv),d13cba,first(timev),last(timev),length(timev));
            d13cbmuᵣ = fillnans(d13cbmuᵣ,150);
            if isnan(d13cbmuᵣ[300])
                d13cbmuᵣ[first(findall(x->isnan(x),d13cbmuᵣ)):300] .= d13cbmuᵣ[first(findall(x->isnan(x),d13cbmuᵣ))-1]
            end
            llᵣ = normpdf_ll(temp,temperror,muᵣ) + normpdf_ll(d13cvals,d13cerror,d13cmuᵣ) + normpdf_ll(d13cbvals,d13cberror,d13cbmuᵣ);
        end

        # is this allowed?
        if log(rand()) < (llᵣ-ll)
            counter = 0
            ll = llᵣ
            logco2vals .= logco2valsᵣ  
            logsvals .= logsvalsᵣ  
            logexpvals .= logexpvalsᵣ
            co2_step_sigma = min(abs(randamplitude),1);
            so2_step_sigma = min(abs(randamplitudes),1)
            exp_step_sigma = min(abs(randamplitudeexp),0.1)
            mu = muᵣ
            d13cmu = d13cmuᵣ
            d13cbmu = d13cbmuᵣ
        else
            counter += 1
        end
        if counter >= 8
            halfwidthc = max(halfwidthc*0.90,0.01);
            halfwidths = max(halfwidths*0.90,0.01);
            halfwidthexp = max(halfwidthexp*0.90,0.01);
        end
        # update the latest values
        lldist[i] = ll;
        co2dist[:,i] = logco2vals;
        sdist[:,i] = logsvals;
        expdist[:,i] = logexpvals;
        tempwsulfarray[:,i] = mu;
        d13carray[:,i] = d13cmu;
        d13cbarray[:,i] = d13cbmu;
        step_sigma_co2_array[i] = co2_step_sigma;
        step_sigma_so2_array[i] = so2_step_sigma;
    end

    # collate all the final values
    all_ll_dist = MPI.Gather(lldist, 0, comm)
    all_co2_dist = MPI.Gather(co2dist, 0, comm)
    all_s_dist = MPI.Gather(sdist, 0, comm)
    all_exp_dist = MPI.Gather(expdist,0,comm)
    all_temps = MPI.Gather(tempwsulfarray, 0, comm)
    all_d13c = MPI.Gather(d13carray,0,comm)
    all_d13cb = MPI.Gather(d13cbarray,0,comm)
    all_co2_step = MPI.Gather(step_sigma_co2_array,0,comm)
    all_so2_step = MPI.Gather(step_sigma_so2_array,0,comm)
    # write them to csv files
    if rank == 0
        writedlm("$loscdir/all_ll_dist.csv",all_ll_dist,',')
        writedlm("$loscdir/all_co2_dist.csv",all_co2_dist,',')
        writedlm("$loscdir/all_s_dist.csv",all_s_dist,',')
        writedlm("$loscdir/all_exp_dist.csv",all_exp_dist,',')
        writedlm("$loscdir/all_temps.csv",all_temps,',')
        # writedlm("$loscdir/all_co2_step.csv",all_co2_step,',')
        # writedlm("$loscdir/all_so2_step.csv",all_so2_step,',')
        writedlm("$loscdir/all_d13c.csv",all_d13c,',')
        writedlm("$loscdir/all_d13cb.csv",all_d13cb,',')
    end  
end       
MPI.Finalize()
