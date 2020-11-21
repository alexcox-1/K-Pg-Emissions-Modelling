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
    system("mkdir -p $prefix")
    system("cp -r $loscdir $prefix")
    cd(prefix)
    bsrtemps = importdataset("LoscarParallel/tempdatabsr.csv",',');
    timev = bsrtemps["time"];
    temp = bsrtemps["temp"];
    temperror = bsrtemps["temperror"];

    # get the time to start at zero, and be in years
    timev .= (timev .- minimum(timev)) .* 1000000;

    # in this time frame, time goes from 0 to 1500 kyr and
    # the kpg boundary is at 500 kyr.


    # we now have a length 300 time array. let's fill it with
    # co2 and so2 emissions.
    # characteristic Pg/y will be 0.01 - 0.1
    # change these to log
    co2vals = zeros(300) .+ 0.04;
    svals = zeros(300) .+ 0.01;
    logco2vals = log.(co2vals);
    logsvals = log.(svals);

    co2doublingrate = 3;

    # do a loscar run!

    tmv,pco2,loscartemp = runloscarp(timev,co2vals,svals,co2doublingrate);
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
        nanmean!(mu,vec(tmv),loscartempwsulf,first(timev),last(timev),length(timev));
        mu = fillnans(mu,50);

        ll = normpdf_ll(temp,temperror,mu);
    end
    numiter = 250;
    num_per_exchange = 1;
    ## monte carlo loop
    # perturb one of the co2 vals and one of the svals
    # work with co2 vals and svals in logspace
    # for parallel, do some number of iterations, collect answers, 
    # then continue 

    logco2valsᵣ = copy(logco2vals);
    logsvalsᵣ = copy(logsvals);
    lldist = Array{Float64,1}(undef,numiter);
    co2dist = Array{Float64,2}(undef,length(logco2vals),numiter);
    sdist = Array{Float64,2}(undef,length(logsvals),numiter);
    tempwsulfarray = Array{Float64,2}(undef,length(mu),numiter);
    # create a record of what other MPI tasks have right now
    all_log_co2 = Array{Float64}(undef, length(logco2vals), ntasks);
    all_log_s = Array{Float64}(undef, length(logsvals), ntasks);
    all_lls = Array{Float64}(undef,ntasks);
    # do the mcmc
    # set the std of the proposal amplitude distribution
    co2_step_sigma = 0.15;
    so2_step_sigma = 0.15;
    @inbounds for i = 1:numiter
        # update current prediction
        copyto!(logco2valsᵣ,logco2vals);
        copyto!(logsvalsᵣ,logsvals);
        # Exchange proposals, sometimes
        if i % num_per_exchange == 0
            # Exchange current proposals across all MPI tasks
            MPI.Allgather!(logco2vals, all_log_co2, comm)
            MPI.Allgather!(logsvals, all_log_s, comm)
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
        end
        # choose which indices to perturb
        randsigma = rand()*length(co2vals)/100
        randsigma2 = rand()*length(co2vals)/100
        randsigma3 = rand()*length(co2vals)/100
        randmu = rand()*length(co2vals)
        randmu2 = rand()*length(co2vals)
        randmu3 = rand()*length(co2vals)
        randamplitude = randn()*2.9*co2_step_sigma
        randamplitude2 = randn()*2.9*co2_step_sigma
        randamplitude3 = randn()*2.9*co2_step_sigma
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
        randamplitudes = randn()*2.9*so2_step_sigma # multiplied by some value related to the last perturbation (2.9*last amplitude)
        randamplitude2s = randn()*2.9*so2_step_sigma
        randamplitude3s = randn()*2.9*so2_step_sigma
        for j=1:length(svals)
            logsvalsᵣ[j] += randamplitudes * normpdf(randmus, randsigmas, j)
            logsvalsᵣ[j] += randamplitude2s * normpdf(randmu2s, randsigma2s, j)
            logsvalsᵣ[j] += randamplitude3s * normpdf(randmu3s, randsigma3s, j)
        end
        # run loscar with the new values
        tmv,pco2,loscartemp = runloscarp(timev,exp.(logco2valsᵣ),exp.(logsvalsᵣ),co2doublingrate);
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
            mu = Array{Float64,1}(undef,length(temp));
            nanmean!(mu,vec(tmv),loscartempwsulf,first(timev),last(timev),length(timev));
            mu = fillnans(mu,5);

            llᵣ = normpdf_ll(temp,temperror,mu);
        end

        # is this allowed?
        if log(rand()) < (llᵣ-ll)
            ll = llᵣ
            logco2vals .= logco2valsᵣ  
            logsvals .= logsvalsᵣ  
            co2_step_sigma = (randamplitude + randamplitude2 + randamplitude3) / 3;
            so2_step_sigma = (randamplitudes + randamplitudes2 + randamplitudes3) / 3;
        end
        # update the latest values
        lldist[i] = ll;
        co2dist[:,i] = logco2vals;
        sdist[:,i] = logsvals;
        tempwsulfarray[:,i] = mu;

    end

    all_ll_dist = MPI.Gather(lldist, 0, comm)
    all_co2_dist = MPI.Gather(co2dist, 0, comm)
    all_s_dist = MPI.Gather(sdist, 0, comm)
    all_temps = MPI.Gather(tempwsulfarray, 0, comm)

    if rank == 0
        writedlm("$loscdir/all_ll_dist.csv",all_ll_dist,',')
        writedlm("$loscdir/all_co2_dist.csv",all_co2_dist,',')
        writedlm("$loscdir/all_s_dist.csv",all_s_dist,',')
        writedlm("$loscdir/all_temps.csv",all_temps,',')
    end  
end       
MPI.Finalize()