# given a converged solution, explore the solution space on 1024 cores
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
        # get the time to start at zero, and be in years
        timev .= (timev .- minimum(timev)) .* 1000000;
    
        # in this time frame, time goes from 0 to 1500 kyr and
        # the kpg boundary is at 500 kyr.
    
    
        # we now have a length 300 time array. let's fill it with
        # co2 and so2 emissions.
        # characteristic Pg/y will be 0.01 - 0.1
        # change these to log
        co2vals = readdlm("co2_soln.csv");
        svals = readdlm("s_soln.csv");
        logco2vals = log.(co2vals);
        logsvals = log.(svals);

        o2doublingrate = 3;

    # do a loscar run!

    tmv,pco2,loscartemp, d13csa = runloscarp(timev,co2vals,svals,co2doublingrate);
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
        d13cmu = Array{Float64,1}(undef,length(temp));
        nanmean!(d13cmu,vec(tmv),d13csa,first(timev),last(timev),length(timev));
        d13cmu = fillnans(d13cmu,50);
        ll = normpdf_ll(temp,temperror,mu) + normpdf_ll(d13cvals,d13cerror,d13cmu);
    end

    numiter = 10;

    logco2valsᵣ = copy(logco2vals);
    logsvalsᵣ = copy(logsvals);
    lldist = Array{Float64,1}(undef,numiter);
    co2dist = Array{Float64,2}(undef,length(logco2vals),numiter);
    sdist = Array{Float64,2}(undef,length(logsvals),numiter);
    tempwsulfarray = Array{Float64,2}(undef,length(mu),numiter);
    d13carray = Array{Float64,2}(undef,length(mu),numiter);
    # create a record of what other MPI tasks have right now
    all_log_co2 = Array{Float64}(undef, length(logco2vals), ntasks);
    all_log_s = Array{Float64}(undef, length(logsvals), ntasks);
    all_lls = Array{Float64}(undef,ntasks);

    @inbounds for i = 1:numiter
        (rank == 0) && @warn "Iteration $i"
        print("Iteration $i")
        # update current prediction
        copyto!(logco2valsᵣ,logco2vals);
        copyto!(logsvalsᵣ,logsvals);
        # randomly adjust a few co2 and s values by 0.1 log units
        randamplitude = 0
        randamplitudes = 0
             # modify co2 vals
        randamplitude = randn()*0.1
            for j=rand(1:300,3)
                 logco2valsᵣ[j] += randamplitude 
            end
         # modify s vals
        randamplitudes = randn()*0.1
            for j=rand(1:300,3)
                 logsvalsᵣ[j] += randamplitudes 
            end

        tmv,pco2,loscartemp, d13csa = runloscarp(timev,exp.(logco2valsᵣ),exp.(logsvalsᵣ),co2doublingrate);
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
            nanmean!(muᵣ,vec(tmv),loscartempwsulf,first(timev),last(timev),length(timev));
            muᵣ = fillnans(muᵣ,5);
            # include log likelihood of d13C
            d13cmuᵣ = Array{Float64,1}(undef,length(temp));
            nanmean!(d13cmuᵣ,vec(tmv),d13csa,first(timev),last(timev),length(timev));
            d13cmuᵣ = fillnans(d13cmuᵣ,50);
            llᵣ = normpdf_ll(temp,temperror,muᵣ) + normpdf_ll(d13cvals,d13cerror,d13cmuᵣ);
        end
        if log(rand()) < (llᵣ-ll)
            counter = 0
            ll = llᵣ
            logco2vals .= logco2valsᵣ  
            logsvals .= logsvalsᵣ  
            co2_step_sigma = min(abs(randamplitude),1);
            so2_step_sigma = min(abs(randamplitudes),1)
            mu = muᵣ
            d13cmu = d13cmuᵣ
        end
        # update the values
        lldist[i] = ll;
        co2dist[:,i] = logco2vals;
        sdist[:,i] = logsvals;
        tempwsulfarray[:,i] = mu;
        d13carray[:,i] = d13cmu;
    end

        # collate all the final values
        all_ll_dist = MPI.Gather(lldist, 0, comm)
        all_co2_dist = MPI.Gather(co2dist, 0, comm)
        all_s_dist = MPI.Gather(sdist, 0, comm)
        all_temps = MPI.Gather(tempwsulfarray, 0, comm)
        all_d13c = MPI.Gather(d13carray,0,comm)

    if rank == 0
        writedlm("$loscdir/ll_dist.csv",all_ll_dist,',')
        writedlm("$loscdir/co2_dist.csv",all_co2_dist,',')
        writedlm("$loscdir/s_dist.csv",all_s_dist,',')
        writedlm("$loscdir/temps.csv",all_temps,',')
        writedlm("$loscdir/d13c.csv",all_d13c,',')
    end  
    
end
MPI.Finalize()