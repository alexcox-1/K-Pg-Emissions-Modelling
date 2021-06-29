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
    # add benthic d13c values
    bsrd13cb = importdataset("LoscarParallel/d13cdatabenthicbsr.csv",',')
    d13cbvals = bsrd13cb["d13cval"];
    d13cberror = bsrd13cb["d13cerror"]
    # add an error in quadrature
    d13cberror = sqrt.((d13cberror.^2) .+ 0.3^2)
    # get the time to start at zero, and be in years
    timev .= (timev .- minimum(timev)) .* 1000000;

    expvals = ones(300);
    logexpvals = log.(expvals);
    co2doublingrate = 3;

    # do a loscar run!

    tmv,pco2,loscartemp, d13csa, d13cba = runloscarp(timev,zeros(300),zeros(300),expvals,co2doublingrate);

    if isnan(tmv[1])
        ll = NaN
    else

        d13cbmu = Array{Float64,1}(undef,300);
        nanmean!(d13cbmu,vec(tmv),d13cba,first(timev),last(timev),length(timev));
        d13cbmu = fillnans(d13cbmu,150);
        if isnan(d13cbmu[300])
            d13cbmu[first(findall(x->isnan(x),d13cbmu)):300] .= d13cbmu[first(findall(x->isnan(x),d13cbmu))-1]
        end

        ll = normpdf_ll(d13cbvals,d13cberror,d13cbmu);

    end

    numiter = 175;
    num_per_exchange = 1;
    ## monte carlo loop
    # perturb one of the co2 vals and one of the svals
    # work with co2 vals and svals in logspace
    # for parallel, do some number of iterations, collect answers, 
    # then continue 
    logexpvalsᵣ = copy(logexpvals);
    lldist = Array{Float64,1}(undef,numiter);
    expdist = Array{Float64,2}(undef,length(logexpvals),numiter);
    d13cbarray = Array{Float64,2}(undef,300,numiter);
    all_log_exp = Array{Float64}(undef, length(logexpvals), ntasks);
    all_lls = Array{Float64}(undef,ntasks);
    exp_step_sigma = 0.01;
    halfwidthexp = 0.25;
    counter = 0;
    @inbounds for i = 1:1
        (rank == 0) && @warn "Iteration $i"
        print("Iteration $i")
        copyto!(logexpvalsᵣ,logexpvals);
        if i % num_per_exchange == 0
            # Exchange current proposals across all MPI tasks
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
            logexpvalsᵣ .= view(all_log_exp, :, chosen)
        end

        randhalfwidthexp = halfwidthexp * rand()*length(expvals)

        randmuexp = rand()*length(expvals)

        randamplitudeexp = randn()*exp_step_sigma*2.9

        for j=1:length(expvals)
            logexpvalsᵣ[j] += randamplitudeexp * ((randmuexp-randhalfwidthexp)<j<(randmuexp+randhalfwidthexp))

        end

        logexpvalsᵣ[logexpvalsᵣ .> 0.47] .= 0.47
        logexpvalsᵣ[logexpvalsᵣ .< -1.386] .= -1.386

        tmv,pco2,loscartemp, d13csa, d13cba = runloscarp(timev,zeros(300),zeros(300),exp.(logexpvalsᵣ),co2doublingrate);

        if all(isnan.(tmv))
            llᵣ = NaN
        else

            d13cbmuᵣ = Array{Float64,1}(undef,300);
            nanmean!(d13cbmuᵣ,vec(tmv),d13cba,first(timev),last(timev),length(timev));
            d13cbmuᵣ = fillnans(d13cbmuᵣ,150);
            if isnan(d13cbmuᵣ[300])
                d13cbmuᵣ[first(findall(x->isnan(x),d13cbmuᵣ)):300] .= d13cbmuᵣ[first(findall(x->isnan(x),d13cbmuᵣ))-1]
            end
            llᵣ = normpdf_ll(d13cbvals,d13cberror,d13cbmuᵣ);
        end

        if log(rand()) < (llᵣ-ll)
            counter = 0
            ll = llᵣ
            logexpvals .= logexpvalsᵣ
            exp_step_sigma = min(abs(randamplitudeexp),0.1)
            d13cbmu = d13cbmuᵣ
        else
            counter += 1
        end
        if counter >= 8
            halfwidthexp = max(halfwidthexp*0.90,0.01);
        end
        # update the latest values
        lldist[i] = ll;
        expdist[:,i] = logexpvals;
        d13cbarray[:,i] = d13cbmu;
    end

    @inbounds for i = 2:175
        (rank == 0) && @warn "Iteration $i"
        print("Iteration $i")
        copyto!(logexpvalsᵣ,logexpvals);
        # perturb three random indices

        changeindices = rand(1:300,3);

        logexpvalsᵣ[changeindices] .+= (randn(3) .* 2.9 .* exp_step_sigma)

        logexpvalsᵣ[logexpvalsᵣ .> 0.47] .= 0.47
        logexpvalsᵣ[logexpvalsᵣ .< -1.386] .= -1.386

        tmv,pco2,loscartemp, d13csa, d13cba = runloscarp(timev,zeros(300),zeros(300),exp.(logexpvalsᵣ),co2doublingrate);

        if all(isnan.(tmv))
            llᵣ = NaN
        else

            d13cbmuᵣ = Array{Float64,1}(undef,300);
            nanmean!(d13cbmuᵣ,vec(tmv),d13cba,first(timev),last(timev),length(timev));
            d13cbmuᵣ = fillnans(d13cbmuᵣ,150);
            if isnan(d13cbmuᵣ[300])
                d13cbmuᵣ[first(findall(x->isnan(x),d13cbmuᵣ)):300] .= d13cbmuᵣ[first(findall(x->isnan(x),d13cbmuᵣ))-1]
            end
            llᵣ = normpdf_ll(d13cbvals,d13cberror,d13cbmuᵣ);
        end

        if log(rand()) < (llᵣ-ll)
            ll = llᵣ
            logexpvals .= logexpvalsᵣ
            d13cbmu = d13cbmuᵣ
        end
        # update the latest values
        lldist[i] = ll;
        expdist[:,i] = logexpvals;
        d13cbarray[:,i] = d13cbmu;

    end

        # collate all the final values
    all_ll_dist = MPI.Gather(lldist, 0, comm)
    all_exp_dist = MPI.Gather(expdist,0,comm)
    all_d13cb = MPI.Gather(d13cbarray,0,comm)
    # write them to csv files
    if rank == 0
        writedlm("$loscdir/exp_ll_dist.csv",all_ll_dist,',')
        writedlm("$loscdir/exp_exp_dist.csv",all_exp_dist,',')
        writedlm("$loscdir/all_d13cb.csv",all_d13cb,',')
    end  
end       
MPI.Finalize()