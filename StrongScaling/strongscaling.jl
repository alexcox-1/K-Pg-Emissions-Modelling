## do the loscar mcmc ac.gr@ 11/12/20
using StatGeochem, MPI, Statistics, DelimitedFiles
MPI.Init()



let
# Get MPI properties
    comm = MPI.COMM_WORLD
    ntasks = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    loscdir = "/dartfs-hpc/scratch/alex/K-Pg-Emissions-Modelling/StrongScaling"
    # Test stdout
    print("Hello from $rank of $ntasks processors!\n")
    
    # a sample 'temperature' distribution, 300 elements long, 0 to 5 degrees
    temp = importdataset("strongscalingtemp.csv",',')
    temp = temp["temp"]
    bsrd13c = importdataset("d13cdatabsr.csv",',')
    d13cvals = bsrd13c["d13cval"][1:300];
    d13cerror = bsrd13c["d13cerror"][1:300];
    # add an error in quadrature given LOSCAR uncertainties (assumed 0.3)
    d13cerror .= sqrt.((d13cerror.^2) .+ 0.3^2)
    # add benthic d13c values
    bsrd13cb = importdataset("d13cdatabenthicbsr.csv",',')
    d13cbvals = bsrd13cb["d13cval"][1:300];
    d13cberror = bsrd13cb["d13cerror"][1:300];
    # add an error in quadrature
    d13cberror = sqrt.((d13cberror.^2) .+ 0.3^2)
    numiter = 15000;
    num_per_exchange = 2;
    ll_dist_array = Array{Float64,2}(undef,numiter,1);
    muarray = Array{Float64,2}(undef,300,numiter);
    all_temps = Array{Float64,3}(undef,300,numiter,ntasks);
    all_co2_dist = Array{Float64,3}(undef,300,numiter,ntasks);
    all_s_dist = Array{Float64,3}(undef,300,numiter,ntasks);
    all_ll_dist = Array{Float64,2}(undef,numiter,ntasks)
    for k = 1:1
	    (rank == 0) && println("Iteration $k")
        ## monte carlo loop
        # perturb one of the co2 vals and one of the svals
        # work with co2 vals and svals in logspace
        # for parallel, do some number of iterations, collect answers, 
        # then continue 
	    temperror = zeros(300) .+ 0.3;
    	co2vals = zeros(300) .+ 0.04;
    	svals = zeros(300) .+ 0.01;
        expvals = ones(300);
    
    	mu = (((0.1 ./ co2vals) .+ 0.04) .* 2) .- (exp.(20 .* svals) .- 1) .+ (0.5 ./ expvals);

        d13cmu = (exp.(10 .* svals) .- 0.5) .+ (1 ./ expvals) .- (co2vals ./ 0.08);

        d13cbmu = ((co2vals ./ 0.15) * 2.5) .- (1 ./ expvals);

   	    ll = normpdf_ll(temp,temperror,mu) + normpdf_ll(d13cvals,d13cerror,d13cmu) + normpdf_ll(d13cbvals,d13cberror,d13cbmu);

    	logco2vals = log.(co2vals);
    	logsvals = log.(svals);
        logexpvals = log.(expvals);
        logco2valsᵣ = copy(logco2vals);
        logsvalsᵣ = copy(logsvals);
        logexpvalsᵣ = copy(logexpvals);
        lldist = Array{Float64,1}(undef,numiter);
        co2dist = Array{Float64,2}(undef,length(logco2vals),numiter);
        sdist = Array{Float64,2}(undef,length(logsvals),numiter);
        expdist = Array{Float64,2}(undef,length(logexpvals),numiter);
        #muarray = Array{Float64,2}(undef,length(mu),numiter);
        # create a record of what other MPI tasks have right now
        ll_dist = Array{Float64,1}(undef,numiter);
        all_log_co2 = Array{Float64,2}(undef, length(logco2vals),ntasks);
        all_log_s = Array{Float64,2}(undef, length(logsvals),ntasks);
        all_log_exp = Array{Float64,2}(undef,length(logexpvals),ntasks)
        all_lls = Array{Float64}(undef,ntasks);
        step_sigma_co2_array = Array{Float64,1}(undef,numiter);
        step_sigma_so2_array = Array{Float64,1}(undef,numiter);
        # which of the loops are we in?
        trialnumber = k
        # do the mcmc
        # set the std of the proposal amplitude distribution
        co2_step_sigma = 0.1;
        so2_step_sigma = 0.1;
        exp_step_sigma = 0.01;
        halfwidthc = 1;
        halfwidths = 1;
        halfwidthexp = 0.5;
        @inbounds for i = 1:5000
           # print("Iteration $i")
            # update current prediction
            copyto!(logco2valsᵣ,logco2vals);
            copyto!(logsvalsᵣ,logsvals);
            copyto!(logexpvalsᵣ,logexpvals);
            # Exchange proposals, sometimes
            if i % num_per_exchange == 0 && i > 1
                # Exchange current proposals across all MPI tasks
                MPI.Allgather!(logco2vals, all_log_co2, comm)
                MPI.Allgather!(logsvals, all_log_s, comm)
                MPI.Allgather!(logexpvals, all_log_exp, comm)
                MPI.Allgather!(lldist[i-1:i-1], all_lls, comm)
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
            # choose which indices to perturb, and perturb it trialnumber times
            # initialize the amplitudes
            randamplitude = 0
            randamplitudes = 0
            randamplitudeexp = 0
            # modify co2 vals
            randhalfwidth = halfwidthc * rand()*length(co2vals)
            randmu = rand()*length(co2vals)
            randamplitude = randn()*co2_step_sigma*2.9
            for j=1:length(co2vals)
                logco2valsᵣ[j] += randamplitude * ((randmu-randhalfwidth)<j<(randmu+randhalfwidth))
            end
        # modify s vals
            randhalfwidths = halfwidths * rand()*length(svals)
            randmus = rand()*length(svals)
            randamplitudes = randn()*so2_step_sigma*2.9
            for j=1:length(svals)
                logsvalsᵣ[j] += randamplitudes * ((randmus-randhalfwidths)<j<(randmus+randhalfwidths))
            end

            randhalfwidthexp = halfwidthexp * rand()*length(svals)
            randmuexp = rand()*length(expvals)
            randamplitudeexp = randn()*exp_step_sigma*2.9
            for j=1:length(svals)
                logexpvalsᵣ[j] += randamplitudeexp * ((randmuexp-randhalfwidthexp)<j<(randmus+randhalfwidthexp))
            end
            # run loscar with the new values
            
            mu = (((0.1 ./ exp.(logco2valsᵣ)) .+ 0.04) .* 2) .- (exp.(20 .* exp.(logsvalsᵣ)) .- 1) .+ (0.5 ./ exp.(logexpvalsᵣ));

            d13cmu = (exp.(10 .* exp.(logsvalsᵣ)) .- 0.5) .+ (1 ./ exp.(logexpvalsᵣ)) .- (exp.(logco2valsᵣ) ./ 0.08);
    
            d13cbmu = ((exp.(logco2valsᵣ) ./ 0.15) * 2.5) .- (1 ./ exp.(logexpvalsᵣ));
            llᵣ = normpdf_ll(temp,temperror,mu) + normpdf_ll(d13cvals,d13cerror,d13cmu) + normpdf_ll(d13cbvals,d13cberror,d13cbmu);
            counter = 0
            # is this allowed?
            if log(rand()) < (llᵣ-ll)
                ll = llᵣ
                logco2vals .= logco2valsᵣ  
                logsvals .= logsvalsᵣ  
                logexpvals .= logexpvalsᵣ
                co2_step_sigma = max(min(abs(randamplitude),1),0.01);
                so2_step_sigma = max(min(abs(randamplitudes),1),0.01);
                exp_step_sigma = max(min(abs(randamplitudeexp),0.5),0.01)
            else
                counter += 1
            end
            if counter >= 20
                halfwidthc *= 0.9;
                halfwidths *= 0.9;
                halfwidthexp *= 0.9;
            end
            # update the latest values
            lldist[i] = ll;
            co2dist[:,i] = logco2vals;
            sdist[:,i] = logsvals;
            expdist[:,i] = logexpvals;
            muarray[:,i] = mu;
            step_sigma_co2_array[i] = co2_step_sigma;
            step_sigma_so2_array[i] = so2_step_sigma;

        end
        @inbounds for i = 5001:15000
            # print("Iteration $i")
             # update current prediction
             copyto!(logco2valsᵣ,logco2vals);
             copyto!(logsvalsᵣ,logsvals);
             copyto!(logexpvalsᵣ,logexpvals);
             # Exchange proposals, sometimes
             # choose which indices to perturb, and perturb it trialnumber times
             # initialize the amplitudes
             randamplitude = 0
             randamplitudes = 0
             randamplitudeexp = 0
             # modify co2 vals
             randamplitude = randn()*0.1
             for j=rand(1:300,2)
                 logco2valsᵣ[j] += randamplitude 
             end
         # modify s vals
             randamplitudes = randn()*0.1
             for j=rand(1:300,2)
                 logsvalsᵣ[j] += randamplitudes 
             end

             randamplitudeexp = randn()*0.01
             for j=rand(1:300,2)
                 logexpvalsᵣ[j] += randamplitudeexp 
             end
             # run loscar with the new values
             
             mu = (((0.1 ./ exp.(logco2valsᵣ)) .+ 0.04) .* 2) .- (exp.(20 .* exp.(logsvalsᵣ)) .- 1) .+ (0.5 ./ exp.(logexpvalsᵣ));

             d13cmu = (exp.(10 .* exp.(logsvalsᵣ)) .- 0.5) .+ (1 ./ exp.(logexpvalsᵣ)) .- (exp.(logco2valsᵣ) ./ 0.08);
    
             d13cbmu = ((exp.(logco2valsᵣ) ./ 0.15) * 2.5) .- (1 ./ exp.(logexpvalsᵣ));
             llᵣ = normpdf_ll(temp,temperror,mu) + normpdf_ll(d13cvals,d13cerror,d13cmu) + normpdf_ll(d13cbvals,d13cberror,d13cbmu);
             counter = 0
             # is this allowed?
             if log(rand()) < (llᵣ-ll)
                 ll = llᵣ
                 logco2vals .= logco2valsᵣ  
                 logsvals .= logsvalsᵣ  
                 logexpvals .= logexpvalsᵣ
             end
             # update the latest values
             lldist[i] = ll;
             co2dist[:,i] = logco2vals;
             sdist[:,i] = logsvals;
             expdist[:,i] = logexpvals;
             muarray[:,i] = mu;
             step_sigma_co2_array[i] = co2_step_sigma;
             step_sigma_so2_array[i] = so2_step_sigma;
         end
        # collate all the final values

	    (rank == 0) && println("Printing lldist from Rank $rank")
        all_ll_dist = MPI.Gather(lldist, 0, comm)
        (rank == 0) && (all_ll_dist = reshape(all_ll_dist,numiter,ntasks))
	    # (rank == 0) && (all_ll_dist = mean(all_ll_dist,dims=2))
        # (rank == 0) && (ll_dist_array[:,k] = all_ll_dist)
	    (rank == 0) && println("Printing all_ll_dist from Rank $rank")
	    #(rank == 0) && print(all_ll_dist)
        #all_co2_dist = MPI.Gather(co2dist, 0, comm)
        #all_s_dist = MPI.Gather(sdist, 0, comm)
        #all_temps = MPI.Gather(muarray, 0, comm)
        #all_d13c = MPI.Gather(d13carray,0,comm)
        #all_co2_step = MPI.Gather(step_sigma_co2_array,0,comm)
        #all_so2_step = MPI.Gather(step_sigma_so2_array,0,comm)
	
    end
    # write them to csv files
    if rank == 0
	
        print("About to write files!")
        writedlm("$loscdir/ll_dist_2_$ntasks.csv",all_ll_dist,',')
        #writedlm("$loscdir/all_co2_dist_$ntasks.csv",all_co2_dist,',')
        #writedlm("$loscdir/all_s_dist_$ntasks.csv",all_s_dist,',')
        #writedlm("$loscdir/all_temps_$ntasks.csv",all_temps,',')
        #writedlm("$loscdir/all_co2_step.csv",all_co2_step,',')
        #writedlm("$loscdir/all_so2_step.csv",all_so2_step,',')
        #writedlm("$loscdir/all_d13c.csv",all_d13c,',')
    end  
end       
MPI.Finalize()
