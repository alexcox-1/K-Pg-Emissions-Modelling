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
    numiter = 100000;
    num_per_exchange = 1;
    ll_dist_array = Array{Float64,2}(undef,numiter,5);
    muarray = Array{Float64,2}(undef,300,numiter);
    for k = 1:5
	    (rank == 0) && println("Iteration $k")
        ## monte carlo loop
        # perturb one of the co2 vals and one of the svals
        # work with co2 vals and svals in logspace
        # for parallel, do some number of iterations, collect answers, 
        # then continue 
	    temperror = zeros(300) .+ 0.3;
    	co2vals = zeros(300) .+ 0.04;
    	svals = zeros(300) .+ 0.01;
    
    	mu = (((co2vals ./ 0.1) .+ 0.04) .* 2) .- (exp.(20 .* svals) .- 1)

   	    ll = normpdf_ll(temp,temperror,mu) 
    	logco2vals = log.(co2vals);
    	logsvals = log.(svals);
        logco2valsᵣ = copy(logco2vals);
        logsvalsᵣ = copy(logsvals);
        lldist = Array{Float64,1}(undef,numiter);
        co2dist = Array{Float64,2}(undef,length(logco2vals),numiter);
        sdist = Array{Float64,2}(undef,length(logsvals),numiter);
        #muarray = Array{Float64,2}(undef,length(mu),numiter);
        # create a record of what other MPI tasks have right now
        all_log_co2 = Array{Float64}(undef, length(logco2vals), ntasks);
        all_log_s = Array{Float64}(undef, length(logsvals), ntasks);
        all_lls = Array{Float64}(undef,ntasks);
        step_sigma_co2_array = Array{Float64,1}(undef,numiter);
        step_sigma_so2_array = Array{Float64,1}(undef,numiter);
        # which of the loops are we in?
        trialnumber = k
        # do the mcmc
        # set the std of the proposal amplitude distribution
        co2_step_sigma = 0.1;
        so2_step_sigma = 0.1;
        @inbounds for i = 1:numiter
           # print("Iteration $i")
            # update current prediction
            copyto!(logco2valsᵣ,logco2vals);
            copyto!(logsvals,logsvals);
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
            randhalfwidth = rand()*length(co2vals)/100

            randmu = rand()*length(co2vals)
            randamplitude = randn()*co2_step_sigma*trialnumber

            for j=1:length(co2vals)
                logco2valsᵣ[j] += randamplitude * ((randmu-randhalfwidth)<j<(randmu+randhalfwidth))

            end
            randhalfwidths = rand()*length(svals)/100

            randmus = rand()*length(svals)

            randamplitudes = randn()*so2_step_sigma*trialnumber

            for j=1:length(svals)
                logsvalsᵣ[j] += randamplitudes * ((randmus-randhalfwidths)<j<(randmus+randhalfwidths))

            end
            # run loscar with the new values
            
            mu = (((exp.(logco2valsᵣ) ./ 0.1) .+ 0.04) .* 2) .- (exp.(20 .* exp.(logsvalsᵣ)) .- 1)
            llᵣ = normpdf_ll(temp,temperror,mu)

            # is this allowed?
            if log(rand()) < (llᵣ-ll)
                ll = llᵣ
                logco2vals .= logco2valsᵣ  
                logsvals .= logsvalsᵣ  
                co2_step_sigma = min(abs(randamplitude),1);
                so2_step_sigma = min(abs(randamplitudes),1)
            end
            # update the latest values
            lldist[i] = ll;
            co2dist[:,i] = logco2vals;
            sdist[:,i] = logsvals;
            muarray[:,i] = mu;
            step_sigma_co2_array[i] = co2_step_sigma;
            step_sigma_so2_array[i] = so2_step_sigma;
        end

        # collate all the final values

	    (rank == 0) && println("Printing lldist from Rank $rank")
        all_ll_dist = MPI.Gather(lldist, 0, comm)
        (rank == 0) && (all_ll_dist = reshape(all_ll_dist,numiter,ntasks))
	    (rank == 0) && (all_ll_dist = mean(all_ll_dist,dims=2))
        (rank == 0) && (ll_dist_array[:,k] = all_ll_dist)
	    (rank == 0) && println("Printing all_ll_dist from Rank $rank")
	    #(rank == 0) && print(all_ll_dist)
        #all_co2_dist = MPI.Gather(co2dist, 0, comm)
        #all_s_dist = MPI.Gather(sdist, 0, comm)
        all_temps = MPI.Gather(muarray, 0, comm)
        #all_d13c = MPI.Gather(d13carray,0,comm)
        #all_co2_step = MPI.Gather(step_sigma_co2_array,0,comm)
        #all_so2_step = MPI.Gather(step_sigma_so2_array,0,comm)
	
    end
    # write them to csv files
    if rank == 0
	
        print("About to write files!")
        writedlm("$loscdir/ll_dist_1_$ntasks.csv",ll_dist_array,',')
        #writedlm("$loscdir/all_co2_dist.csv",all_co2_dist,',')
        #writedlm("$loscdir/all_s_dist.csv",all_s_dist,',')
        writedlm("$loscdir/all_temps_$ntasks.csv",muarray[:,50000:100000],',')
        #writedlm("$loscdir/all_co2_step.csv",all_co2_step,',')
        #writedlm("$loscdir/all_so2_step.csv",all_so2_step,',')
        #writedlm("$loscdir/all_d13c.csv",all_d13c,',')
    end  
end       
MPI.Finalize()
