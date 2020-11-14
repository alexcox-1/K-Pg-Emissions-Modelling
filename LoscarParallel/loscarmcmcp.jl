## do the loscar mcmc ac.gr@ 11/12/20
using StatGeochem, MPI
MPI.Init()

include("runloscarp.jl")
include("fillnans.jl")


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
bsrtemps = importdataset("/LoscarParallel/tempdatabsr.csv",',');
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
co2vals = zeros(300) .+ 0.01;
svals = zeros(300) .+ 0.01;
logco2vals = log.(co2vals);
logsvals = log.(svals);

co2doublingrate = 3;

# do a loscar run!

tmv,pco2,loscartemp = runloscarp(timev,co2vals,svals,co2doublingrate);

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
mu = fillnans(mu,5);

ll = normpdf_ll(temp,temperror,mu);

numiter = 10000;
num_per_exchange = 100;
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
        map!(x -> x*rand(), all_lls, all_lls)
        chosen = argmax(all_lls)
        logco2valsᵣ .= view(all_log_co2, :, chosen)
        logsvalsᵣ .= view(all_log_s, :, chosen)
    end
    # choose which indices to perturb
    randsindex = rand(1:300,30);
    randcindex = rand(1:300,30);
    # perturb these indices
    logco2valsᵣ[randcindex] .+= (randn.(30)./20);
    logsvalsᵣ[randsindex] .+= (randn.(30)./20);
    # run loscar with the new values
    tmv,pco2,loscartemp = runloscarp(timev,exp.(logco2valsᵣ),exp.(logsvalsᵣ),co2doublingrate);
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

    # is this allowed?
    if log(rand()) < (llᵣ-ll)
        ll = llᵣ
        logco2vals .= logco2valsᵣ  
        logsvals .= logsvalsᵣ  
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
    writedlm(all_ll_dist,"$loscdir/all_ll_dist.csv",',')
    writedlm(all_co2_dist,"$loscdir/all_co2_dist.csv",',')
    writedlm(all_s_dist,"$loscdir/all_s_dist.csv",',')
    writedlm(all_temps,"$loscdir/all_temps.csv",',')
end  
    
MPI.Finalize()