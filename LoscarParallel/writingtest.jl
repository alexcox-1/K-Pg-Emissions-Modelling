## do the loscar mcmc ac.gr@ 11/12/20
using DelimitedFiles, StatGeochem, MPI, Statistics
MPI.Init()



let
# Get MPI properties
    comm = MPI.COMM_WORLD
    ntasks = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    loscdir = "/dartfs-hpc/scratch/alex/K-Pg-Emissions-Modelling/LoscarParallel"
    # Test stdout
    print("Hello from $rank of $ntasks processors!\n")
    testarray = 1:1:100; 
    # a sample 'temperature' distribution, 300 elements long, 0 to 5 degrees
    # write them to csv files
    lldist = randn(100)
    testarray = MPI.Gather(lldist,0,comm)
    println("Printing Gathered lldist")
    print(testarray)
    if rank == 0
	print("About to write files!")
        writedlm("$loscdir/test_array.csv",testarray,',')
    end  
end       
MPI.Finalize()
