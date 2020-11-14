using MPI
MPI.Init()
# Find out who and where we are
size = MPI.Comm_size(MPI.COMM_WORLD)
rank = MPI.Comm_rank(MPI.COMM_WORLD)
print("Hello from $rank of $size processors!\n")
# Finalize
MPI.Finalize()