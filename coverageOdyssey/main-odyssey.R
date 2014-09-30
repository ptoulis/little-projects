# read arguments.
args <- as.numeric(commandArgs(trailingOnly = TRUE))

if(length(args) != 5) {
 stop("Not correct no. of args")
}

true.mu = args[1]
true.sigma = args[2]
N = args[3] # no.samples
nreps = args[4] # repetitions/node
job.id = args[5] # job id

stopifnot(true.sigma > 0, N > 0, nreps > 0)

 # main code.
out.coverage = c()

for(i in 1:nreps) {
  # sample dataset
  y = rnorm(N, mean=true.mu, sd=true.sigma)
  # Construct the CI
  se = 1.96 * true.sigma / sqrt(N) 
  CI = c(mean(y) - se, mean(y) + se)
  # add 0,1 based on coverage
  is.covered = (true.mu >= CI[1] & true.mu <= CI[2])
  out.coverage = c(out.coverage, as.numeric(is.covered)) 
}

# out.coverage should have length nreps, (0, 1, ..)
# Store the vector.
print(job.id)
save(out.coverage, file=sprintf("odyssey/coverage_%d.rda", job.id))




