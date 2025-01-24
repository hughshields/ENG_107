############################################
##  file: normal_sampling.R
##   ENG 107 Problem Set 2 Part 1 Solutions
##   - determine the mean and the 95 percentile from a known univariate normal 
##     distribution with a mean of zero and a standard deviation of one with 
##     your estimated uncertainties
######################################################
##  author: Hugh Shields
##  copyright by the author
##  distributed under the GNU general public license
##  https://www.gnu.org/licenses/gpl.html
##  no warranty (see license details at the link above)
##################################################
# version 1: last changes: Jan. 24. 2025 by Hugh Shields
# contact: nikolaus.h.shields.gr@dartmouth.edu
##################################################
# sources:
# - R help files accessed through R-studio for syntax
# - discussions during ENG 107 class time
# - An introduction to R (2010) by Longhow Lam
# - coin-example.R from ENG 107 Canvas
####################################################
# how to run:
# - save the file in a directory
# - go to the directory with this file
# - open R
# - type 'source(normal_sampling.R)'
# - read the printed results in the console and
#   open the new pdf file to analyze the results
############################################

# Clear any existing variables and plots. 
rm(list = ls())
graphics.off()

# Varying Number of Samples #############################

# define a seed for reproducibility
set.seed(0)

# create sample size vector
number_samples = exp(seq(log(1e2), log(1e7), length.out = 6))

# do the simulation for each number of samples and print results
print('Running test for varying sample size')
for (i in seq_along(number_samples)) {
  n = number_samples[i]
  sim=rnorm(n, mean = 0, sd = 1)
  
  # calculate the mean and (5,95) percentiles of each simulation and print
  sim_mean = mean(sim)
  sim_percentiles = quantile(sim, probs = c(0.05, 0.95))
  cat(sprintf("n = %.f, mean = %.5f, (5th,95th) percentiles = (%.5f, %.5f)\n",
              n, sim_mean, sim_percentiles[1], sim_percentiles[2]))
}

# Varying Seeds and Number of Samples ######################

# define a seed for reproducibility
seeds = 1:500

# define the number of samples
sim_means = matrix(NA, nrow = length(seeds), ncol = length(number_samples))

# plot outputs
pdf(file="seed_means.pdf",10,6.17)

# run the simulation for lots of seeds
print('Running test for 500 seeds and varying sample size')
for (i in seq_along(seeds)) {
  set.seed(seeds[i])
  
  # do the simulation for each number of samples for seed i
  #######################
  for (j in seq_along(number_samples)) {
    n = number_samples[j]
    sim=rnorm(n, mean = 0, sd = 1)
    sim_means[i,j]=mean(sim)
  }
  if (i==1) {
    plot(number_samples, sim_means[i,], log = "x", type = "b", ylim = c(-0.3, 0.3),
         xlab = "Number of Samples (log scale)", ylab = "Simulation Means",
         main = "Number of Samples vs Simulation Means for 500 Seeds")
  } else {
    par(new = TRUE)
    plot(number_samples, sim_means[i,], log = "x", type = "b", ylim = c(-0.3, 0.3),
         axes = FALSE, xlab = "", ylab = "")
  }
}
dev.off()

# print the results
for (i in seq_along(number_samples)) {
  # calculate the mean and (5,95) percentiles for each n over 500 seeds
  n = number_samples[i]
  n_mean = mean(sim_means[,i])
  n_percentiles = quantile(sim_means[,i], probs = c(0.05, 0.95))
  cat(sprintf("n = %.f, mean = %.5e, (5th,95th) percentiles = (%.5e, %.5e)\n",
              n, n_mean, n_percentiles[1], n_percentiles[2]))
}
