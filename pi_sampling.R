############################################
##  file: pi_sampling.R
##   ENG 107 Problem Set 2 Part 2 Solutions
##   - determine the value of pi with your estimated uncertainties
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
# - type 'source(pi_sampling.R)'
# - read the printed results in the console and
#   open the new pdf file to analyze the results
############################################

# Clear any existing variables and plots. 
rm(list = ls())
graphics.off()

# Varying Seeds and Number of Samples ######################

# define a seed for reproducibility
seeds = 1:100

# create sample size vector
number_samples = exp(seq(log(1e2), log(1e7), length.out = 6))

# create a matrix to hold the calculated pi values
pi_calc = matrix(NA, nrow = length(seeds), ncol = length(number_samples))

# plot outputs
pdf(file="pi_means.pdf",10,6.17)

# run the simulation for lots of seeds
print('Running test for 100 seeds and varying sample size')
for (i in seq_along(seeds)) {
  set.seed(seeds[i])
  
  # calculate pi for seed i for each number of samples
  #######################
  for (j in seq_along(number_samples)) {
    n = number_samples[j]
    x = runif(n, min = -0.5, max = 0.5)
    y = runif(n, min = -0.5, max = 0.5)
    r = sqrt(x^2+y^2)
    area_rat = length(r[r<=0.5])/length(r)
    pi = area_rat/(0.5^2)
    pi_calc[i,j]=pi
  }
  if (i==1) {
    plot(number_samples, pi_calc[i,], log = "x", type = "b", 
         ylim = c(pi-0.35, pi+0.35),
         xlab = "Number of Samples (log scale)", ylab = "Calculated Pi Value",
         main = "Number of Samples vs Calculated Pi Value for 100 Seeds")
  } else {
    par(new = TRUE)
    plot(number_samples, pi_calc[i,], log = "x", type = "b", 
         ylim = c(pi-0.35, pi+0.35),
         axes = FALSE, xlab = "", ylab = "")
  }
}
dev.off()

# print the results
for (i in seq_along(number_samples)) {
  # calculate the mean and (5,95) percentiles for each n over 100 seeds
  n = number_samples[i]
  n_mean = mean(pi_calc[,i])
  n_percentiles = quantile(pi_calc[,i], probs = c(0.05, 0.95))
  cat(sprintf("n = %.f, mean pi value = %.5f, 
              (5th,95th) percentiles = (%.5f, %.5f)\n",
              n, n_mean, n_percentiles[1], n_percentiles[2]))
}
