############################################
##  file: ps3_shields.R
##  ENG 107 Problem Set 3 Solutions
######################################################
##  author: Hugh Shields
##  copyright by the author
##  distributed under the GNU general public license
##  https://www.gnu.org/licenses/gpl.html
##  no warranty (see license details at the link above)
##################################################
# version 1: last changes: Feb. 7, 2025 by Hugh Shields
# contact: nikolaus.h.shields.gr@dartmouth.edu
##################################################
# citations:
# - R help files accessed through R-studio for syntax
# - discussions during ENG 107 class time
# - An introduction to R (2010) by Longhow Lam
# - example codes from ENG 107 Canvas
# - truncnorm package documentation 
#   (https://cran.r-project.org/web/packages/truncnorm/truncnorm.pdf)
####################################################
# how to run:
# - save the file in a directory
# - go to the directory with this file
# - open R
# - type 'install.packages("truncnorm")'
# - type 'source(ps3_shields.R)'
# - read the printed results in the console and
#   open the new pdf file to analyze the results
############################################

# Clear any existing variables and plots. 
rm(list = ls())
graphics.off()

# define a seed for reproducibility
set.seed(0)

# 2. plot PDF of useable fuel
pdf(file="PDF_of_usable_fuel.pdf", 6,4)
fuel = seq(-36, 104, length=10000)
fuel.mean = 34
fuel.sd = 20
prob = dnorm(fuel, mean=fuel.mean, sd=fuel.sd)
plot(fuel, prob, type="l", lwd=1, xlab="Fuel (L)", ylab = "Probability Density",
     main = "PDF of Avaliable Fuel", xlim = c(-32, 100))
abline(v=34,col="red",lty=1,lwd=4)
dev.off()

# print the probability of negative fuel
neg_prob = pnorm(0, mean = fuel.mean, sd = fuel.sd)
print(sprintf("2. The probability of negative fuel in the tank is %.2f%%", neg_prob*100))

# 4. Grid method to determine posterior
fuel.min = 0
fuel.max = 182
grid = seq(fuel.min-50, fuel.max+50, by = 0.01)
# create a uniformly distributed prior between [0,128] and normalize to sum to 1
prior = ifelse(grid < 0 | grid > 128, 0, 1 / ((128/0.01)+1))
# specify likelihood from above
likelihood = dnorm(grid, mean=fuel.mean, sd=fuel.sd)
# generate normalized posterior
posterior = (prior * likelihood)/sum(prior * likelihood)

# plot grid posterior over original PDF
pdf(file="grid_posterior.pdf", 6,4)
plot(fuel, prob, type="l", lwd=1, xlab="Fuel (L)", ylab = "Probability Density",
     main = "Distribution of Avaliable Fuel", xlim = c(-32, 100))
par(new = TRUE)
plot(grid, posterior, type="l", col="blue", lwd=1, axes = FALSE, xlab = "", ylab = "",
     xlim = c(-32, 100))
legend("topright", legend = c("Original PDF", "Grid Posterior"), col = c("black", "blue"), lwd = 2)
dev.off()

# 5. Bayesian Monte Carlo to determine posterior
# sample N values from the prior
N = 100000
prior_sample = runif(N, min = fuel.min, max = fuel.max)
# compute the likelihood values from the prior sampling
likelihood_values = dnorm(prior_sample, mean = fuel.mean, sd = fuel.sd)
weights = likelihood_values/sum(likelihood_values)
# resample from the prior sample using the weights (with replacement)
posterior_samples = sample(prior_sample, size = N, replace = TRUE, prob = weights)

# plot histogram (probability=TRUE shows density) of BMC results and overplot original PDF
pdf(file="BMC_posterior.pdf", 6,4)
hist(posterior_samples, breaks=50, probability=TRUE, xlab="Fuel (L)", 
     ylab = "Probability Density", main = "Distribution of Avaliable Fuel",
     col="lightblue", border="blue", xlim=c(-32, 100), ylim=c(0, 0.022))
par(new = TRUE)
plot(fuel, prob, type="l", lwd=1, xlab="", ylab = "", main = "",
     xlim = c(-32, 100), ylim=c(0, 0.022))
legend("topright", legend = c("Original PDF", "BMC Posterior"), col = c("black", "lightblue"), lwd = 2)
dev.off()

# 7. Computing distribution of remaining flight time
N = 100000
# install and load truncnorm package to generate samples from the fuel distribution
library(truncnorm)
# sample the available fuel distribution (truncated normal distribution)
fuel_sample = rtruncnorm(N, mean = fuel.mean, sd = fuel.sd, a = fuel.min, b = fuel.max)
# sample the rate of fuel consumption distribution
rate.mean = 18/60
rate.sd = 2/60
rate_sample = rnorm(N, mean=rate.mean, sd=rate.sd)

# compute available flight time
time_sample = fuel_sample/rate_sample

# plot histogram of time when fuel runs out
pdf(file="fuel_depletion_dist.pdf", 6,4)
hist(time_sample, breaks=100, probability=TRUE, xlab="Fuel Depletion Time (min)", 
     ylab = "Probability Density", main = "Distribution Fuel Depletion Times",
     col="grey", border="black")
# plot lines at 100 and 130 minutes
abline(v=130,col="green",lty=1,lwd=4)
abline(v=100,col="red",lty=1,lwd=4)
dev.off()

# 7.a compute the probability of arriving at the airport 100 min away with 30 
# min of fuel reserve
successes = time_sample[time_sample>=130]
success_rate = length(successes)/N
print(sprintf("7a. The probability of arriving with 30 minutes of fuel reserve is %.2f%%", success_rate*100))

# 7.a compute the probability of running out of fuel befor arriving at the 
# airport 100 min away
failures = time_sample[time_sample<100]
failure_rate = length(failures)/N
print(sprintf("7b. The probability of running out of fuel is %.2f%%", failure_rate*100))

