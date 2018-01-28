# Bayes classifier Script
# Input - The proteins and their PAS, mutational, methylation features
# Somnath Tagore

# Install packages
install.packages("cluster")
install.packages("fpc")
install.packages("mclust")
install.packages("ggplot2")
install.packages("maptools")
#install.packages("fviz_cluster")
install.packages("vegan")
#install.packages("ordipointlabel")

# Time calculation
normal_start_time <- proc.time()
# Read in the data and sort it
#dat <- as.vector(read.table("data_1.txt",header=F))
#dat <- read.csv("ST-PAS-og.csv", sep=",")

# Input file
dat <- read.csv("NBC_test_data.csv", sep=",")
namebank <- dat[,3]
namebank
#dat <- sort(dat[,11])
dat <- dat[,11]
dat

# Set up number of simulations and number of parameters to save
reps <- 48
save <- 48
dist_vals <- rep(NA,reps)

# Draw values from the prior distributions
tmp_mean <- runif(reps, -20, 20)
tmp_variance <- runif(reps, 0, 50)

# Loop through the prior values, simulate data and compare to observed data
for(i in 1:reps){
  tmp_dat <- rnorm(length(dat), tmp_mean[i], sqrt(tmp_variance[i]))
  tmp_dat <- sort(tmp_dat)
  dist_vals[i] <- dist(rbind(dat,tmp_dat))
}

# Sort distance values (summary stat)
dist_indexes <- sort(dist_vals, index.return=T)

# Get the top values based on the number of parameter values you want to save
save_indexes <- dist_indexes$ix[1:save]

# Get the corresponding parameter values with the lowest distances from the observed data
saved_means <- tmp_mean[save_indexes]
saved_variances <- tmp_variance[save_indexes]
normal_total_time <- proc.time() - normal_start_time

# Store in a data frame and print as a knitr table
res <- data.frame(Parameter=c("Normal Mean","Normal Variance"),
                  Mean=c(mean(saved_means), mean(saved_variances)),
                  SD=c(sd(saved_means),sd(saved_variances)))
knitr::kable(res, digits=3)
#par(mfrow=c(1,2))
par(mfrow=c(1,1))
normal_start_time <- proc.time()

# K-Means Clustering with 5 clusters
fit <- kmeans(saved_means, 3)
fit

# Centroid Plot against 1st 2 discriminant functions
###fit$cluster
row.names(dat)
namebank
library(fpc)
library(maptools)
require("vegan")

plot(saved_means, col =(fit$cluster +1),  main="blca OG-TS Clusters",pch=20,cex=0.3)
identify(x=saved_means, y = NULL, labels = namebank, pos = FALSE,
        n = length(saved_means), plot = TRUE, atpen = FALSE, offset = 0.5,
         tolerance = 0.25,cex=0.7)
