
################################### MULTIVARIATE PROJECT ############################################

library(biotools)
library(Amelia)
library(MVN)
library(mvoutlier)
library(cluster)
library(HSAUR)
library(rgl)


econ <- read.csv("C:/Users/Alberto/Google Drive/St Andrews/Multivariate/Final Practical/replace/replace/newwdi.csv",
                 sep = ";")
econ$Time <- NULL
sum(is.na(econ))/(dim(econ)[1]*dim(econ)[2])
country <- econ$Country.Name
econ$Country.Name <- NULL
# We still have some NAs: would be a waste to not include them. Moreover, the NAs are just 4% of
# the values, hence no big influence.
# Estimate those values through Imputation

# Try to log(base 10) the bigger values (or other transformation)
econ[,3] <- log(econ[,3], base= 10)


imp.econ <- amelia(econ, parallel = "snow", ncpus = 6, m=1)
summary(imp.econ)

# Some EDA and Assumption checking
# Linearity
pairs(imp.econ$imputations$imp1)
plot(imp.econ$imputations$imp1[,6],imp.econ$imputations$imp1[,13],
     xlab = "Labor force participation rate female", ylab = "Population ages 0-14", col= "red",
     bg= "red", pch= 20)
# Equality of var-cov matrices

# Multivariate Normality
centroid <- apply(imp.econ$imputations$imp1, 2, mean)
mahal <- sort(mahalanobis(imp.econ$imputations$imp1, centroid, cov(imp.econ$imputations$imp1)))
quant <- (1:length(mahal))/(length(mahal)+1)
chi <- qchisq(quant, 143) # ASK TO FIND DEGREES OF FREEDOM
plot(chi, mahal)
abline(20, 1, col= "red")

hzTest(imp.econ$imputations$imp1, qqplot= T) # Not MVN
# Look for symmetry
apply(imp.econ$imputations$imp1, 2, hist) # Most variables are highly skewed


# Outliers: to check them after the PCA is done
chisq.plot(imp.econ$imputations$imp1) # A lot of outliers based on mahalanobis distance


# Execute PCA


test <- princomp(imp.econ$imputations$imp1, cor= F) # Results vary significantly depending on cov/cor

no.gdp <- imp.econ$imputations$imp1[, -c(2,3)]
test1 <- princomp(no.gdp, cor= F)

# Screeplot
screeplot(test, type= "l") # choose first 3 in both
# Kaiser
eig.mean<- mean(test$sdev^2) 
which(test$sdev^2>=eig.mean) # Choose 6 in cor= F, 9 in cor= T


# Try a simple plot
par(mfrow=c(1,1))

# 1 vs 2
par(pty="s")
plot(test$scores[,1],test$scores[,2],
     ylim=range(test$scores[,1]),
     xlab="PC1",ylab="PC2",type="n",lwd=2)
text(test$scores[,1],test$scores[,2],
     labels=abbreviate(country),cex=0.7,lwd=2)
# 1 vs 3
plot(test$scores[,1],test$scores[,3],
     ylim=range(test$scores[,1]),
     xlab="PC1",ylab="PC3",type="n",lwd=2)
text(test$scores[,1],test$scores[,3],
     labels=abbreviate(country),cex=0.7,lwd=2)

# 2 vs 3
plot(test$scores[,2],test$scores[,3],
     ylim=range(test$scores[,2]),
     xlab="PC2",ylab="PC3",type="n",lwd=2)
text(test$scores[,2],test$scores[,3],
     labels=abbreviate(country),cex=0.7,lwd=2)


# Try 3D plot
plot3d(test$scores[,1],test$scores[,2], test$scores[,3], type = "n")
text3d(test$scores[,1],test$scores[,2], test$scores[,3], abbreviate(country))


# Try Bubbleplot with GDP as external
plot(test1$scores[,1],test1$scores[,2], xlab="PC1",ylab="PC2", cex= econ[,2], main="EconData, 1st and 2nd PC, 
     by GDP growth", type= "p") # GDP MUST BE ALREADY IN LOG UNITS!!


plot(test1$scores[,1],test1$scores[,2],
     ylim=range(test1$scores[,1]),
     xlab="PC1",ylab="PC2",type="n",lwd=2)
text(test1$scores[,1],test1$scores[,2],
     labels=abbreviate(country),cex=0.7,lwd=2)

# ========================================================================================================= #
# CLUSTERING


km <- kmeans(imp.econ$imputations$imp1, 2)
dissE <- daisy(imp.econ$imputations$imp1) 
dE2 <- dissE^2

# Plot 1: gives back how "correctly" each observation is clustered
sk2   <- silhouette(km$cl, dE2)
plot(sk2)


# Plot 2: scatterplot with cluster colors. Try different pairs of covariate to find a reasonable interpretation
plot(econ$SG.GEN.PARL.ZS, econ$SE.SEC.ENRR.FE, col=c("red","blue","black")[km$cluster]) # Interesting! No relationship!
plot(econ$SP.POP.DPND, econ$NY.GDP.MKTP.PP.CD, col=c("red","blue","black")[km$cluster])
abline(a= 10, b=0, lty= 3); abline(v=62, lty= 3)

# Trying clustering on the PCA results
km <- kmeans(test$scores[,1:3], 2)
pairs(test$scores[,1:3], col=c("red","blue","black")[km$cluster])

# ANOVA for each variable
new.dat <- cbind(imp.econ$imputations$imp1, cluster= km$cluster)
f <- NULL
for(i in 1:33){
  d <- new.dat[,i]
  res <- summary(aov(d~cluster, data = new.dat))
  f[i] <- res[[1]]$F[1]
}

round(f)
names(f) <- colnames(new.dat)
sort(f, decreasing = T)


