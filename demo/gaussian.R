# LEAST SQUARES EXAMPLE

cat("Running least squares regression example.\n")
# create some data

N <- 1000
X1 <- runif(N)
X2 <- 2*runif(N)
X3 <- factor(sample(letters[1:4],N,replace=T))
X4 <- ordered(sample(letters[1:6],N,replace=T))
X5 <- factor(sample(letters[1:3],N,replace=T))
X6 <- 3*runif(N)
mu <- c(-1,0,1,2)[as.numeric(X3)]

SNR <- 10 # signal-to-noise ratio
Y <- X1**1.5 + 2 * (X2**.5) + mu
sigma <- sqrt(var(Y)/SNR)
Y <- Y + rnorm(N,0,sigma)

X1[sample(1:N,size=100)] <- NA
X3[sample(1:N,size=300)] <- NA

# random weights if you want to experiment with them
# w <- rexp(N)
# w <- N*w/sum(w)
w <- rep(1,N)

data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)

# fit initial model
gbm1 <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
            data=data,                   # dataset
            weights=w,
            var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
            distribution="gaussian",     # bernoulli, adaboost, gaussian, poisson, and coxph available
            n.trees=100,                 # number of trees
            shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
            interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
            bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
            train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
            n.minobsinnode = 10,         # minimum number of obs needed in each node
            keep.data=T)

# check performance and do another 100 iterations
best.iter <- gbm.perf(gbm1,method="OOB")
gbm2 <- gbm.more(gbm1,100)

best.iter <- gbm.perf(gbm2,method="OOB")
while(gbm2$n.trees - best.iter < 10)
{    
   # do 100 more iterations
   gbm2 <- gbm.more(gbm2,100)
   best.iter <- gbm.perf(gbm2,plot.it=F,method="OOB")
}

# plot the performance
best.iter <- gbm.perf(gbm2,method="OOB")  # returns out-of-bag estimated best number of trees
best.iter <- gbm.perf(gbm2,method="test") # returns test set estimate of best number of trees

# plot variable influence
summary(gbm2,n.trees=1)         # based on the first tree
summary(gbm2,n.trees=best.iter) # based on the estimated best number of trees

# print the first and last trees
print(pretty.gbm.tree(gbm2,1))
print(pretty.gbm.tree(gbm2,gbm1$n.trees))

print(gbm2$c.splits[1:3])

# make some new data
N <- 1000
X1 <- runif(N)
X2 <- 2*runif(N)
X3 <- ordered(sample(letters[1:4],N,replace=T))
X4 <- factor(sample(letters[1:6],N,replace=T))
X5 <- factor(sample(letters[1:3],N,replace=T))
X6 <- 3*runif(N)
mu <- c(-1,0,1,2)[as.numeric(X3)]

Y <- X1**1.5 + 2 * (X2**.5) + mu
Y <- Y + rnorm(N,0,sigma)

data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
print(data2[1:10,])

# predict on the new data using "best" number of trees
f.predict <- predict.gbm(gbm2,data2,best.iter) # f.predict will be on the canonical scale (logit,log,etc.)

print(f.predict[1:10])
# least squares error
print(sum((data2$Y-f.predict)^2))

# create marginal plots
# plot variable X1,X2,X3 after "best" iterations
par(mfrow=c(1,3))
plot.gbm(gbm2,1,best.iter)
plot.gbm(gbm2,2,best.iter)
plot.gbm(gbm2,3,best.iter)
par(mfrow=c(1,1))
plot.gbm(gbm2,1:2,best.iter) # contour plot of variables 1 and 2 after "best" number iterations
plot.gbm(gbm2,2:3,best.iter) # lattice plot of variables 2 and 3 after "best" number iterations
plot.gbm(gbm2,3:4,best.iter) # lattice plot of variables 2 and 3 after "best" number iterations

plot.gbm(gbm2,c(1,2,6),best.iter,cont=20) # 3-way plots
plot.gbm(gbm2,1:3,best.iter)
plot.gbm(gbm2,2:4,best.iter)
plot.gbm(gbm2,3:5,best.iter)
