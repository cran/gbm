# LOGISTIC REGRESSION EXAMPLE

cat("Running logistic regression example.\n")
# create some data
N <- 1000
X1 <- runif(N)
X2 <- runif(N)
X3 <- factor(sample(letters[1:4],N,replace=T))
mu <- c(-1,0,1,2)[as.numeric(X3)]

p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
Y <- rbinom(N,1,p)

# random weights if you want to experiment with them
w <- rexp(N)
w <- N*w/sum(w)

data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)

# fit initial model
gbm1 <- gbm(Y~X1+X2+X3,                # formula
            data=data,                 # dataset
            weights=w,
            var.monotone=c(0,0,0),     # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
            distribution="bernoulli",
            n.trees=1000,              # number of trees
            shrinkage=0.001,           # shrinkage or learning rate, 0.001 to 0.1 usually work
            interaction.depth=3,       # 1: additive model, 2: two-way interactions, etc
            bag.fraction = 0.5,        # subsampling fraction, 0.5 is probably best
            train.fraction = 0.5,      # fraction of data for training, first train.fraction*N used for training
            n.minobsinnode = 10)       # minimum total weight needed in each node

# check performance and do another 100 iterations
best.iter <- gbm.perf(gbm1,best.iter.calc="OOB")
gbm2 <- gbm.more(gbm1,100,data=data,weights=w)

best.iter <- gbm.perf(gbm2,best.iter.calc="OOB")
while(gbm2$n.trees - best.iter < 10)
{    
   # do 100 more iterations
   gbm2 <- gbm.more(gbm2,100,data=data,weights=w)
   best.iter <- gbm.perf(gbm2,plot.it=F,best.iter.calc="OOB")
}

# plot the performance
best.iter <- gbm.perf(gbm2,best.iter.calc="OOB")  # returns out-of-bag estimated best number of trees
best.iter <- gbm.perf(gbm2,best.iter.calc="test") # returns test set estimate of best number of trees

# plot variable influence
summary(gbm2,n.trees=1)         # based on the first tree
summary(gbm2,n.trees=best.iter) # based on the estimated best number of trees

# create marginal plots
# plot variable X1,X2,X3 after "best" iterations
par(mfrow=c(1,3))
plot.gbm(gbm2,1,best.iter)
plot.gbm(gbm2,2,best.iter)
plot.gbm(gbm2,3,best.iter)
par(mfrow=c(1,1))
plot.gbm(gbm2,1:2,best.iter) # contour plot of variables 1 and 2 after "best" number iterations
plot.gbm(gbm2,2:3,best.iter) # lattice plot of variables 2 and 3 after "best" number iterations

# 3-way plot
plot.gbm(gbm2,1:3,best.iter)

# print the first and last trees
print(pretty.gbm.tree(gbm2,1))
print(pretty.gbm.tree(gbm2,gbm2$n.trees))

# make some new data
N <- 1000
X1 <- runif(N)
X2 <- runif(N)
X3 <- factor(sample(letters[1:4],N,replace=T))
mu <- c(-1,0,1,2)[as.numeric(X3)]

p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
Y <- rbinom(N,1,p)
data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)

# predict on the new data using "best" number of trees
f.predict <- predict.gbm(gbm2,data2,best.iter) # f.predict will be on the canonical scale (logit,log,etc.)
p.pred <- 1/(1+exp(-f.predict)) # transform to probability scale for logistic regression

# calibration plot for logistic regression - well calibrated means a 45 degree line
par(mfrow=c(1,1))
calibrate.plot(Y,p.pred)

# logistic error
sum(data2$Y*f.predict - log(1+exp(f.predict)))
