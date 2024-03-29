#' Test the \code{gbm} package.
#' 
#' Run tests on \code{gbm} functions to perform logical checks and
#' reproducibility.
#' 
#' The function uses functionality in the \code{RUnit} package. A fairly small
#' validation suite is executed that checks to see that relative influence
#' identifies sensible variables from simulated data, and that predictions from
#' GBMs with Gaussian, Cox or binomial distributions are sensible,
#' 
#' @aliases validate.gbm test.gbm test.relative.influence
#' @return An object of class \code{RUnitTestData}. See the help for
#' \code{RUnit} for details.
#' @note The test suite is not comprehensive.
#' @author Harry Southworth
#' @seealso \code{\link{gbm}}
#' @keywords models
#' @examples
#' 
#' # Uncomment the following lines to run - commented out to make CRAN happy
#' #library(RUnit)
#' #val <- validate.texmex()
#' #printHTMLProtocol(val, "texmexReport.html")

#' @export
test.gbm <- function(){
    # Based on example in R package
    # Gaussian example

    ############################################################################
    ## test Gaussian distribution gbm model
    set.seed(123)

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

    # create a bunch of missing values
    X1[sample(1:N,size=100)] <- NA
    X3[sample(1:N,size=300)] <- NA

    w <- rep(1,N)

    data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)

    # fit initial model
    gbm1 <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
                data=data,                   # dataset
                var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                distribution="gaussian",     # bernoulli, adaboost, gaussian, poisson, coxph, or
                                            # list(name="quantile",alpha=0.05) for quantile regression
                n.trees=2000,                 # number of trees
                shrinkage=0.005,             # shrinkage or learning rate, 0.001 to 0.1 usually work
                interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc
                bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                train.fraction = 0.5,        # fraction of data for training, first train.fraction*N used for training
                n.minobsinnode = 10,         # minimum number of obs needed in each node
                keep.data=TRUE,
                cv.folds=10)                 # do 10-fold cross-validation

    # Get best model
    best.iter <- gbm.perf(gbm1,method="cv", plot.it=FALSE)   # returns cv estimate of best number of trees

    set.seed(223)
    # make some new data
    N <- 1000
    X1 <- runif(N)
    X2 <- 2*runif(N)
    X3 <- factor(sample(letters[1:4],N,replace=TRUE))
    X4 <- ordered(sample(letters[1:6],N,replace=TRUE))
    X5 <- factor(sample(letters[1:3],N,replace=TRUE))
    X6 <- 3*runif(N)
    mu <- c(-1,0,1,2)[as.numeric(X3)]

    # Actual underlying signal
    Y <- X1**1.5 + 2 * (X2**.5) + mu

    # Want to see how close predictions are to the underlying signal; noise would just interfere with this
    # Y <- Y + rnorm(N,0,sigma) 
    data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)

    # predict on the new data using "best" number of trees
    f.predict <- predict(gbm1,data2,best.iter) # f.predict will be on the canonical scale (logit,log,etc.)

    # Base the validation tests on observed discrepancies
    RUnit::checkTrue(abs(mean(data2$Y-f.predict)) < 0.01, msg="Gaussian absolute error within tolerance")
    RUnit::checkTrue(sd(data2$Y-f.predict) < sigma , msg="Gaussian squared erroor within tolerance")

    ############################################################################
    ## test coxph distribution gbm model
    ## COX PROPORTIONAL HAZARDS REGRESSION EXAMPLE

    cat("Running cox proportional hazards regression example.\n")
    # create some data
    set.seed(2)
    N <- 3000
    X1 <- runif(N)
    X2 <- runif(N)
    X3 <- factor(sample(letters[1:4],N,replace=T))
    mu <- c(-1,0,1,2)[as.numeric(X3)]

    f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)
    tt.surv <- rexp(N,exp(f))
    tt.cens <- rexp(N,0.5)
    delta <- as.numeric(tt.surv <= tt.cens)
    tt <- apply(cbind(tt.surv,tt.cens),1,min)

    # throw in some missing values
    X1[sample(1:N,size=100)] <- NA
    X3[sample(1:N,size=300)] <- NA

    # random weights if you want to experiment with them
    w <- rep(1,N)

    data <- data.frame(tt=tt,delta=delta,X1=X1,X2=X2,X3=X3)

    # fit initial model
    gbm1 <- gbm(Surv(tt,delta)~X1+X2+X3,   # formula
                data=data,                 # dataset
                weights=w,
                var.monotone=c(0,0,0),     # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                distribution="coxph",
                n.trees=3000,              # number of trees
                shrinkage=0.001,           # shrinkage or learning rate, 0.001 to 0.1 usually work
                interaction.depth=3,       # 1: additive model, 2: two-way interactions, etc
                bag.fraction = 0.5,        # subsampling fraction, 0.5 is probably best
                train.fraction = 0.5,      # fraction of data for training, first train.fraction*N used for training
                cv.folds = 5,              # do 5-fold cross-validation
                n.minobsinnode = 10,       # minimum total weight needed in each node
                keep.data = TRUE)

    best.iter <- gbm.perf(gbm1,method="test", plot.it=FALSE) # returns test set estimate of best number of trees

    # make some new data
    set.seed(2)
    N <- 1000
    X1 <- runif(N)
    X2 <- runif(N)
    X3 <- factor(sample(letters[1:4],N,replace=T))
    mu <- c(-1,0,1,2)[as.numeric(X3)]

    f <- 0.5*sin(3*X1 + 5*X2^2 + mu/10)  # -0.5 <= f <= 0.5 via sin fn.
    tt.surv <- rexp(N,exp(f))
    tt.cens <- rexp(N,0.5)

    data2 <- data.frame(tt=apply(cbind(tt.surv,tt.cens),1,min),
                        delta=as.numeric(tt.surv <= tt.cens),
                        f=f,
                        X1=X1,X2=X2,X3=X3)

    # predict on the new data using "best" number of trees
    # f.predict will be on the canonical scale (logit,log,etc.)
    f.predict <- predict(gbm1, newdata = data2, n.trees = best.iter)

    #plot(data2$f,f.predict)
    # Use observed sd
    RUnit::checkTrue(sd(data2$f - f.predict) < 0.4, msg="Coxph: squared error within tolerance")

    ############################################################################
    ## Test bernoulli distribution gbm model
    
    set.seed(1)

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
                n.trees=3000,              # number of trees
                shrinkage=0.001,           # shrinkage or learning rate, 0.001 to 0.1 usually work
                interaction.depth=3,       # 1: additive model, 2: two-way interactions, etc
                bag.fraction = 0.5,        # subsampling fraction, 0.5 is probably best
                train.fraction = 0.5,      # fraction of data for training, first train.fraction*N used for training
                cv.folds=5,                # do 5-fold cross-validation
                n.minobsinnode = 10)       # minimum total weight needed in each node

    best.iter.test <- gbm.perf(gbm1,method="test", plot.it=FALSE) # returns test set estimate of best number of trees

    best.iter <- best.iter.test

    # make some new data
    set.seed(2)
    N <- 1000
    X1 <- runif(N)
    X2 <- runif(N)
    X3 <- factor(sample(letters[1:4],N,replace=T))
    mu <- c(-1,0,1,2)[as.numeric(X3)]

    p <- 1/(1+exp(-(sin(3*X1) - 4*X2 + mu)))
    Y <- rbinom(N,1,p)
    data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)

    # predict on the new data using "best" number of trees
    # f.predict will be on the canonical scale (logit,log,etc.)
    f.1.predict <- predict(gbm1,data2, n.trees=best.iter.test)

    # compute quantity prior to transformation
    f.new = sin(3*X1) - 4*X2 + mu

    # Base the validation tests on observed discrepancies
    RUnit::checkTrue(sd(f.new - f.1.predict) < 1.0 )
    
    invisible()
}

################################################################################
########################### test.relative.influence() ##########################
###########################                           ##########################


#' @export
test.relative.influence <- function(){
    # Test that relative.influence really does pick out the true predictors
    set.seed(1234)
    X1 <- matrix(nrow=1000, ncol=50)
    X1 <- apply(X1, 2, function(x) rnorm(1000)) # Random noise
    X2 <- matrix(nrow=1000, ncol=5)
    X2 <- apply(X2, 2, function(x) c(rnorm(500), rnorm(500, 3))) # Real predictors
    cls <- rep(c(0, 1), ea=500) # Class 
    X <- data.frame(cbind(X1, X2, cls))
    mod <- gbm(cls ~ ., data= X, n.trees=1000, cv.folds=5,
                shrinkage=.01, interaction.depth=2)
    ri <- rev(sort(relative.influence(mod)))
    wh <- names(ri)[1:5]
    res <- sum(wh %in% paste("V", 51:55, sep = ""))
    RUnit::checkEqualsNumeric(res, 5, msg="Testing relative.influence identifies true predictors")
}

################################################################################
################################ validate.gbm() ################################
################################                ################################


#' @export
validate.gbm <- function () {
    wh <- (1:length(search()))[search() == "package:gbm"]
    tests <- objects(wh)[substring(objects(wh), 1, 5) == "test."]
    
    # Create temporary directory to put tests into
    sep <- if (.Platform$OS.type == "windows") "\\" else "/"
    
    dir <- file.path(tempdir(), "gbm.tests", fsep = sep)
    
    dir.create(dir)
    
    for (i in 1:length(tests)) {
        str <- paste(dir, sep, tests[i], ".R", sep = "")
        dump(tests[i], file = str)
    }
    res <- RUnit::defineTestSuite("gbm", dirs = dir, testFuncRegexp = "^test.+", 
                                  testFileRegexp = "*.R")
    cat("Running gbm test suite.\nThis will take some time...\n\n")
    RUnit::runTestSuite(res)
}

