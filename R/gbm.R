.First.lib <- function(lib, pkg)
{
     library.dynam("gbm", pkg, lib)
     require(survival)
     require(mgcv)
     require(modreg)
     require(lattice)     
}


predict.gbm <- function(object,newdata,n.trees,...)
{
   X <- model.frame(delete.response(object$Terms),
                                    newdata,
                                    na.action=na.pass)
   cRows <- dim(X)[1]
   cCols <- dim(X)[2]

   for(i in 1:cCols)
   {
      if(is.factor(X[,i]))
      {
         j <- match(levels(X[,i]),object$var.levels[[i]])
         if(any(is.na(j)))
         {
            stop(paste("New levels for variable ",
                        object$var.names[i],": ",
                        levels(X[,i])[is.na(j)],sep=""))
         }
         X[,i] <- as.numeric(X[,i])-1
      }
   }   

   X <- as.vector(unlist(X))
   if(missing(n.trees) || (n.trees > object$n.trees))
   {
      n.trees <- object$n.trees
      warning("Number of trees not specified or exceeded number fit so far. Using ",n.trees,".")
   }

   predF <- .Call("gbm_pred",
                  X=as.double(X),
                  cRows=as.integer(cRows),
                  cCols=as.integer(cCols),
                  n.trees=as.integer(n.trees),
                  initF=object$initF,
                  trees=object$trees,
                  c.split=object$c.split,
                  var.type=as.integer(object$var.type),
                  PACKAGE = "gbm")

   return(predF)
}


plot.gbm <- function(x,
                     i.var=1,
                     n.trees=x$n.trees,
                     continuous.resolution=100,
                     return.grid = FALSE,
                     ...)
{

   if((min(i.var)<1) || (max(i.var)>length(x$var.names)))
   {
      warning("i.var must be between 1 and ",length(x$var.names))
   }
   if(n.trees > x$n.trees)
   {
      warning(paste("n.trees exceeds the number of trees in the model, ",x$n.trees,
                    ". Plotting using ",x$n.trees," trees.",sep=""))
      n.trees <- x$n.trees
   }

   if(length(i.var) > 3)
   {
     warning("plot.gbm creates up to 3-way interaction plots.\nplot.gbm will only return the plotting data structure.") 
     return.grid = TRUE
   }

   # generate grid to evaluate gbm model
   grid.levels <- vector("list",length(i.var))
   for(i in 1:length(i.var))
   {
      # continuous
      if(is.numeric(x$var.levels[[i.var[i]]]))
      {
         grid.levels[[i]] <- seq(min(x$var.levels[[i.var[i]]]),
                                 max(x$var.levels[[i.var[i]]]),
                                 length=continuous.resolution)
      }
      # categorical or ordered
      else
      {
         grid.levels[[i]] <- as.numeric(factor(x$var.levels[[i.var[i]]],
                                               levels=x$var.levels[[i.var[i]]]))-1
      }
   }
   X <- expand.grid(grid.levels)
   names(X) <- paste("x",1:length(i.var),sep="")
   
   # evaluate at each data point
   X$y <- .Call("gbm_plot",
                X=as.double(data.matrix(X)),
                cRows=as.integer(nrow(X)),
                cCols=as.integer(ncol(X)),
                i.var=as.integer(i.var-1),
                n.trees=as.integer(n.trees),
                initF=as.double(x$initF),
                trees=x$trees,
                c.splits=x$c.splits,
                var.type=as.integer(x$var.type),
                PACKAGE = "gbm")
                
   # transform categorical variables back to factors
   f.factor <- rep(FALSE,length(i.var))
   for(i in 1:length(i.var))
   {
      if(!is.numeric(x$var.levels[[i.var[i]]]))
      {
         X[,i] <- factor(x$var.levels[[i.var[i]]][X[,i]+1],
                         levels=x$var.levels[[i.var[i]]])
         f.factor[i] <- TRUE
      }
   }
   
   if(return.grid)
   {
      names(X)[1:length(i.var)] <- x$var.names[i.var]
      return(X)
   }   


   # create the plots
   if(length(i.var)==1)
   {
      if(!f.factor)
      {
         j <- order(X)
         plot(X$x1,X$y,
              type="l",
              xlab=x$var.names[i.var],
              ylab=paste("f(",x$var.names[i.var],")",sep=""),...)
      }
      else
      {
         plot(X$x1,
              X$y,
              xlab=x$var.names[i.var],
              ylab=paste("f(",x$var.names[i.var],")",sep=""),
              ...)
      }
   }
   else if(length(i.var)==2)
   {
      if(!f.factor[1] && !f.factor[2])
      {
         levelplot(y~x1*x2,data=X,
                     xlab=x$var.names[i.var[1]],
                     ylab=x$var.names[i.var[2]],...)
      }
      else if(f.factor[1] && !f.factor[2])
      {
         xyplot(y~x2|x1,data=X,
                xlab=x$var.names[i.var[1]],
                ylab=paste("f(",x$var.names[i.var[1]],",",x$var.names[i.var[2]],")",sep=""),
                type="l",
                ...)
      }
      else if(!f.factor[1] && f.factor[2])
      {
         xyplot(y~x1|x2,data=X,
                xlab=x$var.names[i.var[2]],
                ylab=paste("f(",x$var.names[i.var[1]],",",x$var.names[i.var[2]],")",sep=""),
                type="l",
                ...)
      }
      else
      {
         stripplot(y~x1|x2,data=X,
                   xlab=x$var.names[i.var[2]],
                   ylab=paste("f(",x$var.names[i.var[1]],",",x$var.names[i.var[2]],")",sep=""),
                   ...)
      }
   }
   else if(length(i.var)==3)
   {
      i <- order(f.factor)
      X.new <- X[,i]
      X.new$y <- X$y
      names(X.new) <- names(X)
      
      # 0 factor, 3 continuous
      if(sum(f.factor)==0)
      {
         X.new$x3 <- equal.count(X.new$x3)
         levelplot(y~x1*x2|x3,data=X.new,
                     xlab=x$var.names[i.var[i[1]]],
                     ylab=x$var.names[i.var[i[2]]],...)
      }
      # 1 factor, 2 continuous
      else if(sum(f.factor)==1)
      {
         levelplot(y~x1*x2|x3,data=X.new,
                     xlab=x$var.names[i.var[i[1]]],
                     ylab=x$var.names[i.var[i[2]]],...)
      }
      # 2 factors, 1 continuous
      else if(sum(f.factor)==2)
      {
         xyplot(y~x1|x2*x3,data=X.new,
                type="l",
                xlab=x$var.names[i.var[i[1]]],
                ylab=paste("f(",paste(x$var.names[i.var[1:3]],collapse=","),")",sep=""),
                ...)
      }
      # 3 factors, 0 continuous
      else if(sum(f.factor)==3)
      {
         stripplot(y~x1|x2*x3,data=X.new,
                   xlab=x$var.names[i.var[i[1]]],
                   ylab=paste("f(",paste(x$var.names[i.var[1:3]],collapse=","),")",sep=""),
                   ...)
      }
   }
}


gbm.more <- function(object,
                     n.new.trees = 100,
                     data = NULL,
                     weights = NULL)
{
   if(is.null(object$data) && is.null(data))
   {
      stop("keep.data was set to FALSE on original gbm call and argument 'data' is NULL")
   }
   else if(is.null(object$data))
   {
      m <- eval(object$m, parent.frame())
      Terms <- attr(m, "terms")
      a <- attributes(Terms)
      
      Y <- as.vector(model.extract(m, response))
      Offset <- model.extract(m,offset)
      X <- model.frame(delete.response(Terms),
                       data,
                       na.action=na.pass)
      
      w <- weights
      if(length(w)==0) w <- rep(1, nrow(X))
      w <- w*length(w)/sum(w) # normalize to N
      
      if(is.null(Offset) || (Offset==0))
      {
         Offset <- NA
      }
      Misc <- NA

      if(object$distribution == "coxph")
      {
         Misc <- as.numeric(Y)[-(1:cRows)]
         Y <- as.numeric(Y)[1:cRows]
   
         # reverse sort the failure times to compute risk sets on the fly
         i.train <- order(-Y[1:object$nTrain])
         i.test <- order(-Y[(object$nTrain+1):cRows]) + object$nTrain
         i.timeorder <- c(i.train,i.test)
   
         Y <- Y[i.timeorder]
         Misc <- Misc[i.timeorder]
         X <- X[i.timeorder,]
         w <- w[i.timeorder]
         if(!is.na(Offset)) Offset <- Offset[i.timeorder]
         object$fit <- object$fit[i.timeorder]
      }

      # create index upfront... subtract one for 0 based order
      if(ncol(X) > 1)
      {
         X.order <- apply(X[1:object$nTrain,],2,order,na.last=FALSE)-1
      }
      else
      {
         X.order <- order(X[1:object$nTrain,],na.last=FALSE)-1
      }
      X <- data.matrix(X)
      cRows <- nrow(X)
      cCols <- ncol(X)
   }
   else
   {
      Y           <- object$data$Y
      X           <- object$data$X
      X.order     <- object$data$X.order
      Offset      <- object$data$Offset
      Misc        <- object$data$Misc
      w           <- object$data$w
      cRows <- length(Y)
      cCols <- length(X)/cRows
      if(object$distribution == "coxph")
      {
         i.timeorder <- object$data$i.timeorder
         object$fit <- object$fit[i.timeorder]
      }
   }
   
   X <- as.vector(X)
   
   gbm.fit <- .Call("gbm",
                    Y = as.double(Y),
                    Offset = as.double(Offset),
                    X = as.double(X),
                    X.order = as.integer(X.order),
                    weights = as.double(w),
                    Misc = as.double(Misc),
                    cRows = as.integer(cRows),
                    cCols = as.integer(cCols),
                    var.type = as.integer(object$var.type),
                    var.monotone = as.integer(object$var.monotone),
                    distribution = as.character(object$distribution),
                    n.trees = as.integer(n.new.trees),
                    interaction.depth = as.integer(object$interaction.depth),
                    n.minobsinnode = as.integer(object$n.minobsinnode),
                    shrinkage = as.double(object$shrinkage),
                    bag.fraction = as.double(object$bag.fraction),
                    nTrain = as.integer(object$nTrain),
                    fit.old = as.double(object$fit),
                    n.cat.splits.old = length(object$c.splits),
                    n.trees.old = as.integer(object$n.trees),
                    PACKAGE = "gbm")
   names(gbm.fit) <- c("initF","fit","train.error","valid.error",
                       "oobag.improve","trees","c.splits")
   
   gbm.fit$initF         <- object$initF 
   gbm.fit$train.error   <- c(object$train.error, gbm.fit$train.error)
   gbm.fit$valid.error   <- c(object$valid.error, gbm.fit$valid.error)
   gbm.fit$oobag.improve <- c(object$oobag.improve, gbm.fit$oobag.improve)
   gbm.fit$trees         <- c(object$trees, gbm.fit$trees)
   gbm.fit$c.splits      <- c(object$c.splits, gbm.fit$c.splits)
   
   gbm.fit$n.trees        <- length(gbm.fit$trees)
   gbm.fit$distribution   <- object$distribution
   gbm.fit$train.fraction <- object$train.fraction    
   gbm.fit$shrinkage      <- object$shrinkage
   gbm.fit$bag.fraction   <- object$bag.fraction
   gbm.fit$var.type       <- object$var.type
   gbm.fit$var.monotone   <- object$var.monotone
   gbm.fit$var.names      <- object$var.names
   gbm.fit$interaction.depth       <- object$interaction.depth
   gbm.fit$n.minobsinnode <- object$n.minobsinnode
   gbm.fit$nTrain         <- object$nTrain
   gbm.fit$Terms          <- object$Terms
   gbm.fit$var.levels     <- object$var.levels
      
   if(object$distribution == "coxph")
   {
      gbm.fit$fit[i.timeorder] <- gbm.fit$fit
   }
   if(!is.null(object$data))
   {
      gbm.fit$data <- object$data
   }
   else
   {
      gbm.fit$data <- NULL
      gbm.fit$m <- object$m
   }

   class(gbm.fit) <- "gbm"
   return(gbm.fit)
}


gbm <- function(formula = formula(data),
                distribution = "bernoulli",
                data = list(),
                weights,
                var.monotone = NULL,
                n.trees = 100,
                interaction.depth = 1,
                n.minobsinnode = 10,
                shrinkage = 0.1,
                bag.fraction = 1.0,
                train.fraction = 0.5,
                keep.data = TRUE)
{
   call <- match.call()
   m <- match.call(expand = FALSE)
   m$distribution <- m$var.monotone <- m$n.trees <- NULL
   m$interaction.depth <- m$n.minobsinnode <- m$shrinkage <- NULL
   m$bag.fraction <- m$train.fraction <- m$keep.data <- NULL
   m[[1]] <- as.name("model.frame")
   m$na.action <- na.pass
   m.keep <- m
   
   m <- eval(m, parent.frame())
   
   Terms <- attr(m, "terms")
   a <- attributes(Terms)
   
   Y <- model.extract(m, response)
   
   Offset <- model.extract(m,offset)
   X <- model.frame(delete.response(Terms),
                    data,
                    na.action=na.pass)
   
   cRows <- nrow(X)
   cCols <- ncol(X)
   
   # check dataset size
   if(cRows*train.fraction*bag.fraction <= 2*n.minobsinnode+1)
   {
      stop("The dataset size is too small or subsampling rate is too large: cRows*train.fraction*bag.fraction <= n.minobsinnode")
   }
   if(interaction.depth < 1)
   {
      stop("interaction.depth must be at least 1.")
   }
   
   w <- model.extract(m, weights)
   if(length(w)==0) w <- rep(1, nrow(m))
   else if(any(w < 0)) stop("negative weights not allowed")
   
   w <- w*length(w)/sum(w) # normalize to N
   var.names <- a$term.labels
   # get the character name of the response variable
   response.name <- dimnames(attr(terms(formula),"factors"))[[1]][1]
   if(any(is.na(Y))) stop("Missing values are not allowed in the response, ",response.name)
   
   if(is.null(Offset) || (Offset==0))
   {
      Offset <- NA
   }
   Misc <- NA
   
   # setup variable types
   var.type <- rep(0,cCols)
   var.levels <- vector("list",cCols)
   for(i in 1:length(var.type))
   {
      if(is.ordered(X[,i]))
      {
         var.levels[[i]] <- levels(X[,i])
         X[,i] <- as.numeric(X[,i])-1
         var.type[i] <- 0      
      }
      else if(is.factor(X[,i]))
      {
         var.levels[[i]] <- levels(X[,i])
         X[,i] <- as.numeric(X[,i])-1
         var.type[i] <- max(X[,i],na.rm=TRUE)+1
      }
      else if(is.numeric(X[,i]))
      {
         var.levels[[i]] <- quantile(X[,i],prob=(0:10)/10,na.rm=TRUE)
      }
      else
      {
         stop("variable ",i,": ",var.names[i]," if not of type numeric, ordered, or factor.")
      }
      
      # check for some variation in each variable
      if(length(unique(var.levels[[i]])) == 1)
      {
         stop("variable ",i,": ",var.names[i]," has no variation.")
      }
   }
   
   nTrain <- as.integer(train.fraction*cRows)
      
   supported.distributions <- 
      c("bernoulli","gaussian","poisson","adaboost","laplace","coxph")
   # check potential problems with the distributions
   if(!is.element(distribution,supported.distributions))
   {
      stop("Distribution ",distribution," is not supported")
   }
   if((distribution == "bernoulli") && !all(is.element(Y,0:1)))
   {
      stop("Bernoulli requires the response to be in {0,1}")
   }
   if((distribution == "poisson") && any(Y<0))
   {
      stop("Poisson requires the response to be positive")
   }
   if((distribution == "poisson") && any(Y != trunc(Y)))
   {
      stop("Poisson requires the response to be a positive integer")
   }
   if((distribution == "adaboost") && !all(is.element(Y,0:1)))
   {
      stop("This version of AdaBoost requires the response to be in {0,1}")
   }
   if((distribution == "laplace") && (length(unique(w)) > 1))
   {
      stop("This version of gbm for the Laplace loss lacks a weighted median. For now the weights must be constant.")
   }   
   if(distribution == "coxph")
   {
      if(class(Y)!="Surv")
      {
         stop("Outcome must be a survival object Surv(time,failure)")
      }
      if(attr(Y,"type")!="right")
      {
         stop("gbm() currently only handles right censored observations")
      }
      Misc <- as.numeric(Y)[-(1:cRows)]
      Y <- as.numeric(Y)[1:cRows]
   
      # reverse sort the failure times to compute risk sets on the fly
      i.train <- order(-Y[1:nTrain])
      n.test <- cRows - nTrain
      i.test <- order(-Y[(nTrain+1):cRows]) + nTrain
      i.timeorder <- c(i.train,i.test)

      Y <- Y[i.timeorder]
      Misc <- Misc[i.timeorder]
      X <- X[i.timeorder,]
      w <- w[i.timeorder]
      if(!is.na(Offset)) Offset <- Offset[i.timeorder]

   }

   # create index upfront... subtract one for 0 based order
   if(ncol(X) > 1)
   {
      X.order <- apply(X[1:nTrain,],2,order,na.last=FALSE)-1
   }
   else
   {
      X.order <- order(X[1:nTrain,],na.last=FALSE)-1   
   }
   X <- as.vector(data.matrix(X))
   predF <- rep(0,length(Y))
   train.error <- rep(0,n.trees)
   valid.error <- rep(0,n.trees)
   oobag.improve <- rep(0,n.trees)
   
   if(is.null(var.monotone)) var.monotone <- rep(0,cCols)
   else if(length(var.monotone)!=cCols) 
   {
      stop("Length of var.monotone != number of predictors")
   }
   else if(!all(is.element(var.monotone,-1:1))) 
   {
      stop("var.monotone must be -1, 0, or 1")
   }
   fError <- FALSE

   gbm.fit <- .Call("gbm",
                    Y=as.double(Y),
                    Offset=as.double(Offset),
                    X=as.double(X),
                    X.order=as.integer(X.order),
                    weights=as.double(w),
                    Misc=as.double(Misc),
                    cRows=as.integer(cRows),
                    cCols=as.integer(cCols),
                    var.type=as.integer(var.type),
                    var.monotone=as.integer(var.monotone),
                    distribution=as.character(distribution),
                    n.trees=as.integer(n.trees),
                    interaction.depth=as.integer(interaction.depth),
                    n.minobsinnode=as.integer(n.minobsinnode),
                    shrinkage=as.double(shrinkage),
                    bag.fraction=as.double(bag.fraction),
                    nTrain=as.integer(nTrain),
                    fit.old=as.double(NA),
                    n.cat.splits.old=as.integer(0),
                    n.trees.old=as.integer(0),
                    PACKAGE = "gbm")
   names(gbm.fit) <- c("initF","fit","train.error","valid.error",
                       "oobag.improve","trees","c.splits")
                       
   gbm.fit$bag.fraction <- bag.fraction
   gbm.fit$distribution <- distribution
   gbm.fit$interaction.depth <- interaction.depth
   gbm.fit$n.minobsinnode <- n.minobsinnode
   gbm.fit$n.trees <- length(gbm.fit$trees)
   gbm.fit$nTrain <- nTrain
   gbm.fit$response.name <- response.name
   gbm.fit$shrinkage <- shrinkage
   gbm.fit$Terms <- Terms
   gbm.fit$train.fraction <- train.fraction    
   gbm.fit$var.levels <- var.levels
   gbm.fit$var.monotone <- var.monotone
   gbm.fit$var.names <- var.names
   gbm.fit$var.type <- var.type
      
   if(distribution == "coxph")
   {
      gbm.fit$fit[i.timeorder] <- gbm.fit$fit
   }

   if(keep.data)
   {
      if(distribution == "coxph")
      {
         # put the observations back in order
         gbm.fit$data <- list(Y=Y,X=X,X.order=X.order,Offset=Offset,Misc=Misc,w=w,
                              i.timeorder=i.timeorder)
      }
      else
      {
         gbm.fit$data <- list(Y=Y,X=X,X.order=X.order,Offset=Offset,Misc=Misc,w=w)
      }
   }
   else
   {
      gbm.fit$data <- NULL
      gbm.fit$m <- m.keep
   }

   class(gbm.fit) <- "gbm"
   return(gbm.fit)
}


gbm.perf <- function(object,
            plot.it=TRUE,
            oobag.curve=TRUE,
            overlay=TRUE,
            method = c("OOB","test")[1])
{
     smoother <- NULL
     
     if((method == "OOB") || oobag.curve)
     {
          x <- 1:object$n.trees
          smoother <- loess(object$oobag.improve~x,
                            enp.target=min(max(4,length(x)/10),50))
          smoother$y <- smoother$fitted
          smoother$x <- x

          best.iter.oob <- x[order(-cumsum(smoother$y))[1]]
     }

     if(method == "test")
     {
          if((object$distribution == "gaussian") ||
             (object$distribution == "adaboost") ||
             (object$distribution == "laplace"))
          {
                best.iter.test <- order(object$valid.error)[1]
          }
          else
          {
                best.iter.test <- order(-object$valid.error)[1]      
          }
     }

     if(method == "OOB")
     {
          best.iter <- best.iter.oob
          cat("Out of bag iteration estimate:",best.iter.oob,"\n")
     }
     else if(method == "test")
     {
          best.iter <- best.iter.test
          cat("Test set iteration estimate:",best.iter.test,"\n")
     }
     else stop("method must be OOB or test")

     if(plot.it)
     {
          par(mfrow=c(1,1),mar=c(5,4,4,4)+.1)
          ylab <- switch(substring(object$distribution,1,2),
                                     ga="Squared error loss",
                                     be="Bernoulli log-likelihood",
                                     po="Poisson log-likelihood",
                                     ad="AdaBoost exponential bound",
                                     co="Cox partial log-likelihood")
          if(object$train.fraction==1)
          {
                ylim <- range(object$train.error)
          }
          else
          {
                ylim <- range(object$train.error,object$valid.error)
          }
          
          plot(object$train.error,
                ylim=ylim,
                type="l",
                xlab="Iteration",ylab=ylab)
                
          if(object$train.fraction!=1)
          {
                lines(object$valid.error,col="red")
          }
          
          if(oobag.curve)
          {
                if(overlay)
                {
                     par(new=TRUE)
                     plot(smoother$x,
                             cumsum(smoother$y),
                             col="blue",
                             type="l",
                             xlab="",ylab="",
                             axes=FALSE)
                     axis(4,srt=0)
                     at <- mean(range(smoother$y))
                     mtext(paste("OOB improvement in",ylab),side=4,srt=270,line=2)
                     abline(h=0,col="blue",lwd=2)
                     if(!is.na(best.iter)) abline(v=best.iter,col="blue",lwd=2)
                }
     
                plot(object$oobag.improve,type="l",
                          xlab="Iteration",
                          ylab=paste("OOB change in",ylab))
                lines(smoother,col="red",lwd=2)
                abline(h=0,col="blue",lwd=1)

                abline(v=best.iter,col="blue",lwd=1)
          }
     }
     
     return(best.iter)
}


summary.gbm <- function(object,
                        cBars=length(object$var.names),
                        n.trees=object$n.trees,
                        plotit=TRUE,
                        order=TRUE,
                        ...)
{
   get.rel.inf <- function(obj)
   {
      lapply(split(obj[[6]],obj[[1]]),sum) # 6 - Improvement, 1 - var name
   }

   if(n.trees < 1)
   {
      stop("n.trees must be greater than 0.")
   }
   if(n.trees > object$n.trees)
   {
      warning("Exceeded total number of GBM terms. Results use",object$n.trees,"terms.\n")
      n.trees <- object$n.trees
   }

   temp <- unlist(lapply(object$trees[1:n.trees],get.rel.inf))
   rel.inf.compact <- unlist(lapply(split(temp,names(temp)),sum))
   rel.inf.compact <- rel.inf.compact[names(rel.inf.compact)!="-1"]

   # rel.inf.compact excludes those variable that never entered the model
   # insert 0's for the excluded variables
   rel.inf <- rep(0,length(object$var.names))
   i <- as.numeric(names(rel.inf.compact))+1
   rel.inf[i] <- rel.inf.compact

   if(order)
   {
      i <- order(-rel.inf)
   }
   else
   {
      i <- 1:length(rel.inf)
   }
   if(cBars==0) cBars <- min(10,length(object$var.names))
   if(cBars>length(object$var.names)) cBars <- length(object$var.names)
   
   rel.inf <- 100*rel.inf/sum(rel.inf)
   
   if(plotit)
   {
      barplot(rel.inf[i[cBars:1]],
              horiz=TRUE,
              col=rainbow(cBars,start=3/6,end=4/6),
              names=object$var.names[i[cBars:1]],
              xlab="Relative influence",...)
   }
   return(data.frame(var=object$var.names[i],
                     rel.inf=rel.inf[i]))
}


quantile.rug <- function(x,prob=(0:10)/10,...)
{
     quants <- quantile(x[!is.na(x)],prob=prob)
     if(length(unique(quants)) < length(prob))
     {
          quants <- jitter(quants)
     }
     rug(quants,...)
}

calibrate.plot <- function(y,p,
                           distribution="bernoulli",
                           replace=TRUE,
                           line.par=list(col="black"),
                           shade.col="lightyellow",
                           shade.density=NULL,
                           rug.par=list(side=1),...)
{
   data <- data.frame(y=y,p=p)
   
   if(distribution=="bernoulli")
   {
      family1 = binomial
   } else if(distribution=="poisson")
   {
      family1 = poisson
   } else
   {
      family1 = gaussian
   } 
   gam1 <- gam(y~s(p),data=data,family=family1)
   
   x <- seq(min(p),max(p),length=200)
   yy <- predict(gam1,newdata=data.frame(p=x),se.fit=TRUE,type="response")

   x <- x[!is.na(yy$fit)]
   yy$se.fit <- yy$se.fit[!is.na(yy$fit)]
   yy$fit <- yy$fit[!is.na(yy$fit)]      

   if(replace)
   {
      plot(0,0,
         type="n",
         xlab="Estimated probability",ylab="Actual probability",
         ...)
   }
   if(!is.na(shade.col))
   {
      se.lower <- yy$fit-2*yy$se.fit
      se.upper <- yy$fit+2*yy$se.fit
      se.lower[se.lower < 0] <- 0
      se.upper[se.upper > 1] <- 1
      polygon(c(x,rev(x),x[1]),
               c(se.lower,rev(se.upper),se.lower[1]),
               col=shade.col,
               border=NA,
               density=shade.density)
   }
   lines(x,yy$fit,col=line.par$col)
   quantile.rug(p,side=rug.par$side)
   abline(0,1,col="red")
}



pretty.gbm.tree <- function(object,i.tree=1)
{
   if((i.tree<1) || (i.tree>length(object$trees)))
   {
      stop("i.tree is out of range. Must be less than ",length(object$trees))
   }
   else
   {
      temp <- data.frame(object$trees[[i.tree]])
      names(temp) <- c("SplitVar","SplitCodePred","LeftNode",
                       "RightNode","MissingNode","ErrorReduction",
                       "Weight")
      row.names(temp) <- 0:(nrow(temp)-1)
   }
   return(temp)
}
