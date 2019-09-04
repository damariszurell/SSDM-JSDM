# Codes for: Zurell et al. "Testing species assemblage predictions from stacked and joint species distribution models" 

# Variable selection: only 5 most important, weakly correlated variables are used for subsequent models.
# Univariate variable importance is assessed using GAMs in 5-fold spatial block cross-validation design.

#----------------------------------------------------------------------------------------------------------

# load required packages

library(sperrorest)
library(dismo)
library(mgcv)

#----------------------------------------------------------------------------------------------------------

# Load data. Our datasets were called the following:
aviCH.forest.train		# birds - training data for SDM and JSDM calibration
aviCH.forest.test		# birds - test data for evaluating community predictions
datLFI.train			# trees - training data for SDM and JSDM calibration
datLFI.test				# trees - test data for evaluating community predictions

# Prepare species names and predictor names. We used the following names:
avinames.50presences 	# names of birds with at least 50 presences in study region
lfinames.50presences 	# names of trees with at least 50 presences in study region
avi.pred.names 			# names of all potential predictors for birds
lfi.pred.names 			# names of all potential predictors for trees


#----------------------------------------------------------------------------------------------------------

# prepare blocks for spatially blocked cross-validation
avi.cv.tiles <- partition_tiles(data=aviCH.forest.train[,c('x','y')], nsplit=c(3,2),reassign=T,min_n=200)[[1]]
avi.cv.tilesL <- numeric(nrow(aviCH.forest.train))
for (ti in seq_len(length(avi.cv.tiles))){
	avi.cv.tilesL[avi.cv.tiles[[ti]]$test] <- ti
}


lfi.cv.tiles <- partition_tiles(data= datLFI.train[,c('X','Y')], coords=c('X','Y'),nsplit=c(3,2),reassign=T,min_n=570)[[1]]
lfi.cv.tilesL <- numeric(nrow(datLFI.train))
for (ti in seq_len(length(lfi.cv.tiles))){
	lfi.cv.tilesL[lfi.cv.tiles[[ti]]$test] <- ti
}


#----------------------------------------------------------------------------------------------------------

# Compute univariate variable importance


# helper function for spline formulas
gam.formula <- function(sp, env) {
	formula(paste(sp,'~ s(',env,',k=4)'))
}


# helper function to compute univariate GAMs for different folds of spatial block cv and
# calculate explained deviance (d2) for cross-predictions
compute.univar.cv <- function(sp,env,df,cv.tiles,family='binomial') {
	pred.df <- data.frame(tiles=cv.tiles,pred=NA, pred.null=NA, occ=NA)
	for (n in sort(unique(cv.tiles))) {
		train.df <- df[!cv.tiles ==n,]
		test.df <- df[cv.tiles==n,]
		
		m1 <- gam(formula= gam.formula(sp, env), data=train.df, family=family)
		pred.df[pred.df$tiles==n,'pred'] <- predict(m1,newdata=test.df,type='response')
		pred.df[pred.df$tiles==n,'pred.null'] <- rep(mean(train.df[, sp]), length=nrow(test.df))
		pred.df[pred.df$tiles==n,'occ'] <- test.df[,sp]
	}	
	d2 <- 1 - (calc.deviance(pred.df[,'occ'], pred.df[,'pred'], family=family) / calc.deviance(pred.df[,'occ'], pred.df[,'pred.null'], family=family))
	ifelse(d2<0,0,d2)
}

# compute d2 for cross-predictions
avi.univar.cv.d2 <- data.frame(sapply(avinames.50presences, FUN=function(sp) {sapply(avi.pred.names, FUN= function(env) {compute.univar.cv(sp,env, aviCH.forest.train, avi.cv.tilesL)}) }))
names(avi.univar.cv.d2) <- avinames.50presences
rownames (avi.univar.cv.d2) <- avi.pred.names

lfi.univar.cv.d2 <- data.frame(sapply(lfinames.50presences, FUN=function(sp) {sapply(lfi.pred.names, FUN= function(env) {compute.univar.cv(sp,env, datLFI.train, lfi.cv.tilesL)}) }))
names(lfi.univar.cv.d2) <- lfinames.50presences
rownames (lfi.univar.cv.d2) <- lfi.pred.names


#----------------------------------------------------------------------------------------------------------

# select most important, weakly correlated variables


# select07, adapted from Dormann et al. (2013) Ecography 36: 27-46.
select07 <- function(imp, X, threshold=0.7, method="spearman")
{
    cm <- cor(X, method=method)
    sort.imp <- colnames(X)[order(imp,decreasing=T)]
    
    pairs <- which(abs(cm)>= threshold, arr.ind=T) # identifies correlated variable pairs
    index <- which(pairs[,1]==pairs[,2])           # removes entry on diagonal
    pairs <- pairs[-index,]                        # -"-

    exclude <- NULL
    for (i in seq_len(length(sort.imp)))
    {
    if ((sort.imp[i] %in% row.names(pairs))&
        ((sort.imp[i] %in% exclude)==F)) {
    cv<-cm[setdiff(row.names(cm),exclude),sort.imp[i]]
    cv<-cv[setdiff(names(cv),sort.imp[1:i])]
    exclude<-c(exclude,names(which((abs(cv)>=threshold)))) }
    }

	sort.imp[!(sort.imp %in% unique(exclude)),drop=F]
}


# select07 based on importance provided by d2 averaged over block cv
avi.select07.cv <- apply(avi.univar.cv.d2, 2,FUN=function(x){select07(x, aviCH.forest.train[, avi.pred.names])[1:5]})
lfi.select07.cv <- apply(lfi.univar.cv.d2, 2,FUN=function(x){select07(x, datLFI.train[, lfi.pred.names])[1:5]})

# select07 for global predictor set
avi.select07.cv.glob <- select07(rowSums(avi.univar.cv.d2), aviCH.forest.train[, avi.pred.names])[1:5]
lfi.select07.cv.glob <- select07(rowSums(lfi.univar.cv.d2), datLFI.train[, lfi.pred.names])[1:5]


#----------------------------------------------------------------------------------------------------------

# select the 5 most important, weakly correlated predictors for MEMs

# compute d2 for cross-predictions
aviCH.forest.train$sr <- rowSums(aviCH.forest.train[,avinames.50presences])
avi.mem.univar.cv.d2 <- sapply(avi.pred.names, FUN= function(env) {compute.univar.cv('sr',env, aviCH.forest.train, avi.cv.tilesL, family='poisson')})
avi.select07.cv.mem <- select07(avi.mem.univar.cv.d2, aviCH.forest.train[, avi.pred.names])[1:5]

datLFI.train$sr <- rowSums(datLFI.train[,lfinames.50presences])
lfi.mem.univar.cv.d2 <- sapply(lfi.pred.names, FUN= function(env) {compute.univar.cv('sr',env, datLFI.train, lfi.cv.tilesL, family='poisson')})
lfi.select07.cv.mem <- select07(lfi.mem.univar.cv.d2, datLFI.train[, lfi.pred.names])[1:5]
