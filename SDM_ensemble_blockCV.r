# Codes for: Zurell et al. "Testing species assemblage predictions from stacked and joint species distribution models" 

# Fitting SDMs and evaluating performance based on spatial block cross-validation.

#----------------------------------------------------------------------------------------------------------

# load required packages
library(mgcv)
library(gbm)
library(dismo)
library(randomForest)
library(Hmisc)

#----------------------------------------------------------------------------------------------------------

# Load data. Our datasets were called the following:
aviCH.forest.train		# birds - training data for SDM and JSDM calibration
aviCH.forest.test		# birds - test data for evaluating community predictions
datLFI.train			# trees - training data for SDM and JSDM calibration
datLFI.test				# trees - test data for evaluating community predictions
avi.cv.tilesL			# birds - tiles for 5-fold spatial cross-validation
lfi.cv.tilesL			# trees - tiles for 5-fold spatial cross-validation

# Prepare species names and predictor names. We used the following names:
avinames.50presences 	# names of birds with at least 50 presences in study region
lfinames.50presences 	# names of trees with at least 50 presences in study region
avi.pred.names 			# names of all potential predictors for birds
lfi.pred.names 			# names of all potential predictors for trees


#----------------------------------------------------------------------------------------------------------

# Fit SDMs to presence-absence data

# Function to fit GLM, GAM, GBM and RF
fit.alg <- function(sp,pred,dat) {
	# GLM with linear and quadratic term:
	m.glm <- glm( as.formula(paste(sp,'~',paste(pred,paste0('+ I(',pred,'^2)'),collapse=' + '))), family='binomial', data=dat)
	
	# GAM from mgcv package with k=4
	m.gam <- gam( as.formula(paste(sp,'~',paste(paste0('s(',pred,',k=4)'),collapse=' + '))), family='binomial', data=dat)
	
	# Random Forest
	m.rf <- randomForest( x=dat[,pred], y=as.factor(dat[,sp]), ntree=1000, nodesize=20)
    
	# GBM from dismo package, with automatic adapting of learning rate
	opt.LR <- TRUE
	LR <- 0.01
	while(opt.LR){
		m.brt <- try(gbm.step(data = dat, gbm.x = pred, gbm.y = sp, family = 'bernoulli', tree.complexity = 2, bag.fraction = 0.75, learning.rate = LR, verbose=F, plot.main=F))
        if (class(m.brt) == "try-error" | class(m.brt) == "NULL"){
            LR <- LR/2
        } else
        if(m.brt$gbm.call$best.trees<1000){
			LR <- LR/2
		} else 
		if(m.brt$gbm.call$best.trees>5000){
			LR <- LR*2
		} else { 
			opt.LR <- FALSE}
		}
	
	list(list(glm=m.glm, gam=m.gam, rf=m.rf, brt=m.brt))
}	

# Function to make predictions:
make.preds <- function(model, newdata, family='binomial') {
	if(family=='binomial') {
		switch(class(model)[1],
		glm = predict(model, newdata, type='response'),
		gam = predict(model, newdata, type='response'),
		randomForest = predict(model, newdata, type= 'prob')[,2], 
		gbm = predict.gbm(model, newdata, n.trees=model$gbm.call$best.trees, type="response"))
		} else
	if(family=='poisson') {
		switch(class(model)[1],
		glm = predict(model, newdata, type='response'),
		gam = predict(model, newdata, type='response'),
		randomForest = predict(model, newdata, type= 'response'), 
		gbm = predict.gbm(model, newdata, n.trees=model$gbm.call$best.trees, type="response"))
	}
}

# Function to make cross-validated predictions:
crossval.preds <- function(model, dat, sp, pred, cv.tiles, gbm.dist='bernoulli') {
	
	# Make k-fold data partitions
	kfold <- length(unique(cv.tiles))
	
	cross.val.preds <- data.frame(row = row.names(dat), 
		cross.val.preds = numeric(length = nrow(dat))) 
	
	for(i in seq_len(kfold)){
		cv.train <- dat[cv.tiles!=i,]
		cv.test <- dat[cv.tiles ==i,]
		
		# Because we used the gbm.step() for BRTs, we need a small work-around:
		if (class(model)[1]=='gbm') {
			cv.train.gbm <- cv.train;
			names(cv.train.gbm)[names(cv.train.gbm)==sp] <- 
				model$response.name
			}

		# We update the model for the new training data
		modtmp <- switch(class(model)[1],
			glm = update(model, data=cv.train),
			gam = update(model, data=cv.train),
			randomForest = update(model, data=cv.train),				
			gbm = gbm(model$call, distribution= gbm.dist, 
				data=cv.train.gbm[,c(pred,model$response.name)], 
				n.trees=model$gbm.call$best.trees))
		
		# We make predictions for k-fold:
		if (class(model)[1]=='gbm') {
			cross.val.preds[which(cv.tiles ==i),2] <- 
				predict.gbm(modtmp, cv.test[, pred], 
					n.trees=model$gbm.call$best.trees, type="response")
		} else {
			if(gbm.dist=='bernoulli') {
				cross.val.preds[which(cv.tiles ==i),2] <- make.preds(modtmp, cv.test[, pred])
				} else
			if(gbm.dist=='poisson')	{
				cross.val.preds[which(cv.tiles ==i),2] <- make.preds(modtmp, cv.test[, pred], family=gbm.dist)
				}
			}
		}
	list(cross.val.preds[,2])
	}


# Function to calculate model performances:
calc.eval <- function(dat, sp, predictions, thresh.method='MaxSens+Spec'){
	require(PresenceAbsence)
	
	# Helper functions:
	# True Skill Statistic:
	TSS = function(cmx){
		PresenceAbsence::sensitivity(cmx, st.dev=F) + 
		PresenceAbsence::specificity(cmx, st.dev=F) - 1
		}
	
	thresh.dat <- data.frame(ID=seq_len(nrow(dat)), 
		obs = dat[, sp],
		pred = predictions)
		
	thresh <- optimal.thresholds(DATA= thresh.dat)
	cmx.maxSSS <- cmx(DATA= thresh.dat, threshold=thresh[thresh$Method==thresh.method,2])
	
	data.frame(AUC = PresenceAbsence::auc(thresh.dat, st.dev=F),
		TSS = TSS(cmx.maxSSS), 
		Sens = PresenceAbsence::sensitivity(cmx.maxSSS, st.dev=F),
		Spec = PresenceAbsence::specificity(cmx.maxSSS, st.dev=F),
		thresh = thresh[thresh$Method==thresh.method,2])
}



# Fit models:
avi.models.ind <- sapply(avinames.50presences,FUN=function(sp){fit.alg(sp,pred=avi.select07.cv[,sp], dat= aviCH.forest.train)})
lfi.models.ind <- sapply(lfinames.50presences,FUN=function(sp){fit.alg(sp,pred=lfi.select07.cv[,sp], dat= datLFI.train)})
avi.models.glob <- sapply(avinames.50presences,FUN=function(sp){fit.alg(sp,pred=avi.select07.cv.glob, dat= aviCH.forest.train)})
lfi.models.glob <- sapply(lfinames.50presences,FUN=function(sp){fit.alg(sp,pred=lfi.select07.cv.glob, dat= datLFI.train)})


# Make cross-validated predictions and ensembles:
avi.cross.pred.ind <- sapply(avinames.50presences, FUN=function(sp){ algs=list(sapply(names(avi.models.ind[[sp]]), FUN=function(m){crossval.preds(avi.models.ind[[sp]][[m]], dat= aviCH.forest.train, sp=sp, pred=avi.select07.cv[,sp], cv.tiles= avi.cv.tilesL)})); algs[[1]][['mean.prob']] = rowMeans(data.frame(algs)); algs })
lfi.cross.pred.ind <- sapply(lfinames.50presences, FUN=function(sp){ algs=list(sapply(names(lfi.models.ind[[sp]]), FUN=function(m){crossval.preds(lfi.models.ind[[sp]][[m]], dat= datLFI.train, sp=sp, pred=lfi.select07.cv[,sp], cv.tiles= lfi.cv.tilesL)})); algs[[1]][['mean.prob']] = rowMeans(data.frame(algs)); algs })
avi.cross.pred.glob <- sapply(avinames.50presences, FUN=function(sp){ algs=list(sapply(names(avi.models.glob[[sp]]), FUN=function(m){crossval.preds(avi.models.glob[[sp]][[m]], dat= aviCH.forest.train, sp=sp, pred=avi.select07.cv.glob, cv.tiles= avi.cv.tilesL)})); algs[[1]][['mean.prob']] = rowMeans(data.frame(algs)); algs })
lfi.cross.pred.glob <- sapply(lfinames.50presences, FUN=function(sp){ algs=list(sapply(names(lfi.models.glob[[sp]]), FUN=function(m){crossval.preds(lfi.models.glob[[sp]][[m]], dat= datLFI.train, sp=sp, pred=lfi.select07.cv.glob, cv.tiles= lfi.cv.tilesL)})); algs[[1]][['mean.prob']] = rowMeans(data.frame(algs)); algs })

# Assess cross-validated model performance:
avi.model.perf.ind <- sapply(avinames.50presences, FUN=function(sp){ list(sapply(names(avi.cross.pred.ind[[sp]]), FUN=function(m){calc.eval(dat= aviCH.forest.train, sp=sp, predictions= avi.cross.pred.ind[[sp]][[m]])})) })
lfi.model.perf.ind <- sapply(lfinames.50presences, FUN=function(sp){ list(sapply(names(lfi.cross.pred.ind[[sp]]), FUN=function(m){calc.eval(dat= datLFI.train, sp=sp, predictions= lfi.cross.pred.ind[[sp]][[m]])})) })
avi.model.perf.glob <- sapply(avinames.50presences, FUN=function(sp){ list(sapply(names(avi.cross.pred.glob[[sp]]), FUN=function(m){calc.eval(dat= aviCH.forest.train, sp=sp, predictions= avi.cross.pred.glob[[sp]][[m]])})) })
lfi.model.perf.glob <- sapply(lfinames.50presences, FUN=function(sp){ list(sapply(names(lfi.cross.pred.glob[[sp]]), FUN=function(m){calc.eval(dat= datLFI.train, sp=sp, predictions= lfi.cross.pred.glob[[sp]][[m]])})) })


#----------------------------------------------------------------------------------------------------------

# Fit MEMs to species richness data

# Function to fit MEMs using GLM, GAM, GBM and RF

fit.mem <- function(sp,pred,dat) {
	# GLM with linear and quadratic term:
	m.glm <- glm( as.formula(paste(sp,'~',paste(pred,paste0('+ I(',pred,'^2)'),collapse=' + '))), family='poisson', data=dat)
	
	# GAM from mgcv package with k=4
	m.gam <- gam( as.formula(paste(sp,'~',paste(paste0('s(',pred,',k=4)'),collapse=' + '))), family='poisson', data=dat)
	
	# Random Forest
	m.rf <- randomForest( x=dat[,pred], y=dat[,sp], ntree=1000, nodesize=20)
    
	# GBM from dismo package
	opt.LR <- TRUE
	LR <- 0.01
	while(opt.LR){
		m.brt <- try(gbm.step(data = dat, gbm.x = pred, gbm.y = sp, family = 'poisson', tree.complexity = 2, bag.fraction = 0.75, learning.rate = LR, verbose=F, plot.main=F))
        if (class(m.brt) == "try-error" | class(m.brt) == "NULL"){
            LR <- LR/2
        } else
        if(m.brt$gbm.call$best.trees<1000){
			LR <- LR/2
		} else 
		if(m.brt$gbm.call$best.trees>5000){
			LR <- LR*2
		} else { 
			opt.LR <- FALSE}
		}
	
	list(list(glm=m.glm, gam=m.gam, rf=m.rf, brt=m.brt))
}	


# Fit MEMs:
aviCH.forest.train$sr <- rowSums(aviCH.forest.train[,avinames.50presences])
avi.mem.glob <- fit.mem('sr',pred=avi.select07.cv.glob, dat= aviCH.forest.train)
avi.mem.ind <- fit.mem('sr',pred=avi.select07.cv.mem, dat= aviCH.forest.train)

datLFI.train$sr <- rowSums(datLFI.train[,lfinames.50presences])
lfi.mem.glob <- fit.mem('sr',pred=lfi.select07.cv.glob, dat= datLFI.train)
lfi.mem.ind <- fit.mem('sr',pred=lfi.select07.cv.mem, dat= datLFI.train)

# Make cross-validated predictions and ensembles:
avi.mem.cross.pred.glob <- data.frame(sapply(names(avi.mem.glob), FUN=function(m){crossval.preds(avi.mem.glob[[m]], dat= aviCH.forest.train, sp='sr', pred=avi.select07.cv.glob, cv.tiles= avi.cv.tilesL, gbm.dist='poisson')}))
avi.mem.cross.pred.glob$mean.prob <- rowMeans(avi.mem.cross.pred.glob)
avi.mem.cross.pred.ind <- data.frame(sapply(names(avi.mem.ind), FUN=function(m){crossval.preds(avi.mem.ind[[m]], dat= aviCH.forest.train, sp='sr', pred=avi.select07.cv.mem, cv.tiles= avi.cv.tilesL, gbm.dist='poisson')}))
avi.mem.cross.pred.ind$mean.prob <- rowMeans(avi.mem.cross.pred.ind)

lfi.mem.cross.pred.glob <- data.frame(sapply(names(lfi.mem.glob), FUN=function(m){crossval.preds(lfi.mem.glob[[m]], dat= datLFI.train, sp='sr', pred=lfi.select07.cv.glob, cv.tiles= lfi.cv.tilesL, gbm.dist='poisson')}))
lfi.mem.cross.pred.glob$mean.prob <- rowMeans(lfi.mem.cross.pred.glob)
lfi.mem.cross.pred.ind <- data.frame(sapply(names(lfi.mem.ind), FUN=function(m){crossval.preds(lfi.mem.ind[[m]], dat= datLFI.train, sp='sr', pred=lfi.select07.cv.mem, cv.tiles= lfi.cv.tilesL, gbm.dist='poisson')}))
lfi.mem.cross.pred.ind$mean.prob <- rowMeans(lfi.mem.cross.pred.ind)

# AUC
rcorr.cens(avi.mem.cross.pred.ind$mean.prob, aviCH.forest.train$sr)[1]
rcorr.cens(avi.mem.cross.pred.glob$mean.prob, aviCH.forest.train$sr)[1]
rcorr.cens(lfi.mem.cross.pred.ind$mean.prob, datLFI.train$sr)[1]
rcorr.cens(lfi.mem.cross.pred.glob$mean.prob, datLFI.train$sr)[1]

# correlation
cor(avi.mem.cross.pred.ind$mean.prob, aviCH.forest.train$sr)
cor(avi.mem.cross.pred.glob$mean.prob, aviCH.forest.train$sr)
cor(lfi.mem.cross.pred.ind$mean.prob, datLFI.train$sr)
cor(lfi.mem.cross.pred.glob$mean.prob, datLFI.train$sr)


# make MEM predictions on test data for later community predictions
avi.mem.glob.preds <- data.frame(sapply(names(avi.mem.glob), FUN=function(m){make.preds(avi.mem.glob[[m]], newdata= aviCH.forest.test, family='poisson')}))
avi.mem.glob.preds$mean.prob <- rowMeans(avi.mem.glob.preds)
avi.mem.ind.preds <- data.frame(sapply(names(avi.mem.ind), FUN=function(m){make.preds(avi.mem.ind[[m]], newdata= aviCH.forest.test, family='poisson')}))
avi.mem.ind.preds$mean.prob <- rowMeans(avi.mem.ind.preds)

lfi.mem.glob.preds <- data.frame(sapply(names(lfi.mem.glob), FUN=function(m){make.preds(lfi.mem.glob[[m]], newdata= datLFI.test, family='poisson')}))
lfi.mem.glob.preds$mean.prob <- rowMeans(lfi.mem.glob.preds)
lfi.mem.ind.preds <- data.frame(sapply(names(lfi.mem.ind), FUN=function(m){make.preds(lfi.mem.ind[[m]], newdata= datLFI.test, family='poisson')}))
lfi.mem.ind.preds$mean.prob <- rowMeans(lfi.mem.ind.preds)
