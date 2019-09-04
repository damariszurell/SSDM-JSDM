# Codes for: Zurell et al. "Testing species assemblage predictions from stacked and joint species distribution models" 

# Fitting JSDMs and evaluating performance based on spatial block cross-validation.

#----------------------------------------------------------------------------------------------------------

# load required packages
library(boral)

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

# Prepare input data for boral

avi.dat <- aviCH.forest.train
lfi.dat <- datLFI.train

# Add quadratic terms
avi.select07.cv.glob.sqr <- lfi.select07.cv.glob.sqr <- NULL
for (i in avi.select07.cv.glob) {
    i.new = paste0(i,'.sqr')
    avi.dat[,i.new] = avi.dat[,i]^2
    avi.select07.cv.glob.sqr = c(avi.select07.cv.glob.sqr, i.new)
}
for (i in lfi.select07.cv.glob) {
    i.new = paste0(i,'.sqr')
    lfi.dat[,i.new] = lfi.dat[,i]^2
    lfi.select07.cv.glob.sqr = c(lfi.select07.cv.glob.sqr, i.new)
}

#----------------------------------------------------------------------------------------------------------

# estimate JSDMs on all training data

setwd('jsdms')

# LVM
avi.n.burn <- 20000
avi.n.iter <- 50000
lfi.n.burn <- 50000
lfi.n.iter <- 100000
n.thin <- 50
n.lv <- 5

# Fit pure LVM without environmental predictors
avi.lvm <- boral(y = avi.dat[,avinames.50presences], family='binomial', num.lv=n.lv, row.eff='fixed', save.model=T, mcmc.control=list(n.burnin=avi.n.burn, n.iteration=avi.n.iter, n.thin=n.thin), model.name='lvm.avi.full')
avi.lvm.gew.pvals <- 2*pnorm(abs(unlist(avi.lvm$geweke.diag[[1]])), lower.tail = FALSE)
avi.lvm.conv <- p.adjust(avi.lvm.gew.pvals, method = "holm") 	# if resulting p values are below 0.05, the models might not have converged properly

lfi.lvm <- boral(y = lfi.dat[,lfinames.50presences], family='binomial', num.lv=n.lv, row.eff='fixed', save.model=T, mcmc.control=list(n.burnin=lfi.n.burn, n.iteration=lfi.n.iter, n.thin=n.thin), model.name='lvm.lfi.full')
lfi.lvm.gew.pvals <- 2*pnorm(abs(unlist(lfi.lvm$geweke.diag[[1]])), lower.tail = FALSE)
lfi.lvm.conv <- p.adjust(lfi.lvm.gew.pvals, method = "holm") 	# if resulting p values are below 0.05, the models might not have converged properly


# JSDM including environmental predictors
avi.jsdm <- boral(y = avi.dat[,avinames.50presences], X = avi.dat[,c(avi.select07.cv.glob,avi.select07.cv.glob.sqr)], family='binomial', num.lv=n.lv, save.model=T, mcmc.control=list(n.burnin=avi.n.burn, n.iteration=avi.n.iter, n.thin=n.thin), model.name='jsdm.avi.full')
avi.jsdm.gew.pvals <- 2*pnorm(abs(unlist(avi.jsdm$geweke.diag[[1]])), lower.tail = FALSE)
avi.jsdm.conv <- p.adjust(avi.jsdm.gew.pvals, method = "holm")     # if resulting p values are below 0.05, the models might not have converged properly

lfi.jsdm <- boral(y = lfi.dat[,lfinames.50presences], X = lfi.dat[,c(lfi.select07.cv.glob,lfi.select07.cv.glob.sqr)], family='binomial', num.lv=n.lv, save.model=T, mcmc.control=list(n.burnin=lfi.n.burn, n.iteration=lfi.n.iter, n.thin=n.thin), model.name='jsdm.lfi.full')
lfi.jsdm.gew.pvals <- 2*pnorm(abs(unlist(lfi.jsdm$geweke.diag[[1]])), lower.tail = FALSE)
lfi.jsdm.conv <- p.adjust(lfi.jsdm.gew.pvals, method = "holm")     # if resulting p values are below 0.05, the models might not have converged properly


#----------------------------------------------------------------------------------------------------------

# assess JSDM performance based on 5-fold spatial block cross-validation

# cross.validation avi
avi.cross.pred.jsdm <- data.frame(row = row.names(avi.dat), matrix(nrow=nrow(avi.dat),ncol=avinames.50presences)) 
names(avi.cross.pred.jsdm)[-1] <- avinames.50presences

for(i in seq_len(length(unique(avi.cv.tilesL))){
	cv.train <- avi.dat[avi.cv.tilesL!=i,]
	cv.test <- avi.dat[avi.cv.tilesL ==i,]
	
	# We update the model for the new training data
	avi.jsdm.cv <- boral(y = cv.train[,avinames.50presences], X = cv.train[,c(avi.select07.cv.glob,avi.select07.cv.glob.sqr)], family='binomial', num.lv=n.lv, save.model=T, mcmc.control=list(n.burnin=avi.n.burn, n.iteration=avi.n.iter, n.thin=n.thin), model.name=paste0('jsdm.avi.cv.',i))
	
	# Predictions on test data
	avi.cross.pred.jsdm[which(avi.cv.tilesL ==i),-1] = pnorm(predict(avi.jsdm.cv, newX = cv.test[,c(avi.select07.cv.glob,avi.select07.cv.glob.sqr)], predict.type='marginal', est='median', prob=0.95)[[1]])
}


# cross-validation lfi
lfi.cross.pred.jsdm <- data.frame(row = row.names(lfi.dat), matrix(nrow=nrow(lfi.dat),ncol=lfinames.50presences)) 
names(lfi.cross.pred.jsdm)[-1] <- lfinames.50presences

for(i in seq_len(length(unique(lfi.cv.tilesL))){
	cv.train <- lfi.dat[lfi.cv.tilesL!=i,]
	cv.test <- lfi.dat[lfi.cv.tilesL ==i,]
	
	# We update the model for the new training data
	lfi.jsdm.cv <- boral(y = cv.train[,lfinames.50presences], X = cv.train[,c(lfi.select07.cv.glob,lfi.select07.cv.glob.sqr)], family='binomial', num.lv=n.lv, save.model=T, mcmc.control=list(n.burnin=lfi.n.burn, n.iteration=lfi.n.iter, n.thin=n.thin), model.name=paste0('jsdm.lfi.cv.',i))
	
	# Predictions on test data
	lfi.cross.pred.jsdm[which(lfi.cv.tilesL ==i),-1] = pnorm(predict(lfi.jsdm.cv, newX = cv.test[,c(lfi.select07.cv.glob,lfi.select07.cv.glob.sqr)], predict.type='marginal', est='median', prob=0.95)[[1]])
}


# performance:
avi.model.perf.jsdm <- sapply(avinames.50presences, FUN=function(sp){ calc.eval(dat= aviCH.forest.train, sp=sp, predictions= avi.cross.pred.jsdm[[sp]])})
lfi.model.perf.jsdm <- sapply(lfinames.50presences, FUN=function(sp){ calc.eval(dat= datLFI.train, sp=sp, predictions= lfi.cross.pred.jsdm[[sp]])})
