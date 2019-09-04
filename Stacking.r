# Codes for: Zurell et al. "Testing species assemblage predictions from stacked and joint species distribution models" 

# Fitting SDMs and evaluating performance based on spatial block cross-validation.

#----------------------------------------------------------------------------------------------------------

# load required packages
require(PresenceAbsence)
require(ecospat)

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

# Functions

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


# probability ranking rule - adapted from the ecospat package
SESAM.prr <- function (proba, sr) {
    projSR <- round(round(as.vector(sr[[1]])))
    new.prob.prr <- proba
    dataSSDM_p <- proba
    for (i in 1:nrow(proba)) {
        print(paste("test.prr, processing row ", i, sep = ""))
        SR <- projSR[i]
        if (SR > 0) {
            predcom <- dataSSDM_p[i, ]
            predcom_p <- dataSSDM_p[i, ]
            com <- order(predcom_p, decreasing = TRUE)
            pres <- com[1:SR]
            predcom[, pres] <- 1
            predcom[, -pres] <- 0
        }
        else {
            predcom[, ] <- 0
        }
        new.prob.prr[i, ] <- predcom
    }
    new.prob.prr
}

# Calabrese correction: based on Calabrese et al. (2014) Methods in Ecology and Evolution 23: 99-112.
nLL.Calabrese <- function(par,sr,probs) {
	require(poibin)
	logit = function(x) {x=ifelse(x<0.0001,0.0001,ifelse(x>0.9999,.9999,x));   ;log(x/(1 - x))}
	invlogit = function(x) {exp(x)/(1+exp(x))}

	
	bysite <- function(j) {
		logit.probs <- logit(as.numeric(probs[j,]))
		corr.probs <- invlogit( logit.probs + par[1]*sr[j] + par[2] )
		dp <- dpoibin(sr[j],as.numeric(corr.probs))
		log(ifelse(dp<.0001,.0001,dp))
	}
	- sum(sapply(seq_len(length(sr)),bysite)) 	# optim will perform minimization but we aim for maximum likelihood and thus invert
	}


logit <- function(x) {x=ifelse(x<0.0001,0.0001,ifelse(x>0.9999,.9999,x));   ;log(x/(1 - x))}
invlogit <- function(x) {exp(x)/(1+exp(x))}


# function for stacking
stacking <- function(preds, thresh, mem.preds, obs.sr.train, preds.train){

	# stack probabilities
	prob.df <- data.frame(preds)
	prob.stack <- rowSums(preds)
	
	# stack binaries
	bin.df <- data.frame(sapply(names(data.frame(preds)),FUN=function(sp){ifelse(data.frame(preds)[sp]<thresh[sp],0,1)}))
	bin.stack <- rowSums(bin.df)
	
	# estimated correction parameters based on Calabrese correction
	adj.par.nll <- optim(par=c(0,0), fn=nLL.Calabrese, sr= obs.sr.train, probs= preds.train) 
	
	# correct probabilities using probability stacks on test data
	prob.corr.probsum.df <- data.frame( apply(prob.df,2,FUN=function(x){invlogit(logit(x)+adj.par.nll$par[1]*prob.stack+adj.par.nll$par[2])}))
	prob.corr.probsum.stack <- rowSums(prob.corr.probsum.df)
	
	# correct probabilities using MEMs on test data
	prob.corr.mem.df <- data.frame( apply(prob.df,2,FUN=function(x){invlogit(logit(x)+adj.par.nll$par[1]*mem.preds+adj.par.nll$par[2])}))
	prob.corr.mem.stack <- rowSums(prob.corr.mem.df)
	
	# probability ranking rule, using probability sums as constraint
	prr.probsum.df <- SESAM.prr(prob.df, data.frame(prob.stack) )
	prr.probsum.stack <- rowSums(prr.probsum.df)
	
	# probability ranking rule, using MEMs as constraint
	prr.mem.df <- SESAM.prr(prob.df, data.frame(mem.preds) )
	prr.mem.stack <- rowSums(prr.mem.df)
	
	# probability ranking rule, using sums of corrected probabilities as constraint
	prr.corr.probsum.df <- SESAM.prr(prob.corr.probsum.df, data.frame(prob.corr.probsum.stack) )
	prr.corr.probsum.stack <- rowSums(prr.corr.probsum.df)
	
	# probability ranking rule, using sums of corrected probabilities as constraint (but probabilities corrected based on MEM prediction) 
	prr.corr.mem.df <- SESAM.prr(prob.corr.mem.df, data.frame(prob.corr.mem.stack) )
	prr.corr.mem.stack <- rowSums(prr.corr.mem.df)

	
	return(list(spp=list(prof.df=prob.df, bin.df=bin.df, prob.corr.probsum.df= prob.corr.probsum.df, prob.corr.mem.df= prob.corr.mem.df, prr.corr.probsum.df= prr.corr.probsum.df, prr.corr.mem.df= prr.corr.mem.df, prr.probsum.df= prr.probsum.df, prr.mem.df= prr.mem.df), stack=list(prob.stack=prob.stack, bin.stack=bin.stack, prob.corr.probsum.stack= prob.corr.probsum.stack, prob.corr.mem.stack= prob.corr.mem.stack, prr.corr.probsum.stack= prr.corr.probsum.stack, prr.corr.mem.stack= prr.corr.mem.stack, prr.probsum.stack= prr.probsum.stack, prr.mem.stack= prr.mem.stack)))
}


# evaluate community predictions
comm.eval <- function(thestack, test.obs, n.tir.max=100){
	# check for integer predictions
	is.prob <- sapply(seq_len(length(thestack$stack)),FUN=function(x){sum((thestack$stack[[x]]*1000)%%1000)>0})	
	
	evals.lst  <- lapply(seq_len(length(thestack$spp)),function(i){ecospat.CommunityEval(test.obs, thestack$spp[[i]], is.prob[i], ntir=ifelse(is.prob[i],n.tir.max,1) )})
	names(evals.lst) <- names(thestack$spp)
	return(evals.lst)
}


#----------------------------------------------------------------------------------------------------------

# load cross-validated model performance:
avi.model.perf.glob
lfi.model.perf.glob
avi.model.perf.ind
lfi.model.perf.ind
avi.model.perf.jsdm 
lfi.model.perf.jsdm 

# load predictions on independent test data:
avi.glob.preds 
avi.ind.preds 
lfi.glob.preds 
lfi.ind.preds 
lfi.jsdm.preds
avi.jsdm.preds


# make assemblage and richness predictions
avi.jsdm.stacks <- stacking(preds=avi.jsdm.preds, thresh=apply(avi.model.perf.jsdm,2,FUN=function(x){x$thresh}), mem.preds=avi.mem.ind.preds$mean.prob, obs.sr.train=rowSums(aviCH.forest.train[,avinames.50presences]), preds.train=avi.cross.pred.jsdm[avinames.50presences])
avi.glob.ens.stacks <- stacking(preds=sapply(names(avi.glob.preds),FUN=function(sp){avi.glob.preds[[sp]]$prob.means}), thresh=sapply(names(avi.model.perf.glob),FUN=function(sp){avi.model.perf.glob[[sp]][,'mean.prob']$thresh}), mem.preds=avi.mem.ind.preds$mean.prob, obs.sr.train=rowSums(aviCH.forest.train[,avinames.50presences]), preds.train=sapply(names(avi.cross.pred.glob),FUN=function(sp){avi.cross.pred.glob[[sp]]$mean.prob}))
avi.glob.glm.stacks <- stacking(preds=sapply(names(avi.glob.preds),FUN=function(sp){avi.glob.preds[[sp]]$glm}), thresh=sapply(names(avi.model.perf.glob),FUN=function(sp){avi.model.perf.glob[[sp]][,'glm']$thresh}), mem.preds=avi.mem.ind.preds$mean.prob, obs.sr.train=rowSums(aviCH.forest.train[,avinames.50presences]), preds.train=sapply(names(avi.cross.pred.glob),FUN=function(sp){avi.cross.pred.glob[[sp]]$glm}))
avi.jsdm.stacks.glob.mem <- stacking(preds=avi.jsdm.preds, thresh=apply(avi.model.perf.jsdm,2,FUN=function(x){x$thresh}), mem.preds=avi.mem.glob.preds$mean.prob, obs.sr.train=rowSums(aviCH.forest.train[,avinames.50presences]), preds.train=avi.cross.pred.jsdm[avinames.50presences])
avi.glob.ens.stacks.glob.mem <- stacking(preds=sapply(names(avi.glob.preds),FUN=function(sp){avi.glob.preds[[sp]]$prob.means}), thresh=sapply(names(avi.model.perf.glob),FUN=function(sp){avi.model.perf.glob[[sp]][,'mean.prob']$thresh}), mem.preds=avi.mem.glob.preds$mean.prob, obs.sr.train=rowSums(aviCH.forest.train[,avinames.50presences]), preds.train=sapply(names(avi.cross.pred.glob),FUN=function(sp){avi.cross.pred.glob[[sp]]$mean.prob}))
avi.glob.glm.stacks.glob.mem <- stacking(preds=sapply(names(avi.glob.preds),FUN=function(sp){avi.glob.preds[[sp]]$glm}), thresh=sapply(names(avi.model.perf.glob),FUN=function(sp){avi.model.perf.glob[[sp]][,'glm']$thresh}), mem.preds=avi.mem.glob.preds$mean.prob, obs.sr.train=rowSums(aviCH.forest.train[,avinames.50presences]), preds.train=sapply(names(avi.cross.pred.glob),FUN=function(sp){avi.cross.pred.glob[[sp]]$glm}))
avi.ind.ens.stacks <- stacking(preds=sapply(names(avi.ind.preds),FUN=function(sp){avi.ind.preds[[sp]]$prob.means}), thresh=sapply(names(avi.model.perf.ind),FUN=function(sp){avi.model.perf.ind[[sp]][,'mean.prob']$thresh}), mem.preds=avi.mem.ind.preds$mean.prob, obs.sr.train=rowSums(aviCH.forest.train[,avinames.50presences]), preds.train=sapply(names(avi.cross.pred.ind),FUN=function(sp){avi.cross.pred.ind[[sp]]$mean.prob}))
avi.ind.glm.stacks <- stacking(preds=sapply(names(avi.ind.preds),FUN=function(sp){avi.ind.preds[[sp]]$glm}), thresh=sapply(names(avi.model.perf.ind),FUN=function(sp){avi.model.perf.ind[[sp]][,'glm']$thresh}), mem.preds=avi.mem.ind.preds$mean.prob, obs.sr.train=rowSums(aviCH.forest.train[,avinames.50presences]), preds.train=sapply(names(avi.cross.pred.ind),FUN=function(sp){avi.cross.pred.ind[[sp]]$glm}))

lfi.jsdm.stacks <- stacking(preds=lfi.jsdm.preds, thresh=apply(lfi.model.perf.jsdm,2,FUN=function(x){x$thresh}), mem.preds=lfi.mem.ind.preds$mean.prob, obs.sr.train=rowSums(datLFI.train[,lfinames.50presences]), preds.train=lfi.cross.pred.jsdm[lfinames.50presences])
lfi.glob.ens.stacks <- stacking(preds=sapply(names(lfi.glob.preds),FUN=function(sp){lfi.glob.preds[[sp]]$prob.means}), thresh=sapply(names(lfi.model.perf.glob),FUN=function(sp){lfi.model.perf.glob[[sp]][,'mean.prob']$thresh}), mem.preds=lfi.mem.ind.preds$mean.prob, obs.sr.train=rowSums(datLFI.train[,lfinames.50presences]), preds.train=sapply(names(lfi.cross.pred.glob),FUN=function(sp){lfi.cross.pred.glob[[sp]]$mean.prob}))
lfi.glob.glm.stacks <- stacking(preds=sapply(names(lfi.glob.preds),FUN=function(sp){lfi.glob.preds[[sp]]$glm}), thresh=sapply(names(lfi.model.perf.glob),FUN=function(sp){lfi.model.perf.glob[[sp]][,'glm']$thresh}), mem.preds=lfi.mem.ind.preds$mean.prob, obs.sr.train=rowSums(datLFI.train[,lfinames.50presences]), preds.train=sapply(names(lfi.cross.pred.glob),FUN=function(sp){lfi.cross.pred.glob[[sp]]$glm}))
lfi.jsdm.stacks.glob.mem <- stacking(preds=lfi.jsdm.preds, thresh=apply(lfi.model.perf.jsdm,2,FUN=function(x){x$thresh}), mem.preds=lfi.mem.glob.preds$mean.prob, obs.sr.train=rowSums(datLFI.train[,lfinames.50presences]), preds.train=lfi.cross.pred.jsdm[lfinames.50presences])
lfi.glob.ens.stacks.glob.mem <- stacking(preds=sapply(names(lfi.glob.preds),FUN=function(sp){lfi.glob.preds[[sp]]$prob.means}), thresh=sapply(names(lfi.model.perf.glob),FUN=function(sp){lfi.model.perf.glob[[sp]][,'mean.prob']$thresh}), mem.preds=lfi.mem.glob.preds$mean.prob, obs.sr.train=rowSums(datLFI.train[,lfinames.50presences]), preds.train=sapply(names(lfi.cross.pred.glob),FUN=function(sp){lfi.cross.pred.glob[[sp]]$mean.prob}))
lfi.glob.glm.stacks.glob.mem <- stacking(preds=sapply(names(lfi.glob.preds),FUN=function(sp){lfi.glob.preds[[sp]]$glm}), thresh=sapply(names(lfi.model.perf.glob),FUN=function(sp){lfi.model.perf.glob[[sp]][,'glm']$thresh}), mem.preds=lfi.mem.glob.preds$mean.prob, obs.sr.train=rowSums(datLFI.train[,lfinames.50presences]), preds.train=sapply(names(lfi.cross.pred.glob),FUN=function(sp){lfi.cross.pred.glob[[sp]]$glm}))
lfi.ind.ens.stacks <- stacking(preds=sapply(names(lfi.ind.preds),FUN=function(sp){lfi.ind.preds[[sp]]$prob.means}), thresh=sapply(names(lfi.model.perf.ind),FUN=function(sp){lfi.model.perf.ind[[sp]][,'mean.prob']$thresh}), mem.preds=lfi.mem.ind.preds$mean.prob, obs.sr.train=rowSums(datLFI.train[,lfinames.50presences]), preds.train=sapply(names(lfi.cross.pred.ind),FUN=function(sp){lfi.cross.pred.ind[[sp]]$mean.prob}))
lfi.ind.glm.stacks <- stacking(preds=sapply(names(lfi.ind.preds),FUN=function(sp){lfi.ind.preds[[sp]]$glm}), thresh=sapply(names(lfi.model.perf.ind),FUN=function(sp){lfi.model.perf.ind[[sp]][,'glm']$thresh}), mem.preds=lfi.mem.ind.preds$mean.prob, obs.sr.train=rowSums(datLFI.train[,lfinames.50presences]), preds.train=sapply(names(lfi.cross.pred.ind),FUN=function(sp){lfi.cross.pred.ind[[sp]]$glm}))


#----------------------------------------------------------------------------------------------------------
#community evaluation using ecospat.CommunityEval

mystacks <- grep('\\.stacks',ls(),value=T)
mystacks.avi <- grep('avi',mystacks,value=T)
mystacks.lfi <- grep('lfi',mystacks,value=T)

avi.comm.eval <- sapply(mystacks.avi,FUN=function(s){list(comm.eval(get(s),aviCH.forest.test[,avinames.50presences]))})
lfi.comm.eval <- sapply(mystacks.lfi,FUN=function(s){list(comm.eval(get(s), datLFI.test[,lfinames.50presences]))})


