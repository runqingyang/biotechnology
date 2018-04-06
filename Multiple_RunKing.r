###################################################################################################
#                                                                                                 #
#           Multi-RunKing: Genome-wide multi-locus as-sociation study based on R/glmnet           #
#                                                                                                 #
###################################################################################################

##################################################################################
#                                    Subroutine                                  #
##################################################################################

meansub<-function(x){
  mm<-mean(x,na.rm=TRUE)
  x[is.na(x)]<-mm
  x
}

SimQTL <- function(snpid,nqtl){
    qid <- sample(snpid,nqtl)
    eff <- rgamma(nqtl, shape=0.4, scale=1.66)
    sig <- sample(c(1,-1),nqtl,replace=T)
    eff <- eff*sig
	cbind(qid, eff)
}

posi_Choice3 <- function(dp,ws) { 
	if (nrow(dp) <= 2) return(dp)
	if (nrow(dp) >  2) {
		Pos <- c()	
		for(j in 1:nqtl){
			sigpos<-abs(dp[,1]-qtlinfo[j,1])
			sig <- which(sigpos<=ws)
			if(length(sig)!=0){
				b1 <- dp[sig,]
				if(length(b1)==3){b1 <- matrix(b1,1)}
				sigpos <- which.min(b1[,3])
				Pos <- rbind(Pos,b1[sigpos,])
			}
		}
	}
	Pos <- unique(Pos)
	return(Pos)
}

Power_Proximity_Mse <- function(dp,threshold,wt){
	Pow=c()
	re_nn <- c()
	re_mse <- c()
	for(i in 1:length(threshold)){
		sig <- which(dp[,3]<=threshold[i])
		if(length(sig)==1){
			sig0 <- matrix(dp[sig,],1)
		}else{
			sig0 <- dp[sig,]
		}
		sl <- length(sig0[,1])
		pow=0
		mse1=0
		mse2=0
		sig1 <- c()
		if(sl!=0){
			for(j in 1:sl){		
				a <- abs(sig0[j,1]-qtlinfo[,1])
				sigpos0 <- min(a)
				sigp1 <- which.min(a)
				if(sigpos0<=wt){ 
					pow=pow+1
					sig1 <- c(sig1,sig0[j,1]) #if sl==0, similarity will report an error	
					mse1=mse1+(sig0[j,2]-qtlinfo[sigp1,2])^2
					mse2=mse2+abs((sig0[j,2]-qtlinfo[sigp1,2])/qtlinfo[sigp1,2])
				}
				mse <- c(mse1/pow,mse2/pow)				
			}
		}else{
			mse=c(NA,NA)
		}
		n <- sig1%in%qtlinfo[order(-abs(qtlinfo[,2])),1][1:length(sig1)]
		nn <- length(which(n==TRUE))/length(sig1)
		Pow <- c(Pow,pow)
		re_mse <- rbind(re_mse,mse)
		re_nn <- c(re_nn,nn)
	}
	PPM <- list(rpow=Pow,rmse=re_mse,rnn=re_nn)
	return(PPM)
}

GBLUP <- function(ynew,gnew){
	nobs <- nrow(gnew)
	pre_loglike <- 10^10
	ynew <- V050*ynew[,1]
	for(i in 1:length(hh)){
		gnew1 <- V050[,i]*gnew
		fit=fastLmPure(y = ynew[,i], X = gnew1)
		eff <- fit$coefficients		
		resi <- ynew[,i]-gnew1%*%as.matrix(eff)
		ve <- sum(resi^2)/(nobs-ncol(gnew)-1)
		loglike <- logV0[i]+nobs*log(ve)
		if(loglike > pre_loglike){
			break
		}else{
			pre_eff <- eff 
			h2 <- hh[i] 
			pre_ve <- ve 
			pre_loglike <- loglike
			}
	}
	fva <- list(eff=pre_eff,h2=h2,ve=pre_ve,loglike=pre_loglike)
	return(c(fva))
}

fastlm <- function(i,ynew1,gnew1){
	fit <- fastLmPure(y=ynew1[,i],X=gnew1)
	eff <- fit$coefficients
	resi <- ynew1[,i]-(gnew1[,1]*eff[1]+gnew1[,2]*eff[2])
	ve <- sum(resi^2)/(nobs-2)
	ster <- fit$stderr[2]
	t <- eff[2]/ster*sqrt(ve)/sqrt(ve0)
	p <- 2*pnorm(abs(t),lower.tail=FALSE)
	loglike0 <- logV0[i]+nobs*log(ve)
	fastlm <- c(eff[2],p,loglike0)
}

fast_lmm <- function(ynew1,gnew,threshold0){
	gnew1 <- V050[,hpos]*gnew
	f10 <- fastlm(hpos,ynew1,gnew1)
	if(f10[2] < threshold0){
		pre_eff <- f10[1]
		pre_p <- f10[2]
		pre_loglike <- f10[3]
		for(j in 1:(hpos-1)){
			i <- hpos-j
			gnew1 <- V050[,i]*gnew
			f10 <- fastlm(i,ynew1,gnew1)
			loglike <- f10[3]
			if(loglike > pre_loglike) break
			if(i==1) cat("No Heritability!","\n")
			pre_loglike <- loglike
			pre_eff <- f10[1]
			pre_p <- f10[2]			
		}
		hhk <- hh[i+1]
		f10 <- c(pre_eff,pre_p,pre_loglike)
			
		if(hhk==hh0){
			pre_eff <- f10[1]
			pre_p <- f10[2]
			pre_loglike <- f10[3]
			for(j in 1:(length(hh)-hpos)){
				i <- hpos+j
				gnew1 <- V050[,i]*gnew
				f10 <- fastlm(i,ynew1,gnew1)
				loglike <- f10[3]
				if(loglike > pre_loglike) break
				if(i==length(hh)) cat("No Solution!","\n")
				pre_loglike <- loglike
				pre_eff <- f10[1]
				pre_p <- f10[2]					
			}
			f10 <- c(pre_eff,pre_p,pre_loglike)
		}	
    }
	return(f10)	
}

##################################################################################
#                            Load packages and data                              #
##################################################################################

library(data.table)
library(glmnet)
library(RcppArmadillo)
library(Rcpp)
library(snow)
library(survival)
dir <- getwd()
setwd(dir)
g <- fread("/data1/mydata/g_human_300k.txt",nrows = -1,sep=" ",stringsAsFactors=FALSE)
g <- as.matrix(g)
nobs <- nrow(g)
nmar <- ncol(g)
trait <- "Continous" #"Binary";"Survival"
h2   <- c(0.6)
nqtl <- c(100)
threshold<-10^(-seq(10,2,-0.2))
wt<-0
delt <- 10

############################################################
#                     Spectral decompose                   #
############################################################
A <- scale(g)
B <- t(A)
cl <- makeSOCKcluster(rep("localhost",20))
G <- parMM(cl, A, B)
stopCluster(cl)
G <- G/nmar
diag(G) <- diag(G) + 0.001
g1 <- A
eig <- eigen(G)
sg <- eig$values
ug <- eig$vectors
step <- 1
hh <- seq(0,0.999,0.001*step)
hh <- round(hh,3)
sg <- matrix(rep(sg,length(hh)),nobs)
V <- t(t(sg)*hh/(1-hh)) + 1
logV0 <- apply(log(V),2,sum)
V050 <- 1/sqrt(V)
gnewb0 <- t(ug)%*%as.matrix(rep(1,nobs))
cl <- makeSOCKcluster(rep("localhost",20))
gnew <- parMM(cl,t(ug),g1)
stopCluster(cl)

##################################################################################
#                                  Execute program                               #
##################################################################################


############################################################
#         Simulating phenotypes for Continous trait        #
############################################################	
	if(trait=="Continous"){
		qtlinfo <- SimQTL(c(1:nmar),nqtl)
		TQP <- qtlinfo[,1]
		tbv <- g1[,TQP]%*%qtlinfo[,2]
		VarG <- var(tbv)
		Ve=5.0
		Va=Ve*h2/(1-h2)
		tbv1=tbv%*%sqrt(Va/VarG)
		qtlinfo[,2]=qtlinfo[,2]%*%sqrt(Va/VarG)
		err <- rnorm(nobs,mean=0,sd=sqrt(Ve))
		y   <- tbv1+err
		y   <- y-min(y)+0.1
		ynew <- t(ug)%*%y
	}
	
############################################################
#          Simulating phenotypes for Binary trait          #
############################################################
	if(trait=="Binary"){
		qtlinfo <- SimQTL(1:nmar,nqtl) 
		TQP <- qtlinfo[,1]
		tbv <- g1[,TQP]%*%qtlinfo[,2]
		incidence <- 0.5
		vart <- var(tbv)
		ve <- 5.0
		c0 <- pi/sqrt(3)
		Va <- ve*h2/(1-h2)
		qtlinfo[,2] <- qtlinfo[,2]%*%sqrt(Va/vart)  
		fx1 <- tbv%*%sqrt(Va/vart) 
		thr <- log(incidence/(1-incidence)) #logit model的阈值
		fx <- thr+fx1-mean(fx1)                
		px <- exp(c0*fx)
		px <- px/(1+px)
		y <- rbinom(length(px),1,px)
		ynew <- t(ug)%*%y	
	}
	
############################################################
#         Simulating phenotypes for Survival traits        #
############################################################	
	if(trait=="Survival"){	
		qtlinfo <- SimQTL(1:nmar,nqtl) 
		TQP <- qtlinfo[,1]
		tbv <- g[,TQP]%*%qtlinfo[,2]
		vart <- var(tbv)
		Ve <- 5.0
		Va <- Ve*h2/(1-h2)
		qtlinfo[,2] <- qtlinfo[,2]%*%sqrt(Va/vart)  
		fx1 <- tbv%*%sqrt(Va/vart)
		lam <- 0.0000001
		v <- 6
		u <- runif(nobs) 
		y <- (-log(u)/(lam*exp(fx1)))^(1/v)
		status <- rep(1,nobs)
		ynew <- t(ug)%*%log(y)
	}
	
############################################################
#             GBLUP(by spectral transformation)            #
############################################################					
	fva0 <- GBLUP(ynew,as.matrix(gnewb0))
	hh0 <- fva0$h2
	ve0 <- fva0$ve
	
############################################################
#                           FaSTLMM                        #
############################################################	
	hpos <- which(hh==hh0)
	ynew1 <- V050*ynew[,1]
	dp <- matrix(nrow=nmar,ncol=4)
	for(j in 1:nmar){
		fv <- fast_lmm(ynew1,cbind(gnewb0,gnew[,j]),0.05)
		dp[j,] <- c(j,fv[c(3,1,2)])
	}
	dp_FaSTLMM <- posi_Choice3(dp[,-2],delt)
	PPM_FaSTLMM <- Power_Proximity_Mse(dp_FaSTLMM,threshold,wt)	
	
############################################################
#                        Multi_RunKing                     #
############################################################	
	require(doMC)
	registerDoMC(32)
	if(trait=="Continous") m1 <- cv.glmnet(g1,y,parallel=TRUE)
	if(trait=="Binary") m1 <- cv.glmnet(g1,y,family='binomial',parallel=TRUE)
	if(trait=="Survival"){
		surv<-Surv(y,as.matrix(status))
		m1 <- cv.glmnet(g,surv,family='cox',maxit = 1000,parallel=TRUE)
	}
	b1 <- summary(coef(m1))[,c(1,3)]
	if(trait=="Continous"|trait=="Binary"){
		d  <- b1[-1,1] -1
		bb <- b1[-1,2]
	}
	if(trait=="Survival"){
		d  <- b1[,1] 
		bb <- b1[,2]
	}
	if(length(d)==1) bb <- as.matrix(bb)
	dp <- cbind(d,d,bb)
	dp[,3] <- 1-abs(dp[,3])
	d3 <- posi_Choice3(dp,delt)
	dd <- d3[1]
	if(length(d3)>2){dd <- d3[,1]}
	xn1<-as.matrix(g[,dd])
	colnames(xn1)<-paste("V",1:length(dd),sep="")
	if(trait=="Continous") fit0<-lm(y~xn1)
	if(trait=="Binary") fit0<-glm(y~xn1,family = 'binomial')
	if(trait=="Survival") fit0<-coxph(surv~xn1)
	coeffi1<-summary(fit0)$coefficients
	names(dd)<-paste("xn1V",1:length(dd),sep="")
	index<-dd[intersect(names(dd),rownames(coeffi1))]
	d_Las<-as.matrix(index)
	if(trait=="Continous"|trait=="Binary"){
		eff_pv <- coeffi1[-1,c(1,4)]
		dp_MR <- cbind(d_Las,eff_pv)
	}
	if(trait=="Survival"){
		eff_pv <- coeffi1[,c(1,5)]
		dp_MR <- cbind(d_Las,eff_pv)
	}
	PPM_MR <- Power_Proximity_Mse(dp_MR,threshold,wt)	
	
############################################################
#                          LMM_LASSO                       #
############################################################
	if(trait=="Continous"){	
		dd <- as.numeric(d_Las)
		gnewk <- cbind(gnewb0,gnew[,dd])*V050[,hpos]	
		fit <- fastLmPure(y = ynew1[,hpos], X =gnewk)
		eff <- fit$coefficients
		resi <- ynew1[,hpos]-gnewk%*%eff
		ve <- sum(resi^2)/(nobs-length(dd)-1)
		ster <- fit$stderr[-1]
		t <- eff[-1]/ster*sqrt(ve)/sqrt(ve0)
		p <- 2*pnorm(abs(t),lower.tail=FALSE)
		dp_LMMLasso <- cbind(dd,eff[-1],p)
		PPM_LMMLasso <- Power_Proximity_Mse(dp_LMMLasso,threshold,wt)	
	}	

	
	
	
	
	