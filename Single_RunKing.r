###################################################################################################
#                                                                                                 #
#   Single-RunKing: Bare-bones regression scan for genome-wide mixed model association analysis   #
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

Power_Mse <- function(dp,threshold,wt){
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
					sig1 <- c(sig1,sig0[j,1]) 	
					mse1=mse1+(sig0[j,2]-qtlinfo[sigp1,2])^2
					mse2=mse2+abs((sig0[j,2]-qtlinfo[sigp1,2])/qtlinfo[sigp1,2])
				}
				mse <- c(mse1/pow,mse2/pow)				
			}
		}else{
			mse=c(NA,NA)
		}    
		Pow <- c(Pow,pow)
		re_mse <- rbind(re_mse,mse)
	}
	PM <- list(rpow=Pow,rmse=re_mse)
	return(PM)
}

plink <- function(y,g){
	pbp <- matrix(nrow=ncol(g),ncol=3)
	for(i in 1:ncol(g)){
		fit <- fastLmPure(y=y,X=cbind(1,g[,i]))
		pbp[i,] <- c(i,fit$coefficients[2],fit$stderr[2])
	}
	pbp[,3] <- pbp[,2]/pbp[,3]
	p <- 2*pnorm(abs(pbp[,3]),lower.tail=FALSE)
	pbp <- cbind(pbp,p)
	return(pbp)
}

plink1 <- function(y,g){
	nobs <- nrow(g)
	pbp <- matrix(nrow=ncol(g),ncol=3)
	for(i in 1:ncol(g)){
		fit <- fastLmPure(y=y,X=cbind(gnewb0,g[,i]))
		pbp[i,] <- c(i,fit$coefficients[2],fit$stderr[2])
		resi <- y-(cbind(gnewb0,g[,i]))%*%fit$coefficients
		rv <- sum(resi^2)/(nobs-2)
		pbp[i,3] <- pbp[i,2]/pbp[i,3]*sqrt(rv)/sqrt(ve0)		
	}
	p <- 2*pnorm(abs(pbp[,3]),lower.tail=FALSE)
	pbp <- cbind(pbp,p)
	return(pbp)
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
library(snow)
library(RcppArmadillo)
library(MASS)
dir <- getwd()
setwd(dir)
g<-fread("/data/newdata/Ames_GDD_DTS_GD.csv",nrows = -1,sep=",")
g <- t(g)
scp <- read.table("/data/newdata/maize_SCP.txt",sep="\t",header=T)
colnames(g) <- scp[,1]
Var <- apply(g,2,var)
g <- g[,Var>0.05]
gid <- sample(1:ncol(g),300000)
gid <- sort(gid)
g <- g[,gid]
nobs<-nrow(g)
nmar<-ncol(g)
snpnames <- colnames(g)
chrnames <- (scp[,2])[scp[,1]%in%snpnames]
nan_chr <- c(which(chrnames==6),which(chrnames==8))
snpid <- c(1:nmar)[-nan_chr]
cat('position_6_8','\n',length(nan_chr),'\n')
mk <- mean(g)
h2   <- c(0.6)
nqtl <- c(100)
threshold<-10^(-seq(10,2,-0.2))
wt<-0
delt <- 20

############################################################
#                     Spectral decompose                   #
############################################################
A <- scale(g)
B <- t(A)
cl <- makeSOCKcluster(rep("localhost",20))
G <- parMM(cl, A, B)
stopCluster(cl)
G <- G/nmar
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

As <- g1[,sample(1:nmar,10000)]
G2 <- As%*%t(As)
G2 <- G2/ncol(As)
kpr <- prcomp(G2)
kpc <- predict(kpr)[,c(1:3)]

##################################################################################
#                                  Execute program                               #
##################################################################################

############################################################
#         Simulating phenotypes for Continous trait        #
############################################################
	qtlinfo <- SimQTL(snpid,nqtl)
	TQP <- qtlinfo[,1]
	tbv <- g1[,TQP]%*%qtlinfo[,2]
	VarG <- var(tbv)
	Ve <- 5.0
	Va <- Ve*h2/(1-h2)
	qtlinfo[,2] <- qtlinfo[,2]%*%sqrt(Va/VarG)
	tbv1 <- tbv%*%sqrt(Va/VarG)
	err <- rnorm(nobs,mean=0,sd=sqrt(Ve))
	y <- tbv1+err
	y <- y-min(y)+0.1
	kk <- lm(y~kpc)
	resi_pc <- summary(kk)$residual
	ynew <- t(ug)%*%y
	
############################################################
#             GBLUP(by spectral transformation)            #
############################################################		
	fva0 <- GBLUP(ynew,as.matrix(gnewb0))
	hh0 <- fva0$h2
	ve0 <- fva0$ve

############################################################
#              Regresison One by one(PLINK)                #
############################################################	
	dp <- plink(resi_pc,g1)
	pv_Oney <- dp[,4]
	dp <- posi_Choice3(dp[,-3],delt)
	dp_Oney <- posi_Choice3(dp,delt)
	PPM_Oney <- Power_Mse(dp_Oney,threshold,wt)	
	Pva_NQA_Oney <- sort(as.matrix(pv_Oney[nan_chr]))
	Pva_QA_Oney <- sort(pv_Oney[snpid])

############################################################
#                    Single-RunKing                        #
############################################################
	hpos <- which(hh==hh0)
	ynew1 <- V050*ynew[,1]
	dp <- matrix(nrow=nmar,ncol=4)
	for(j in 1:nmar){
		fva <- fast_lmm(ynew1,cbind(gnewb0,gnew[,j]),1)
		dp[j,] <- c(j,fva)
	}
	pv_SR <- dp[,3]
	dp_SR <- posi_Choice3(dp[,-4],delt)
	PPM_SR <- Power_Mse(dp_SR,threshold,wt)	
	Pva_NQA_SR <- sort(as.matrix(pv_SR[nan_chr]))
	Pva_QA_SR <- sort(pv_SR[snpid])

############################################################
#                         EMMAX                            #
############################################################
	sp <- hh0*1000+1
	ynew1 <- V050[,sp]*ynew
	gnew1 <- V050[,sp]*gnew
	dp <- plink1(ynew1,gnew1)
	pv_EMMAX <- dp[,4]
	dp_EMMAX <- posi_Choice3(dp[,-3],delt)
	PPM_EMMAX <- Power_Mse(dp_EMMAX,threshold,wt)
	Pva_NQA_EMMAX <- sort(as.matrix(pv_EMMAX[nan_chr]))
	Pva_QA_EMMAX <- sort(pv_EMMAX[snpid])

############################################################
#                         GRAMMAR                          #
############################################################
	b0 <- as.numeric(fva0$eff)
	detv <- hh0/(1-hh0)
	ebv <- G%*%(detv*solve(G*detv+diag(nobs)))%*%(y-b0)
	dp <- plink(y-ebv,g1)
	pv_GRAMMAR <- dp[,4]
	dp_GRAMMAR <- posi_Choice3(dp[,-3],delt)
	PPM_GRAMMAR <- Power_Mse(dp_GRAMMAR,threshold,wt)	
	Pva_NQA_GRAMMAR <- sort(as.matrix(pv_GRAMMAR[nan_chr]))
	Pva_QA_GRAMMAR <- sort(pv_GRAMMAR[snpid])

