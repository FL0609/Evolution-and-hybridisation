#library(scales)
############  Obtain a file name from a list of parameters
getfile.norep <- function(ff) {
	return(paste0('res_N',ff[1],'_L',ff[2],'_U',formatC(ff[3]),'_k',ff[4],'_a',formatC(ff[5]),'_n',ff[6],'_F',ff[7],'_q',ff[8],'_MVN',ff[9],'_diff',ff[10],'_rep'))
}

############ Obtain squared distance from the optimum for a single genotype
getS <- function(g1,g2,zm,doms,zstart,opt=c(1,0)) 
{ 
	hom <- which( (g1+g2) == 2)
	het <- which( g1 != g2)
	
	m <- rbind(zm[hom,],zm[het,]*doms[het,])
	if(!is.null(dim(m)))
	{	
		z <- zstart+apply(m,MARGIN=2,FUN=sum)
		return(sum(z^2)) # Sum across traits
	}
	z <- (zstart+sum(m))-opt
	return(sum(z^2))
}	
############ The fitness function ############
w <- function(S,k)
{
	return(exp(-0.5*S^(k/2)))	
}
############ This is the main function ############
getcrosses <- function(xp1, xp2, nsamps=nsamps,params,dmax=NULL,opt=c(0,0))
{
  # include only subs before "time=tstop"
  dmax <- ifelse(is.null(dmax),nrow(xp1), dmax+1)
  tstop <- min(xp1[dmax,'fix_time'], xp2[dmax,'fix_time'])
  xp1trim <- xp1[which(xp1$fix_time <= tstop),]
  xp2trim <- xp2[which(xp2$fix_time <= tstop),]
  t1 <- xp1trim[nrow(xp1trim),'fix_time']
  t2 <- xp2trim[nrow(xp2trim),'fix_time']
  if(any(is.na(c(t1,t2)))){print("t1 or t2 is NA")}

  #select columns with muteffects and dom coeff for each trait
  mc1 <- grep('mut_effect',colnames(xp1trim))
  mc2 <- grep('mut_effect',colnames(xp2trim))
  dc1 <- grep('dom',colnames(xp1trim))
  dc2 <- grep('dom',colnames(xp2trim))
  
  n <- length(mc1)
  if(length(mc2)!=n) stop('Inconsistent number of traits in input files!')
  
  # Get the parental trait values, and the effects of each mutation
  zp1 <- apply(xp1trim[,mc1],FUN=sum,MARGIN=2)
  zp2 <- apply(xp2trim[,mc2],FUN=sum,MARGIN=2)
  zm <- rbind(xp2trim[,mc2],-xp1trim[,mc1])
  doms <- rbind(xp2trim[,dc2],1-xp1trim[,dc1])
  

  d <- dim(zm)[1]
  d1 <- dim(xp1trim)[1]
  d2 <- dim(xp2trim)[1]
  cat('fixations = ',d1,' + ',d2,' = ',d,'. ')
   
	# Used in the fitness function
	k <- params['k']

  # (1) Get the fitness for the parental population, and the F1 
	zF1 <- zp1+apply(zm*doms,MARGIN=2,FUN=sum)
 	
	SP1 <- sum((zp1-opt)^2)
 	SP2 <- sum((zp2-opt)^2)
 	SF1 <- sum((zF1-opt)^2)
 	
 	wP1 <- unname(w(SP1,k))
 	wP2 <- unname(w(SP2,k))
 	wF1 <- unname(w(SF1,k))
 	
 	
 	## (2) Create the gametes from the F1 and P1, and stats for BCs and the F2
 	gamete_matrix1  <- matrix(rbinom(nsamps*d,1,0.5),d,nsamps)
 	gamete_matrix2  <- matrix(rbinom(nsamps*d,1,0.5),d,nsamps)
 	gamete_matrixP1 <- matrix(0,d,1)
 	gamete_matrixP2 <- matrix(1,d,1)
 	cg <- col(gamete_matrix1)
 	spgm1F1  <- split(gamete_matrix1,cg)
 	spgm2F1  <- split(gamete_matrix2,cg)
 	spgmP1 <- split(gamete_matrixP1,col(gamete_matrixP1))
 	spgmP2 <- split(gamete_matrixP2,col(gamete_matrixP2))
 	F2ind <- mapply(FUN=rbind,spgm1F1,spgm2F1,SIMPLIFY =F)
 	spgm1F2 <- lapply(F2ind,FUN=function(x) {apply(x,2,FUN=sample,size=1)})
 	spgm2F2 <- lapply(F2ind,FUN=function(x) {apply(x,2,FUN=sample,size=1)})
 	F3ind <- mapply(FUN=rbind,spgm1F2,spgm2F2,SIMPLIFY =F)
 	spgm1F3 <- lapply(F3ind,FUN=function(x) {apply(x,2,FUN=sample,size=1)})
 	spgm2F3 <- lapply(F3ind,FUN=function(x) {apply(x,2,FUN=sample,size=1)})
 	F4ind <- mapply(FUN=rbind,spgm1F3,spgm2F3,SIMPLIFY =F)
 	spgm1F4 <- lapply(F4ind,FUN=function(x) {apply(x,2,FUN=sample,size=1)})
 	spgm2F4 <- lapply(F4ind,FUN=function(x) {apply(x,2,FUN=sample,size=1)})
 	F5ind <- mapply(FUN=rbind,spgm1F4,spgm2F4,SIMPLIFY =F)
 	spgm1F5 <- lapply(F5ind,FUN=function(x) {apply(x,2,FUN=sample,size=1)})
 	spgm2F5 <- lapply(F5ind,FUN=function(x) {apply(x,2,FUN=sample,size=1)})
 	F6ind <- mapply(FUN=rbind,spgm1F5,spgm2F5,SIMPLIFY =F)
 	spgm1F6 <- lapply(F6ind,FUN=function(x) {apply(x,2,FUN=sample,size=1)})
 	spgm2F6 <- lapply(F6ind,FUN=function(x) {apply(x,2,FUN=sample,size=1)})
 	F7ind <- mapply(FUN=rbind,spgm1F6,spgm2F6,SIMPLIFY =F)
 	spgm1F7 <- lapply(F7ind,FUN=function(x) {apply(x,2,FUN=sample,size=1)})
 	spgm2F7 <- lapply(F7ind,FUN=function(x) {apply(x,2,FUN=sample,size=1)})
 	
 	SF2  <- mapply(FUN=getS,spgm1F1,spgm2F1, MoreArgs=list(zm=zm,doms=doms,zstart=zp1))
 	SF3  <- mapply(FUN=getS,spgm1F2,spgm2F2, MoreArgs=list(zm=zm,doms=doms,zstart=zp1))
 	SF4  <- mapply(FUN=getS,spgm1F3,spgm2F3, MoreArgs=list(zm=zm,doms=doms,zstart=zp1))
 	SF5  <- mapply(FUN=getS,spgm1F4,spgm2F4, MoreArgs=list(zm=zm,doms=doms,zstart=zp1))
 	SF6  <- mapply(FUN=getS,spgm1F5,spgm2F5, MoreArgs=list(zm=zm,doms=doms,zstart=zp1))
 	SF7  <- mapply(FUN=getS,spgm1F6,spgm2F6, MoreArgs=list(zm=zm,doms=doms,zstart=zp1))
 	SF8  <- mapply(FUN=getS,spgm1F7,spgm2F7, MoreArgs=list(zm=zm,doms=doms,zstart=zp1))
 	## 	SB1  <- mapply(FUN=getS,spgm1F1,spgmP1, MoreArgs=list(zm=zm,doms=doms,zstart=zp1))
 	## 	SB2  <- mapply(FUN=getS,spgm1F1,spgmP2, MoreArgs=list(zm=zm,doms=doms,zstart=zp1))

  wF2 <- w(SF2,k)
  wF3 <- w(SF3,k)
  wF4 <- w(SF4,k)
  wF5 <- w(SF5,k)
  wF6 <- w(SF6,k)
  wF7 <- w(SF7,k)
  wF8 <- w(SF8,k)
##  wB1 <- w(SB1,k)
##  wB2 <- w(SB2,k)
	
  #### CHANGED - bracket ####
	return(list(wP1=wP1,wP2=wP2,wF1=wF1,wF2=wF2,wF3=wF3,wF4=wF4,wF5=wF5,wF6=wF6,wF7=wF7,wF8=wF8,zm=zm,d=d))
}


##############################################################################################################################
##############################################################################################################################
#####################   make crosses                                                               ###########################
##############################################################################################################################
##############################################################################################################################

# Use this if you want to make crosses from lots of files
##N,L,U,k,s,n,Fd,q,MVN,diff
# myfunction <- function(U=0.01,k=2,s=0.01,n=2,Fd=0,q=0.5,MVN=1,N=1000,E=2,dmax=100,theta=0)
# {
# 	L <- -1 # Reduced recombination not yet implemented...
# 
# 	#setwd("~/Desktop/Hilde-predpaper/Crosses")
# 	#source('01processfilenames.R')
# 
# 	setwd(filedir)
#   
# 	params <- c(N,L,U,k,s,n,Fd,q,MVN,E)
# 	names(params) <- c('N','L','U','k','s','n','F','q','M','E')
# 	
# 	pat <- paste0(sapply(1:length(params),FUN=function(x) {paste0(names(params)[x],as.character(params)[x])}),collapse="_")
# 	pat <- gsub(pat,pattern='1e-04', replacement='0.0001')
# 	
# 	nsamps <- 1e4
# 
# 	lall <- list.files(filedir,full.names=TRUE);l <- lall[grep(lall,pattern=pat)]
# 	cat(length(l), ' out of ',length(lall),' files meet criteria. \n')
# 	
# 	# Read in the files 
# 	mf <- c(1,2)#grep(getfile.norep(params),l)[1:2] #rep3 is better than rep2.
# 	xp1 <- read.table(l[[mf[1]]],header=T)
# 	xp2 <- read.table(l[[mf[2]]],header=T)
# 
# 	# Get the info on crosses
# 	res <- getcrosses(xp1, xp2,nsamps= nsamps,params=params,dmax=dmax)
# 
# 	return(res)
# 
# 
# }
### Run the function...
#res <- myfunction(dmax=10,theta=0,U=0.0001,s=0.01)

N=10 #population size
L=-1 #free recombination
U=0.1 #mutation rate
k=2 #landscape curvature
s=0.001 #mutation effect size in optimal bg
n=2 #nb of traits
Fd=0 #dominance variance
q=0.5 #mean dominance coef
MVN=1 #distribution of mutation effects
E=0 #environmental scenario
nsamps=50 #nb of hybrids per cross type
dmax=100 #maximum nb of subs to consider

params <- c(N,L,U,k,s,n,Fd,q,MVN,E)
names(params) <- c('N','L','U','k','s','n','F','q','M','E')
##############################################################################################################################################

xp1 = read.csv("~/Desktop/Simulation_output/Diploid/Parameter4/P4R1.txt",sep="\t")#substitution file for P1
xp2 = read.csv("~/Desktop/Simulation_output/Diploid/Parameter4/P4R3.txt",sep="\t")#substitution file for P2
  

res <- getcrosses(xp1,xp2,nsamps=nsamps,params,dmax,opt=rep(0,params['n']))
Dip_P4R1R3_raw <- res[grep("w",names(res))]

##### save data for compiled data plot (havnet calculated mean of this replica yet)
save(Dip_P4R1R3_raw, file="Dip_P4R1R3raw.RData")




#####plotting of the individual replica

# alpha = 0.05
# degrees_of_freedom = nsamps - 1
# t_score = qt(p=alpha/2, df=degrees_of_freedom,lower.tail=F)
# calc_data <- sapply(Dip_P2R2R2_raw,FUN=function(x) {c(m=mean(x),sd=sd(x),se=sd(x)/sqrt(nsamps),me=t_score*sd(x)/sqrt(nsamps),up=mean(x)+t_score*sd(x)/sqrt(nsamps),lb=mean(x)-t_score*sd(x)/sqrt(nsamps))})
# 
# plot_data <- cbind(x=c(1:10),as.data.frame(t(calc_data)))
# 
# ####use ggplot to plot the data
# require(ggplot2)
# ggplot(plot_data, aes(x = plot_data$x, y = plot_data$m)) +
#   geom_point(size = 4) +
#   geom_errorbar(aes(ymax = plot_data$up, ymin = plot_data$lb)) +
#   labs(y= "Fitness",x="") +
#   scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10"),labels=c("P1","P2","F1","F2","F3","F4","F5","F6","F7","F8"))
# 
# complete_data <- cbind(phenotype=c("trait"),ploidy=c("2X"),line=c("P2R7R4"),type=c("P1","P2","F1","F2","F3","F4","F5","F6","F7","B1","B2"),as.data.frame(t(calc_data)),n=c("100"))
# library("writexl")
# write_xlsx(complete_data, "~/Desktop/Cross_output/F1-7_output/Dip_F1-7/Dip_P2R7R8.xlsx")
# 
# 


#save(plot_data, file="Dip_P6R7R8.RData")


########original plotting method
# 
# #mean of the list
# wmeans <- c(wP1m=res$wP1,wP2m=res$wP2,wF1m=res$wF1,wF2m=mean(res$wF2),wB1m=mean(res$wB1),wB2m=mean(res$wB2))
# 
# #standard deviation of the list
# wstd <- c(wP1std=0,wP2std=0,wF1std=0,wF2std=sd(res$wF2),wB1std=sd(res$wB1),wB2std=sd(res$wB2))
# #standard error of the list
# wse <- c(wP1se=0,wP2se=0,wF1se=0,wF2se=wstd["wF2std"]/sqrt(nsamps),wB1se=wstd["wB1std"]/sqrt(nsamps),wB2se=wstd["wB2std"]/sqrt(nsamps))
# names(wse)[4] <- "wF2se"
# names(wse)[5] <- "wB1se"
# names(wse)[6] <- "wB2se"
# #determine t-score
# alpha = 0.05
# degrees_of_freedom = nsamps - 1
# t_score = qt(p=alpha/2, df=degrees_of_freedom,lower.tail=F)
# #calculate margin of error
# wme <- c(wP1me=0,wP2me=0,wF1me=0,wF2me=t_score*wse["wF2se"],wB1me=t_score*wse["wB1se"],wB2me=t_score*wse["wB2se"])
# names(wme)[4] <- "wF2me"
# names(wme)[5] <- "wB1me"
# names(wme)[6] <- "wB2me"
# 
# #calculate upper and lower bound
# wub <- c(wP1ub=wmeans["wP1m"],wP2ub=wmeans["wP2m"],wF1ub=wmeans["wF1m"],wF2ub=wmeans["wF2m"]+wme["wF2me"],wB1ub=wmeans["wB1m"]+wme["wB1me"],wB2ub=wmeans["wB2m"]+wme["wB2me"])
# wlb <- c(wP1lb=wmeans["wP1m"],wP2lb=wmeans["wP2m"],wF1lb=wmeans["wF1m"],wF2lb=wmeans["wF2m"]-wme["wF2me"],wB1lb=wmeans["wB1m"]-wme["wB1me"],wB2lb=wmeans["wB2m"]-wme["wB2me"])
# names(wub)[1] <- "wP1ub"
# names(wub)[2] <- "wP2ub"
# names(wub)[3] <- "wF1ub"
# names(wub)[4] <- "wF2ub"
# names(wub)[5] <- "wB1ub"
# names(wub)[6] <- "wB2ub"
# 
# names(wlb)[1] <- "wP1lb"
# names(wlb)[2] <- "wP2lb"
# names(wlb)[3] <- "wF1lb"
# names(wlb)[4] <- "wF2lb"
# names(wlb)[5] <- "wB1lb"
# names(wlb)[6] <- "wB2lb"
# 
# #combine the data needed for plotting
# plot_data <- data.frame(cbind(x=c(1:6),wmeans,wlb,wub))
# 
# 
