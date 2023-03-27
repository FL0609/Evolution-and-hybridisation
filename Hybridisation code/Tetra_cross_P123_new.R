#library(scales)
############  Obtain a file name from a list of parameters
getfile.norep <- function(ff) {
	return(paste0('res_N',ff[1],'_L',ff[2],'_U',formatC(ff[3]),'_k',ff[4],'_a',formatC(ff[5]),'_n',ff[6],'_F',ff[7],'_q',ff[8],'_MVN',ff[9],'_diff',ff[10],'_rep'))
}

############ Obtain squared distance from the optimum for a single genotype
#changed - added g3, g4; deleted hom and het ====
getSP123 <- function(g1,g2,g3,g4,zm12,zm13,zstart1,opt=c(1,0)) 
{ 
  { 
    m1 <- zm12*g1*0.25
    m2 <- zm12*g2*0.25
    m3 <- zm13*g3*0.25
    m4 <- zm13*g4*0.25
    
    z <- zstart1*0.25+apply(m1,MARGIN=2,FUN=sum)+
      zstart1*0.25+apply(m2,MARGIN=2,FUN=sum)+
      zstart1*0.25+apply(m3,MARGIN=2,FUN=sum)+
      zstart1*0.25+apply(m4,MARGIN=2,FUN=sum)
    return(sum(z^2)) # Sum across traits
  }	
}	

getSP234 <- function(g1,g2,g3,g4,zm24,zm34,zstart4,opt=c(1,0)) 
{ 
  { 
    m1 <- zm24*g1*0.25
    m2 <- zm24*g2*0.25
    m3 <- zm34*g3*0.25
    m4 <- zm34*g4*0.25
    
    z <- zstart4*0.25+apply(m1,MARGIN=2,FUN=sum)+
      zstart4*0.25+apply(m2,MARGIN=2,FUN=sum)+
      zstart4*0.25+apply(m3,MARGIN=2,FUN=sum)+
      zstart4*0.25+apply(m4,MARGIN=2,FUN=sum)
    return(sum(z^2)) # Sum across traits
  }	
}	

getSP1234 <- function(g1,g2,g3,g4,zm12,zm34,zstart1,zstart4,opt=c(1,0)) 
{ 
    m1 <- zm12*g1*0.25
    m2 <- zm12*g2*0.25
    m3 <- zm34*g3*0.25
    m4 <- zm34*g4*0.25
    
    z <- zstart1*0.25+apply(m1,MARGIN=2,FUN=sum)+
      zstart1*0.25+apply(m2,MARGIN=2,FUN=sum)+
      zstart4*0.25+apply(m3,MARGIN=2,FUN=sum)+
      zstart4*0.25+apply(m4,MARGIN=2,FUN=sum)
    return(sum(z^2)) # Sum across traits
}	

getSP12 <- function(g1,g2,g3,g4,zm12,zstart1,opt=c(1,0)) 
{ 
  
  m <- rbind(zm12*g1*0.25,zm12*g2*0.25,zm12*g3*0.25,zm12*g4*0.25)
  if(!is.null(dim(m)))
  {	
    z <- zstart1+apply(m,MARGIN=2,FUN=sum)
    return(sum(z^2)) # Sum across traits
  }
  z <- (zstart1+sum(m))-opt
  return(sum(z^2))
}	

getSP13 <- function(g1,g2,g3,g4,zm13,zstart1,opt=c(1,0)) 
{ 
  
  m <- rbind(zm13*g1*0.25,zm13*g2*0.25,zm13*g3*0.25,zm13*g4*0.25)
  if(!is.null(dim(m)))
  {	
    z <- zstart1+apply(m,MARGIN=2,FUN=sum)
    return(sum(z^2)) # Sum across traits
  }
  z <- (zstart1+sum(m))-opt
  return(sum(z^2))
}	

getSP34 <- function(g1,g2,g3,g4,zm34,zstart4,opt=c(1,0)) 
{ 
  
  m <- rbind(zm34*g1*0.25,zm34*g2*0.25,zm34*g3*0.25,zm34*g4*0.25)
  if(!is.null(dim(m)))
  {	
    z <- zstart4+apply(m,MARGIN=2,FUN=sum)
    return(sum(z^2)) # Sum across traits
  }
  z <- (zstart4+sum(m))-opt
  return(sum(z^2))
}	

getSP24 <- function(g1,g2,g3,g4,zm24,zstart4,opt=c(1,0)) 
{ 
  
  m <- rbind(zm24*g1*0.25,zm24*g2*0.25,zm24*g3*0.25,zm24*g4*0.25)
  if(!is.null(dim(m)))
  {	
    z <- zstart4+apply(m,MARGIN=2,FUN=sum)
    return(sum(z^2)) # Sum across traits
  }
  z <- (zstart4+sum(m))-opt
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
  xp3trim <- xp3[which(xp3$fix_time <= tstop),]
  xp4trim <- xp4[which(xp4$fix_time <= tstop),]
  t1 <- xp1trim[nrow(xp1trim),'fix_time']
  t2 <- xp2trim[nrow(xp2trim),'fix_time']
  t3 <- xp2trim[nrow(xp3trim),'fix_time']
  t4 <- xp2trim[nrow(xp4trim),'fix_time']
  if(any(is.na(c(t1,t2)))){print("t1 or t2 is NA")}

  #select columns with muteffects and dom coeff for each trait (from input file?)
  mc1 <- grep('mut_effect',colnames(xp1trim))
  mc2 <- grep('mut_effect',colnames(xp2trim))
  mc3 <- grep('mut_effect',colnames(xp3trim))
  mc4 <- grep('mut_effect',colnames(xp4trim))
  dc1 <- grep('dom',colnames(xp1trim))
  dc2 <- grep('dom',colnames(xp2trim))
  dc3 <- grep('dom',colnames(xp3trim))
  dc4 <- grep('dom',colnames(xp4trim))
  
  n <- length(mc1)
  if(length(mc2)!=n) stop('Inconsistent number of traits in input files!')
  
  # Get the parental trait values, and the effects of each mutation
  zp1 <- apply(xp1trim[,mc1],FUN=sum,MARGIN=2)
  zp2 <- apply(xp2trim[,mc2],FUN=sum,MARGIN=2)
  zp3 <- apply(xp3trim[,mc3],FUN=sum,MARGIN=2)
  zp4 <- apply(xp4trim[,mc4],FUN=sum,MARGIN=2)
  zp4 <- apply(xp4trim[,mc4],FUN=sum,MARGIN=2)
  zm12 <- rbind(xp2trim[,mc2],-xp1trim[,mc1])
  zm13 <- rbind(xp3trim[,mc3],-xp1trim[,mc1])
  zm34 <- rbind(xp3trim[,mc3],-xp4trim[,mc4])
  zm24 <- rbind(xp2trim[,mc2],-xp4trim[,mc4])
  doms12 <- rbind(xp2trim[,dc2],1-xp1trim[,dc1])
  doms13 <- rbind(xp3trim[,dc3],1-xp1trim[,dc1])
  doms34 <- rbind(xp3trim[,dc3],1-xp4trim[,dc4])
  doms24 <- rbind(xp2trim[,dc2],1-xp4trim[,dc4])

  d12 <- dim(zm12)[1]
  d13 <- dim(zm13)[1]
  d34 <- dim(zm34)[1]
  d24 <- dim(zm24)[1]
  d1 <- dim(xp1trim)[1]
  d2 <- dim(xp2trim)[1]
  d3 <- dim(xp3trim)[1]
  d4 <- dim(xp4trim)[1]
  cat('fixations P12 = ',d1,' + ',d2,' = ',d12,'. ',
      'fixations P13 = ',d1,' + ',d3,' = ',d13,'.',
      'fixations P34 = ',d3,' + ',d4,' = ',d34,'.',
      'fixations P24 = ',d2,' + ',d4,' = ',d24,'.')
   
	# Used in the fitness function
	k <- params['k']

  # (1) Get the fitness for the parental population, and the F1 
	#### NOT SURE - zp4+apply? ####
	zF1P12 <- zp1+apply(zm12*doms12,MARGIN=2,FUN=sum)
	zF1P13 <- zp1+apply(zm13*doms13,MARGIN=2,FUN=sum)
	zF1P34 <- zp4+apply(zm34*doms34,MARGIN=2,FUN=sum)
	zF1P24 <- zp4+apply(zm24*doms24,MARGIN=2,FUN=sum)
	
	SP1 <- sum((zp1-opt)^2)
	SP2 <- sum((zp2-opt)^2)
	SP3 <- sum((zp3-opt)^2)
	SP4 <- sum((zp4-opt)^2)
	SF1P12 <- sum((zF1P12-opt)^2)
	SF1P13 <- sum((zF1P13-opt)^2)
	SF1P34 <- sum((zF1P34-opt)^2)
	SF1P24 <- sum((zF1P24-opt)^2)
	
	wP1 <- unname(w(SP1,k))
	wP2 <- unname(w(SP2,k))
	wP3 <- unname(w(SP3,k))
	wP4 <- unname(w(SP4,k))
	wF1P12 <- unname(w(SF1P12,k))
	wF1P13 <- unname(w(SF1P13,k))
	wF1P34 <- unname(w(SF1P34,k))
	wF1P24 <- unname(w(SF1P24,k))

 	
 	## (2) Create the gametes from the F1 and P1, and stats for BCs and the F2
 	#changed - gamete matrix ====
 	gamete_matrix1a  <- matrix(rbinom(nsamps*d12,1,0.5),d12,nsamps)
	gamete_matrix1b  <- matrix(rbinom(nsamps*d12,1,0.5),d12,nsamps)
 	gamete_matrix2a  <- matrix(rbinom(nsamps*d13,1,0.5),d13,nsamps)
 	gamete_matrix2b  <- matrix(rbinom(nsamps*d13,1,0.5),d13,nsamps)
 	gamete_matrix3a  <- matrix(rbinom(nsamps*d34,1,0.5),d34,nsamps)
 	gamete_matrix3b  <- matrix(rbinom(nsamps*d34,1,0.5),d34,nsamps)
 	gamete_matrix1c  <- matrix(rbinom(nsamps*d12,1,0.5),d12,nsamps)
 	gamete_matrix1d  <- matrix(rbinom(nsamps*d12,1,0.5),d12,nsamps)
 	gamete_matrix2c  <- matrix(rbinom(nsamps*d13,1,0.5),d13,nsamps)
 	gamete_matrix2d  <- matrix(rbinom(nsamps*d13,1,0.5),d13,nsamps)
 	gamete_matrix3c  <- matrix(rbinom(nsamps*d34,1,0.5),d34,nsamps)
 	gamete_matrix3d  <- matrix(rbinom(nsamps*d34,1,0.5),d34,nsamps)
 	gamete_matrix4a  <- matrix(rbinom(nsamps*d24,1,0.5),d24,nsamps)
 	gamete_matrix4b  <- matrix(rbinom(nsamps*d24,1,0.5),d24,nsamps)
 	gamete_matrix4c  <- matrix(rbinom(nsamps*d24,1,0.5),d24,nsamps)
 	gamete_matrix4d  <- matrix(rbinom(nsamps*d24,1,0.5),d24,nsamps)

 	gamete_matrixP1 <- matrix(0,d1,1)
 	gamete_matrixP2 <- matrix(1,d2,1)
 	gamete_matrixP3 <- matrix(1,d3,1)
 	gamete_matrixP4 <- matrix(0,d4,1)
 	cg1a <- col(gamete_matrix1a)
 	cg1b <- col(gamete_matrix1b)
 	cg2a <- col(gamete_matrix2a)
 	cg2b <- col(gamete_matrix2b)
 	cg3a <- col(gamete_matrix3a)
 	cg3b <- col(gamete_matrix3b)
 	cg1c <- col(gamete_matrix1c)
 	cg1d <- col(gamete_matrix1d)
 	cg2c <- col(gamete_matrix2c)
 	cg2d <- col(gamete_matrix2d)
 	cg3c <- col(gamete_matrix3c)
 	cg3d <- col(gamete_matrix3d)
 	cg4a <- col(gamete_matrix4a)
 	cg4b <- col(gamete_matrix4b)
 	cg4c <- col(gamete_matrix4c)
 	cg4d <- col(gamete_matrix4d)
  
 	#changed - spgm3&4 ====
 	#changed - added SB3&SB4 NEED TO CHANGE ====
 	
 	spgm1a  <- split(gamete_matrix1a,cg1a)
 	spgm2a  <- split(gamete_matrix2a,cg2a)
 	spgm3a  <- split(gamete_matrix3a,cg3a)
 	spgm4a  <- split(gamete_matrix4a,cg4a)
 	
 	spgm1b  <- split(gamete_matrix1b,cg1b)
 	spgm2b  <- split(gamete_matrix2b,cg2b)
 	spgm3b  <- split(gamete_matrix3b,cg3b)
 	spgm4b  <- split(gamete_matrix4b,cg4b)
 	
 	spgm1c  <- split(gamete_matrix1c,cg1c)
 	spgm2c  <- split(gamete_matrix2c,cg2c)
 	spgm3c  <- split(gamete_matrix3c,cg3c)
 	spgm4c  <- split(gamete_matrix4c,cg4c)
 	
 	spgm1d  <- split(gamete_matrix1d,cg1d)
 	spgm2d  <- split(gamete_matrix2d,cg2d)
 	spgm3d  <- split(gamete_matrix3d,cg3d)
 	spgm4d  <- split(gamete_matrix4d,cg4d)

 	spgmP1 <- split(gamete_matrixP1,col(gamete_matrixP1))
 	spgmP2 <- split(gamete_matrixP2,col(gamete_matrixP2))
 	spgmP3 <- split(gamete_matrixP3,col(gamete_matrixP3))
 	spgmP4 <- split(gamete_matrixP4,col(gamete_matrixP4))
 	SF2P12  <- mapply(FUN=getSP12,spgm1a,spgm1b,spgm1c,spgm1d, MoreArgs=list(zm12=zm12,zstart1=zp1))
 	SF2P13  <- mapply(FUN=getSP13,spgm2a,spgm2b,spgm2c,spgm2d, MoreArgs=list(zm13=zm13,zstart1=zp1))
 	SF2P34  <- mapply(FUN=getSP34,spgm3a,spgm3b,spgm3c,spgm3d, MoreArgs=list(zm34=zm34,zstart4=zp4))
 	SF2P24  <- mapply(FUN=getSP24,spgm4a,spgm4b,spgm4c,spgm4d, MoreArgs=list(zm24=zm24,zstart4=zp4))
 	SF2P123  <- mapply(FUN=getSP123,spgm1a,spgm1b,spgm2a,spgm2b, MoreArgs=list(zm12=zm12,zm13=zm13,zstart1=zp1))
 	SF2P234  <- mapply(FUN=getSP234,spgm4a,spgm4b,spgm3a,spgm3b, MoreArgs=list(zm24=zm24,zm34=zm34,zstart4=zp4))
 	SF2P1234  <- mapply(FUN=getSP1234,spgm1a,spgm1b,spgm3a,spgm3b, MoreArgs=list(zm12=zm12,zm34=zm34,zstart1=zp1,zstart4=zp4))
 	# SB1  <- mapply(FUN=getS,spgm1,spgm2,spgmP1,spgmP1, MoreArgs=list(zm=zm,doms=doms,zstart=zp1))
 	# SB2  <- mapply(FUN=getS,spgm1,spgm2,spgmP2,spgmP2, MoreArgs=list(zm=zm,doms=doms,zstart=zp1))
 	wF2P12 <- w(SF2P12,k)
 	wF2P13 <- w(SF2P13,k)
 	wF2P34 <- w(SF2P34,k)
 	wF2P24 <- w(SF2P24,k)
 	wF2Pabc <- w(SF2P123,k)
 	wF2Pbcd <- w(SF2P234,k)
 	wF2P_all <- w(SF2P1234,k)
  # wB1 <- w(SB1,k)
  # wB2 <- w(SB2,k)
  
 	return(list(wP1=wP1,wP2=wP2,wP3=wP3,wP4=wP4,
 	            wF1P12=wF1P12,wF1P13=wF1P13,wF1P24=wF1P24,wF1P34=wF1P34,
 	            wF2P12=wF2P12,wF2P13=wF2P13,wF2P24=wF2P24,wF2P34=wF2P34,
 	            wF2Pabc=wF2Pabc,wF2Pbcd=wF2Pbcd,wF2P_all=wF2P_all,
 	            zm12=zm12,zm13=zm13,zm34=zm34,zm24=zm24,d12=d12,d13=d13,d34=d34,d24=d24))
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

N= #population size
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

xp1 = read.csv("~/Desktop/Simulation_output/Tetraploid/Parameter6/p6r1.txt",sep="\t")#substitution file for P1
xp2 = read.csv("~/Desktop/Simulation_output/Tetraploid/Parameter6/p6r2.txt",sep="\t")#substitution file for P2
xp3 = read.csv("~/Desktop/Simulation_output/Tetraploid/Parameter6/p6r3.txt",sep="\t")#substitution file for P3
xp4 = read.csv("~/Desktop/Simulation_output/Tetraploid/Parameter6/p6r4.txt",sep="\t")#substitution file for P4
  

res <- getcrosses(xp1,xp2,nsamps=nsamps,params,dmax,opt=rep(0,params['n']))

Tetra_p6r1234_raw <- res[grep("w",names(res))]


##### save data for compiled data plot (havnet calculated mean of this replica yet)
save(Tetra_p6r1234_raw, file="Tetra_p6r1234_raw.RData")

# ##### group data for as hybrid type (havnet calculated mean of this replica yet)
# Ps <- unlist(res[grep("wP",names(res))])
# F1s <- unlist(res[grep("wF1",names(res))])
# F2_2P <- unlist(c((res[grep("wF2P12",names(res))])))
# F2_3P <- unlist(res[grep("wF2Pabc",names(res))])
# F2_4P <- unlist(res[grep("wF2P_all",names(res))])
# 
# Tetra_p2r1234_raw <- list(Ps=unname(Ps),F1s=unname(F1s),F2_2P=unname(F2_2P),F2_3P=unname(F2_3P),F2_4P=unname(F2_4P))

# ######graph plotting of individual replicate
# alpha = 0.05
# degrees_of_freedom = nsamps - 1
# t_score = qt(p=alpha/2, df=degrees_of_freedom,lower.tail=F)
# 
# calc_data <- sapply(Tetra_p2r1234_raw,FUN=function(x) {c(m=mean(x),sd=sd(x),se=sd(x)/sqrt(nsamps),me=t_score*sd(x)/sqrt(nsamps),up=mean(x)+t_score*sd(x)/sqrt(nsamps),lb=mean(x)-t_score*sd(x)/sqrt(nsamps))})
# 
# plot_data <- cbind(x=c(1:15),as.data.frame(t(calc_data)))
# 
# #use ggplot to plot the data
# require(ggplot2)
# ggplot(plot_data, aes(x = plot_data$x, y = plot_data$m)) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(ymax = plot_data$up, ymin = plot_data$lb)) +
#   labs(y= "Fitness",x="") +
#   scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"),
#                    labels=c("P1","P2","P3","P4","F1(A*B)","F1(A*C)","F1(B*D)","F1(C*D)","F2(A*B)","F2(A*C)","F2(B*D)","F2(C*D)","F2(A*B)*(A*C)","F2(B*D)*(C*D)","F2(A*B)*(C*D)"))


# save(plot_data, file="Tetra_p6r5678.RData")


# ####original graph plotting code
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
# #use ggplot to plot the data
# require(ggplot2)
# ggplot(plot_data, aes(x = plot_data$x, y = plot_data$wmeans)) +
#   geom_point(size = 4) +
#   geom_errorbar(aes(ymax = plot_data$wub, ymin = plot_data$wlb)) +
#   labs(y= "fitness",x="") +
#   scale_x_discrete(limits=c("1","2","3","4","5","6"),labels=c("P1","P2","F1","F2","B1","B2"))





