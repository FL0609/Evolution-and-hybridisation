############ Diploid F2 ############
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

plot_data <- data.frame(wF2m=numeric(199),wF2ub=numeric(199),wF2lb=numeric(199))

for(dmax in 1:199) {
  

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
 	spgm1  <- split(gamete_matrix1,cg)
 	spgm2  <- split(gamete_matrix2,cg)
 	spgmP1 <- split(gamete_matrixP1,col(gamete_matrixP1))
 	spgmP2 <- split(gamete_matrixP2,col(gamete_matrixP2))
 	SF2  <- mapply(FUN=getS,spgm1,spgm2, MoreArgs=list(zm=zm,doms=doms,zstart=zp1))
 	SB1  <- mapply(FUN=getS,spgm1,spgmP1, MoreArgs=list(zm=zm,doms=doms,zstart=zp1))
 	SB2  <- mapply(FUN=getS,spgm1,spgmP2, MoreArgs=list(zm=zm,doms=doms,zstart=zp1))

  wF2 <- w(SF2,k)
  wB1 <- w(SB1,k)
  wB2 <- w(SB2,k)
	
  #### CHANGED - bracket ####
	return(list(wP1=wP1,wP2=wP2,wF1=wF1,wF2=wF2,wB1=wB1,wB2=wB2,zm=zm,d=d))
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
nsamps=100 #nb of hybrids per cross type
#dmax=100 #maximum nb of subs to consider
  
  
params <- c(N,L,U,k,s,n,Fd,q,MVN,E)
names(params) <- c('N','L','U','k','s','n','F','q','M','E')

xp1 = read.csv("~/Desktop/Simulation_output/Diploid/Parameter6/P6R7.txt",sep="\t")#substitution file for P1
xp2 = read.csv("~/Desktop/Simulation_output/Diploid/Parameter6/P6R8.txt",sep="\t")#substitution file for P2
  

res <- getcrosses(xp1,xp2,nsamps=nsamps,params,dmax,opt=rep(0,params['n']))

#mean of the list
wF2m <- mean(res$wF2)

#standard deviation of the list
wF2std <- sd(res$wF2)
#standard error of the list
wF2se <- wF2std/sqrt(nsamps)

#determine t-score
alpha = 0.05
degrees_of_freedom = nsamps - 1
t_score = qt(p=alpha/2, df=degrees_of_freedom,lower.tail=F)
#calculate margin of error
wF2me <- t_score*wF2se

#calculate upper and lower bound
wF2ub <- wF2m + wF2me
wF2lb <- wF2m - wF2me

#combine the data needed for plotting

plot_data$wF2m[dmax] <- wF2m
plot_data$wF2ub[dmax] <- wF2ub
plot_data$wF2lb[dmax] <- wF2lb

}





#use ggplot to plot the data
plot_data_Dip_F2 <- data.frame(cbind(x=c(1:199),plot_data))
colors <- c("Diploid" = "black", "Tetraploid" = "red")

require(ggplot2)
ggplot(plot_data_Dip_F2, aes(x = plot_data_Dip_F2$x, y = plot_data_Dip_F2$wF2m, color = "Diploid")) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymax = plot_data$wF2ub, ymin = plot_data$wF2lb, color = "Diploid")) + 
  geom_point(size = 0.5, data = plot_data_Tetra_F2, aes(x = plot_data_Tetra_F2$x, y = plot_data_Tetra_F2$wF2m, color = "Tetraploid")) +
  geom_errorbar(aes(ymax = plot_data_Tetra_F2$wF2ub, ymin = plot_data_Tetra_F2$wF2lb, color = "Tetraploid"))+
  labs(y= "Fitness",x="Dmax",colour="Ploidy")


save(plot_data_Dip_F2, file="plot_data_Dip_F2.RData")
