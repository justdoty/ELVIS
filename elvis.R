#####################################################################################
library(nloptr)
library(pracma)
library(GenSA)
#####Parallelize####################
# cl <- makeCluster(4)
######################################README#########################################
# Routine implementing ELVIS in R
# The main routine the user should call is elvis at the end of this file

# general naming conventions
# un:unobservables
# theta: parameter vector
# nuis: nuisance parameters (profiled over)
# gam: tilting parameter gamma (see paper for details)

# Calculate average of moment function over the unobservables for each individual
# Moment: pointer to moment function
# theta: parameter vector, including nuisance parameters
# gam: tilting parameter gamma
# rep: array controling replications:
# 	rep[1]: number of equilibration steps
# 	rep[2]: number of averaging steps
# Z: Observables
# dimf: number of moment equations. Needs to match the dimension of Moment
#thetain: Initial values of theta to search over
#gamin: Initial values of tilting parameter, gamma to search over (must match the dimension of dimf)
####################ELVIS######################################################
elvisutil <- function(Z, thetain, gamin, rep, Moment=function(U, Z, theta) (Z[,1]+U-Z[,2]*theta[1])*Z[,2],
	guess_un=function(Z, theta, gam) matrix(runif(n)), jump_un=function(U, Z, theta) matrix(runif(n)), 
	rho=function(U, Z, theta) matrix(rep(1, n))){
#######Metropolis Hastings Algorithm############################################
 	avg_mom <- function(thetagam){ 
 		theta <- thetagam[1:length(thetain)]
 		gam <- thetagam[(length(thetain)+1):(length(thetain)+dimf)]
 		r <- -rep[1]+1
 		#Store average as a for accepted draws
 		avg <- matrix(0, nrow=n, ncol=dimf)
 		#Store accepted values of the vector of unobservables
 		#This is useful for MCMC diagnostics as a function of a sequence of
 		#parameter values for theta and gamma (esp. the optimal values)
 		unmat <- matrix(0, nrow=n, ncol=rep[2])
 		#Initial Guess for Unobservables
 		un <- guess_un(Z, theta, gam)
		for (i in r:rep[2]){
			#Proposal distribution 
			try_un <- jump_un(un, Z, theta)
			#Acceptance Probability
			try_dens <- exp((Moment(try_un, Z, theta)%*%gam)-(Moment(un, Z, theta)%*%gam))*
			rho(try_un, Z, theta)/rho(un, Z, theta)
			#Create a set of indices for acceptance rule
			aindex <- which(runif(n) < try_dens)
			#Update unobservables that satisfy acceptance rule
			un[aindex] <- try_un[aindex]
			if (i>0){
				unmat[,i] <- un
				avg <- avg+Moment(un, Z, theta)/rep[2]
			}
		}
		avg[is.na(avg)] <- 0
 		return(avg)
 	}

#############Define MM Objective Function##############################
#Searching for gamma is convex by construction so we use numerical methods
#Analytical jacobian for newton-raphson can be constructed
 	GMM_objective <- function(thetagam){
 		mom <- avg_mom(thetagam)
 		obj <- nrow(mom)*colMeans(mom)%*%inv(var(mom))%*%as.matrix(colMeans(mom))
 		return(-obj)
 	}
 	GMM_profile <- function(theta){
 		#Subplex method
 		gam <- bobyqa(x0=gamin, fn=function(z) GMM_objective(c(theta, z)))
 		gamsoln <- gam$value
 		return(gamsoln)
 		#MM method
 		# gamsol <- newtonsys(Ffun=function(z) GMM_objective(c(theta, z)), x0=gamin, maxiter=500)$zero
 		# return(GMM_objective(c(theta, gamsol)))
 	}
 #Optimize over theta using GMM
 results <- bobyqa(x0=thetain, fn=GMM_profile, lower=thetao-1, upper=thetao+1, control=list(maxeval=100))
 print(results)
 thetasoln <- results$par
 return(thetasoln)
}
######################################MCMC Rules ##################################
#Burn in and integration replications
rep <- c(50,500)
#Number of MC replications
nreps <- 1
#Initial Guess for the unobservable
guess_un <- function(Z, theta, gam){
	return(matrix(runif(n)))
}
#Rule for trial jumps in Metropolis Monte Carlo
jump_un <- function(U, Z, theta){
	return(matrix(runif(n)))
}
#User Specified Measure Rho
rho <- function(U, Z, theta){
	return(matrix(rep(1, n)))
}
###########################################DGP#####################################
#Interval Valued Data Example
###########################Moment Function########################
Moment <- function(U, Z, theta){
	M <- cbind(Z[,1]+U-Z[,2]*theta[1], (Z[,1]+U-Z[,2]*theta[1])*Z[,2], 
	((Z[,1]+U-Z[,2]*theta[1])^2)-theta[2], (Z[,2]^2)*(((Z[,1]+U-Z[,2]*theta[1])^2)-theta[2]))
	return(M)
}
#Sample Size
n <- 250
res <- 1
#True value of Theta[1]
thetao <- c(1, 0.2)
#Matrix to store objective function values over grid
soln <- array(0, dim=c(nreps, length(thetao)))
#Time entire experiment
overall.start.time <- Sys.time()
for (j in 1:nreps){
	print(j)
	loop.time <- proc.time()
	set.seed(123456+j)
	x <- rnorm(n, mean=0, sd=1)
	dy <- 0.5*rnorm(n, mean=0, sd=1)
	yt <- thetao[1]*x+dy
	y <- res*floor(yt/res)
	Z <- cbind(y, x)
	dimf <- dim(Moment(U=runif(n), Z=Z, theta=thetao))[2]
	soln[j,] <- elvisutil(Z=Z, thetain=thetao, gamin=rep(0, times=dimf), rep=rep, Moment=Moment, guess_un=guess_un,
		jump_un=jump_un, rho=rho)
	print(proc.time()-loop.time)
}
mean.theta <- colMeans(soln)
se.theta <- apply(soln, 2, sd)
print(mean.theta)
print(se.theta)
print(Sys.time()-overall.start.time)





