#####################################################################################
library(nloptr)
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
elvisutil <- function(Z, theta, Moment=function(U, Z, theta) (Z[,1]+U-Z[,2]*theta[1])*Z[,2], 
	rep, thetanuisin, gamin){
#######Metropolis Hastings Algorithm############################################
 	avg_mom <- function(theta, nuisgam){ 
 		if (getnbnuis()==0 & is.null(thetanuisin)){
 			gam <- nuisgam
 		} else {
 			thetanuis <- nuisgam[1:length(thetanuisin)]
 			gam <- nuisgam[length(thetanuisin)+1:dimf]
 		}
 		r <- -rep[1]+1
 		avg <- matrix(0, nrow=n, ncol=dimf)
 		#Initial Guess for Unobservables
 		un <- guess_un(Z, theta, gam)
		for (i in r:rep[2]){
			set.seed(123456+i)
			#Proposal distribution 
			try_un <- jump_un(un, Z, theta)
			#Acceptance Probability
			try_dens <- exp((Moment(try_un, Z, theta)%*%gam)-(Moment(un, Z, theta)%*%gam))*
			rho(try_un, Z, theta)/rho(un, Z, theta)
			#Create a set of indices for acceptance rule
			aindex <- which(runif(n) < try_dens)
			#Update unobservables that satisfy acceptance rule
			un[aindex] <- try_un[aindex]
			}
			if (i>0){
				avg <- avg+Moment(un, Z, theta)/rep[2]
			}
 		avg[is.na(avg)] <- 0
 		return(avg)
 	}
 ############Define GMM Objective Function#############################
 	gmm_objective <- function(theta, nuisgam){
 		avg <- avg_mom(theta, nuisgam)
 		gm <- colMeans(avg)
 		V <- (t(avg)%*%avg/n-gm%*%t(gm))
 		return(-t(gm)%*%pinv(V)%*%gm)
 	}
 #######Grid Search over theta (For GMM objective)
 		gamp <- sbplx(x0=c(thetanuisin, gamin), fn=function(z) gmm_objective(theta, z), control=list(ftol_abs=1e-3))$value
 		return(gamp)
}
######################################MCMC Rules ##################################
#Burn in and integration replications
rep <- c(5,10)
#Number of MC replications
nreps <- 2
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
#Indicate the number of nuisance paramters (to be profiled out)
getnbnuis <- function(){
	return(1)
}
###########################################DGP#####################################
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
thetao <- 1
#True value of nuisance parameter Theta[2]
thetanuiso <- 0.2
#Grid for Theta search
theta <- seq(0, 2, by=0.05)
#Matrix to store objective function values over grid
soln <- array(0, dim=c(length(theta), nreps))
#Time entire experiment
overall.start.time <- Sys.time()
for (j in 1:nreps){
	loop.time <- proc.time()
	set.seed(123456+j)
	x <- rnorm(n, mean=0, sd=1)
	dy <- 0.5*rnorm(n, mean=0, sd=1)
	yt <- thetao[1]*x+dy
	y <- res*floor(yt/res)
	Z <- cbind(y, x)
	dimf <- dim(Moment(U=runif(n), Z=Z, theta=c(thetao, thetanuiso)))[2]
	for (t in 1:length(theta)){
		soln[t,j] <- elvisutil(Z=Z, theta=theta[t], Moment=Moment, rep=rep, thetanuisin=0.2, gamin=rep(0, times=dimf))
	}
	print(proc.time()-loop.time)
	if (j==nreps){
		plot(theta, soln[,j], type="l", xlab=expression(theta), ylab=~hat(L[n])(theta), main="Objective Function for Interval-Valued Data Regression")
		abline(v=thetao, lty=3)
	}
}
par <- theta[apply(soln, 2, which.max)]
mean.theta <- mean(par)
se.theta <- sd(par)
print(Sys.time()-overall.start.time)





