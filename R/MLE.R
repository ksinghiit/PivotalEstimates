## ...........block progressive censored for shape-scale family distribution............................

rm(list=ls(all=T))
rm(list=ls(all=TRUE))
start.time <- Sys.time()
library(FixedPoint)

alp <- 1.2; bet <- 1.5;

N <- c(50, 40, 30) # n1, n2, n3
M <- c(45, 33, 25) # s1, s2, s3 


 
x0=0.5;
N0=1000; NN=2000;los=0.05
times=1;
k <- length(M)
R <- list()
for(r in 1:k){
	R[[r]] <- c(rep(0,(M[r]-1)),(N[r]-M[r])) 
}

sam_generation_progressive <- function(nn, mm,al, bt,RR){
	ww <- runif(mm);
	vv <- numeric(0);
	ui <- numeric(0);
	for(i1 in 1:mm){
		vv[i1] <- (ww[i1])^(1/(i1 + sum(RR[(mm-i1+1):mm])));
	}
	for(i2 in 1:mm){
		ui[i2] <- 1-prod(vv[(mm-i2+1):mm]);
	}
	xs <- ((-1/bt)*log(1-ui))^(1/al); #quantile of Weibull distribution
	return(xs);
}
xsamp <- sam_generation_progressive(N[1],M[1], alp, bet, R[[1]]);xsamp; # sample for progressive censoring

smpl_list <- list()
for(sm in 1:k){
	smpl_list[[sm]] <- sam_generation_progressive(N[sm],M[sm],alp,bet,R[[sm]]) + abs(rnorm(M[sm], 0, 0.001))
}
x <- smpl_list; #x;

#-----Psi_i(alp) function------#
ww <- function(al,rr,gt){
	wt <- dwt <- ddwt <- numeric();	
		for(j in 1:k){
			wt[j] <- sum((rr[[j]]+1)*(gt[[j]]^al)) # w_i(alpha);
			dwt[j] <- sum((rr[[j]]+1)*(gt[[j]]^al)*(log(gt[[j]])))# first derivative of w_i(alp);
			ddwt[j] <- sum((rr[[j]]+1)*(gt[[j]]^al)*((log(gt[[j]]))^2))# second derivative of w_i(alp);
		}
	wwt <- list();
	wwt[[1]] <- wt; #W_i
	wwt[[2]] <- dwt;#Dw_i
	wwt[[3]] <- ddwt;#DDw_i
	return(wwt)
}
gtt <- list()
for(i in 1:k){
	gtt[[i]] <- x[[i]];#Weibull distribution
}
gtt;

#....alpha estimates using Fixed point
alp_est_fun <- function(ap,rr,gt,mm){
	wtt <- ww(ap,rr,gt);
	wt <- wtt[[1]];
	dwt <- wtt[[2]];
	smj <- numeric()
	smi1 <- smi2 <- 0;
	for(i in 1:k){
		smj[i] <- 0;
		for(j in 1:mm[i]){
			smj[i] <- smj[i] + (log(gt[[i]][j]));
		}
		smi1 <- smi1 + smj[i];
		smi2 <- smi2 + (mm[i]*(dwt[i]/wt[i]));	
	}
	st <- -((sum(mm))/(smi1-smi2));
	return(st)
}

    # g(ap): returns the updated alpha given current ap
	g_alpha <- function(ap){alp_est_fun(ap, rr = R, gt = gtt, mm = M)}
	
	# fixed point iteration starting from ap0
	mle_alp <- fixed_point_est(gfun=g_alpha, init = 1)   # 1 is initial guess
	mle_alp;
