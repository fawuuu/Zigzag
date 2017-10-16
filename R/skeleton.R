#' next_point_global
#'
#' a function which determines a new skeleton point for the Zig-zag process given global bounds on the directional
#' derivatives
#'
#'@param xi current skeleton position
#'@param theta current direction
#'@param t_flip most recent flip time
#'@param d number of spatial dimensions
#'@param derivatives function which takes position x as argument and outputs evaluation of directional derivatives at x
#'@param bounds vector of global bounds for directional derivatives
#'
#'

# unclear how much this function should include, see somments below in "skeleton"
#if subsampl=TRUE then derivatives as list of vector-valued functions
next_point_global <- function(xi, theta, t_flip, d, derivatives, bounds, subsample){

    flip <- FALSE

    while(flip == FALSE){

      test_taus <- rexp(d, bounds)
      tau <- min(test_taus)
      ind <- which.min(test_taus)

      xi <- xi + tau*theta
      t_flip <- t_flip + tau

      if(subsample == FALSE){

        if(runif(1) < max(0, theta[ind]*derivatives(xi)[ind])/bounds[ind]){
          theta[ind] <- -theta[ind]
          flip <- TRUE
        }

      }else{

        j <- sample(1, 1:length(derivatives))
        if(runif(1) < max(0, theta[ind]*derivatives[[j]](xi)[ind])/bounds[ind]){
          theta[ind] <- -theta[ind]
          flip <- TRUE
        }

      }

    }

    return(list(xi = xi, theta = theta, t_flip = t_flip))
}

#' next_point_hess
#'
#' A function which determines a new skeleton point for the Zig-zag process given Hessian-type bounds on the directional
#' derivatives
#'
#'@param xi current skeleton position
#'@param theta current direction
#'@param t_flip most recent flip time
#'@param d number of spatial dimensions
#'@param derivatives function which takes position x as argument and outputs evaluation of directional derivatives at x
#'@param a parameter defining hessian-type bound, calculated elsewhere
#'@param b parameter defining hessian-type bound, calculated elsewhere
#'
#'
next_point_hess <- function(xi, theta, t_flip, d, derivatives, a, b){

  flip <- FALSE

  while(flip == FALSE){

    test_taus <- -a/b + sqrt((a/b)^2 + 2/b*(-log(runif(d))))
    tau <- min(test_taus)
    ind <- which.min(test_taus)

    a <- a + b*tau
    xi <- xi + tau*theta
    t_flip <- t_flip + tau

    a_temp <- theta[ind]*derivatives(xi)[ind]

    if(runif(1) < max(0, a_temp)/max(0, a[ind])){

      theta[ind] <- -theta[ind]
      flip <- TRUE

    }

    a[ind] <- a_temp*(-1)^as.numeric(flip)

  }

  return(list(xi = xi, theta = theta, t_flip = t_flip, a = a))
}

#' next_point_lipschitz
#'
#'@param xi current skeleton position
#'@param theta current direction
#'@param t_flip most recent flip time
#'@param d number of spatial dimensions
#'@param derivatives list of functions which takes position x as argument and outputs evaluation of directional derivatives at x
#'@param bounds vector of global Lipschitz constants
#'@param xi0 reference point for the bounds
#'@param dervatives_ref vector of directional derivatives evaluated at the reference point
#'@param p distances in the p-norm
#'
#'
next_point_lipschitz <- function(xi, theta, t_flip, d, derivatives, bounds, xi0, derivatives_ref, p, subsample){

  flip <- FALSE
  derivatives_ref = theta * derivatives_ref
  derivatives_ref[which(derivatives_ref<0)] = 0
  if(p == Inf){
    a <- derivatives_ref + bounds*max(abs(xi-xi0))
    b <- bounds
  }else{
    a <- derivatives_ref + bounds*sum(abs(xi-xi0)^p)^(1/p)
    b <- bounds*d^(1/p)
  }

  while(flip == FALSE){

    test_taus <- -a/b + sqrt((a/b)^2 + 2/b*(-log(runif(d))))
    tau <- min(test_taus)
    ind <- which.min(test_taus)

    xi <- xi + tau*theta
    t_flip <- t_flip + tau

    if(subsample == FALSE){

      if(runif(1) < max(0, derivatives(xi)[ind])/max(0, a[ind]+b*tau)){
        theta[ind] <- -theta[ind]
        flip <- TRUE
      }
    }else{
      j <- sample(1, 1:length(derivatives))
      if(runif(1) < max(0, derivatives[[j]](xi)[ind])/max(0, a[ind]+b*tau)){
        theta[ind] <- -theta[ind]
        flip <- TRUE

    }


    if(p == Inf){
      a = derivatives_ref + bounds*max(abs(xi-xi0))
    }else{
      a = derivatives_ref + bounds*sum(abs(xi-xi0)^p)^(1/p)
    }

    }

  }

  return(list(xi = xi, theta = theta, t_flip = t_flip))
}

#' skeleton
#'
#' Function which produces the skeleton of a Zig-zag process
#'
#' @param xi an initial position
#' @param theta an initial direction
#' @param n number of skeleton points to be generated
#' @param dervatives function which takes position x as argument and outputs evaluation of directional derivatives at x
#' @param bounds ...
#' @param bound_type specifies type of bounds eg global, hessian
#' @param subsample boolean value; false by default; incompatible with hessian bounds;
#' @param xi0 reference value for
#' @param p choice of L_p space
#'
#'
skeleton <- function(xi, theta, n, derivatives, bounds, bound_type = "global", subsample = FALSE, xi0 = 0, p = 2){

  #set inital time, no. dimensions and record tables for skeleton
  t_flip <- 0
  d <- length(xi)

  xi_rec <- matrix(0,n,d)
  theta_rec <- matrix(0,n,d)
  t_flip_rec <- numeric(n)

  xi_rec[1, ] <- xi; theta_rec[1, ] <- theta

  #order of while & if loops here should be figured out - pointless to check bound_type at every while iteration,
  #but want to do if(test_point$flip == TRUE) part for every bound type, so also pointless to type this many times -
  #maybe arranging next_point functions better would help here

  #In response: Check whether new point is actually flipped in next_point function so that these actually give the next
  #point we want

  if (bound_type == "global"){

    for (i in 2:n){

    next_point <- next_point_global(xi_rec[i-1,], theta_rec[i-1,], t_flip_rec[i-1], d, derivatives, bounds, subsample)

        xi_rec[i, ] <- next_point$xi
        theta_rec[i, ] <- next_point$theta
        t_flip_rec[i] <- next_point$t_flip

    }
  }


  if (subsample == FALSE && bound_type == "hessian"){

    a = numeric(d); b = numeric(d)
      for(i in 1:d){
          a[i] = theta[i]*derivatives(xi)[i]
          b[i] = sqrt(d)*sqrt(sum(bounds[,i]^2))
      }

    for (i in 2:n){

      next_point <- next_point_hess(xi_rec[i-1,], theta_rec[i-1,], t_flip_rec[i-1], d, derivatives, a, b)

        xi_rec[i, ] <- next_point$xi
        theta_rec[i, ] <- next_point$theta
        t_flip_rec[i] <- next_point$t_flip
        a <- next_point$a

      }
  }

  if (bound_type == "lipschitz"){

    derivatives_ref <- derivatives[[1]](xi_ref)
    for (j in 2:length(derivatives)){
      derivatives_ref = derivatives_ref + derivatives[[j]](xi_ref)
    }

    for (i in 2:n){

      next_point <- next_point_lipschitz(xi_rec[i-1, ], theta_rec[i-1, ], t_flip_rec[i-1], d, derivatives, bounds, xi_ref, derivatives_ref, p, subsample)

      xi_rec[i, ] <- next_point$xi
      theta_rec[i, ] <- next_point$theta
      t_flip_rec[i] <- next_point$t_flip

    }
  }

  return(list(xi=xi_rec, theta=theta_rec, t_flip=t_flip_rec))
}

#' gen_sample
#'
#' Function which generates discrete samples from continuous sample path
#'
#' @param xi an vector of skeleton points
#' @param theta a vector of directions
#' @param t total time of the trajectory
#' @param B number of points to generate
#'
gen_sample = function(xi, theta, t_flip, B){
  sample = matrix(0, nrow = B, ncol = ncol(xi))
  for(i in 1:B){
    ind = max(which(t_flip<(i*max(t_flip)/B)))
    sample[i,] = xi[ind,] + ((i*max(t_flip)/B)-t_flip[ind])*theta[ind,]
  }
  return(sample)
}

#' integrate
#'
#' Function which estimates integral of function using discrete samples
#'
#' @param f a function to integrate
#' @param xi an vector of skeleton points
#' @param theta a vector of directions
#' @param t total time of the trajectory
#' @param B number of points to generate
#'
integrate = function(f, xi, theta, t_flip, B){
  return(mean(apply(gen_sample(xi, theta, t_flip, B), 1, f)))
}

#' ESS
#'
#' Function which estimates the effective sample size for a function h
#'
#' @param h function h taking real values
#' @param xi an vector of skeleton points
#' @param theta a vector of directions
#' @param t total time of the trajectory
#' @param B number of batches
#' @param num number of samples to use per batch
#'
ESS = function(h, xi, theta, t_flip, B, num = round(nrow(xi)/B)){
  y = numeric(B)
  tau = max(t_flip)
  samples = gen_sample(xi, theta, t_flip, B*num)
  for(i in 1:B){
    y[i] = sqrt(tau/B) * mean(apply(samples[((i-1)*num+1):(i*num),], 1, h))
  }
  var_as = var(y)
  mean_h = mean(apply(samples, 1, h))/tau
  var_h = mean(apply(samples, 1, h)^2)/tau
  return(tau*(var_h-mean_h^2)/var_as)
}


#####TESTING
###Cauchy###
del = function(x){
  return(3*x/(1+sum(x^2)))
}# max =1
set.seed(888)
test=skeleton(c(1,1),c(1,1),10000,del,bounds=c(1,1),bound_type="global")
library(graphics)
z=matrix(0,nrow=length(seq(-7.5,7.5,0.1)),ncol=length(seq(-7.5,7.5,0.1)))
for(i in 1:length(seq(-7.5,7.5,0.1))){
  for(j in 1:length(seq(-7.5,7.5,0.1))){
    z[i,j] = log(1/(1+sum(c(seq(-7.5,7.5,0.1)[i],seq(-7.5,7.5,0.1)[j])^2))^(1.5))
  }
}

par(mfrow=c(1,1),mai=c(0.3,0.7,0.1,0.1),oma=c(1,0,0,0))
filled.contour(
  x = seq(-7.5,7.5,length.out=nrow(z)),
  y = seq(-7.5,7.5,length.out=ncol(z)),
  z = z,
  levels = seq(-7.5,0,0.5),
  plot.axes={lines(test$xi,cex=0.5); axis(1, seq(-7,7,2)); axis(2, seq(-7,7,2)); box()}
)

lines(test$xi,type="l",xlab=expression(xi[1]),ylab=expression(xi[2]),xlim=c(-7,7),ylim=c(-7,7))
###Multinormal###
library(MASS)
n=1000
set.seed(888)
data = mvrnorm(n, mu=c(0,0), Sigma=matrix(c(1,0.5,0.5,2),nrow=2))
grad = function(x){
  return(x - matrix(c(8/7,-2/7,-2/7,4/7),nrow=2)%*%(apply(data,2,sum)-n*x))
}
Q = matrix(c(1000*8/7+1,-1000*2/7,-1000*2/7,1000*4/7+1),nrow=2)
test2 = skeleton(c(0,0),c(1,1),10000,grad,bounds=Q,bound_type="hessian")
plot(test2$xi,type="l",xlab=expression(xi[1]),ylab=expression(xi[2]),xlim=c(-0.1,0.2),ylim=c(-0.15,0.15))
library(mvtnorm)
post = function(x){
  return(dmvnorm(x,log=T)+sum(dmvnorm(data,mean=x,sigma=matrix(c(1,0.5,0.5,2),nrow=2),log=T)))
}
z=matrix(0,nrow=length(seq(-0.1,0.2,0.01)),ncol=length(seq(-0.15,0.15,0.01)))
for(i in 1:length(seq(-0.1,0.2,0.01))){
  for(j in 1:length(seq(-0.15,0.15,0.01))){
    z[i,j] = post(c(seq(-0.1,0.2,0.01)[i],seq(-0.15,0.15,0.01)[j]))
  }
}
filled.contour(
  x = seq(-0.1,0.2,length.out=nrow(z)),
  y = seq(-0.15,0.15,length.out=ncol(z)),
  z = z,
  levels = seq(-3120,-3080,5),
  plot.axes={lines(test2$xi,cex=0.5); axis(1, seq(-0.1,0.2,0.05)); axis(2, c(-0.15,-0.1,-0.05,0,0.05,0.1,0.15)); box()}
)
