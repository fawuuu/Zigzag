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
  derivatives_baseline = theta * derivatives_ref
  derivatives_baseline[which(derivatives_baseline<0)] = 0

  if(p == Inf){
    a <- derivatives_baseline + bounds*max(abs(xi-xi0))
    b <- bounds
  }else{
    a <- derivatives_baseline + bounds*sum(abs(xi-xi0)^p)^(1/p)
    b <- bounds*d^(1/p)
  }

  while(flip == FALSE){

    test_taus <- -a/b + sqrt((a/b)^2 + 2/b*(-log(runif(d))))

    tau <- min(test_taus)
    ind <- which.min(test_taus)
    
    xi <- xi + tau*theta
    t_flip <- t_flip + tau
    
    if(subsample == FALSE){
      
      if(runif(1) < max(0, theta[ind] * derivatives(xi)[ind])/max(0, a[ind]+b*tau)){
        theta[ind] <- -theta[ind]
        flip <- TRUE
      }
    }else{
      j <- sample(1, 1:length(derivatives))
      e_ij = derivatives_ref[ind] + derivatives[[j]](xi)[ind] - derivatives[[j]](xi0)[ind]
      if(runif(1) < max(0, theta[ind] * e_ij)/max(0, a[ind]+b*tau)){
        theta[ind] <- -theta[ind]
        flip <- TRUE
        
      }
      
      
      if(p == Inf){
        a = derivatives_baseline + bounds*max(abs(xi-xi0))
      }else{
        a = derivatives_baseline + bounds*sum(abs(xi-xi0)^p)^(1/p)
      }
      
    }
    
  }
  
  return(list(xi = xi, theta = theta, t_flip = t_flip))
}