#' get_samples
#' 
#' A function which produces evenly spaced samples from the skeleton of a Zig-zag process
#' 
#' @param skeleton A Zig-zag skeleton produced by "skeleton" function
#' @param k No. of samples to be drawn
#' 
#'

get_samples <- function(skeleton, k){
  
  xi <- skeleton$xi
  theta <- skeleton$theta
  times <- skeleton$t_flip
  samples <- matrix(0,k,ncol(xi))
  
  t_end <- times[length(times)]
  
  v <- 1:k
  sample_times <- (t_end/k)*v
  
  j <- 1
  for (i in seq_along(sample_times)){
    
    while(sample_times[i] < times[j] | sample_times[i] > times[j+1]){
      j <- j+1
    }
    samples[i, ] <- xi[j, ] + (sample_times[i] - times[j])*theta[j, ]
  }
  samples
}

#' zz_integrate
#' 
#' A function which performs Monte Carlo estimates of integrals wrt a distribution pi
#' 
#' @param f function whose integral is to be estimated
#' @param skeleton skeleton of Zig-zag process with stationary distribution pi 
#' @param averaging "discrete" or "continuous" - specifies type of ergoidc average to take
#' @param k number of samples to use in discrete estimation
#' 
#' 

zz_integrate <- function(f, skeleton, averaging = "discrete", k, ...){
  
  if (averaging == "discrete"){
    
    samples <- get_samples(skeleton, k)
    f_samples <- apply(samples, 1, f)
    if(class(f_samples) == "matrix"){
      return(apply(f_samples, 1, mean))
    }else{
      return(mean(f_samples))
    } 
    
  }
  
  if (averaging == "continuous") {
    
  xi <- skeleton$xi
  theta <- skeleton$theta
  times <- skeleton$t_flip
  
  l <- length(times)
    
  pi_hat <- matrix(0, nrow = l-1, ncol = length(f(xi[1,])))
  
  for(j in 1:length(f(xi[1,]))){
    
    for (i in 1:(l-1)){
      
      f_i <- function(s){
        f(xi[i, ] + (s - times[i])*theta[i, ],...)[j]
      }
      
      lower_lim <- times[i]
      upper_lim <- times[i+1]
      pi_hat[i,j] <- integrate(Vectorize(f_i), lower_lim, upper_lim)$value
    }
    
  }
  
  return(1/times[l]*apply(pi_hat,2,sum))
  }
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
ESS <- function(h, skeleton, B, num = round(nrow(skeleton$xi)/B)){
  
  y <- numeric(B)
  tau <- max(skeleton$t_flip)
  
  samples <- get_samples(skeleton, B*num)

  for(i in 1:B){
    y[i] <- sqrt(tau/B) * mean(apply(as.matrix(samples[((i-1)*num+1):(i*num),]), 1, h))
  }
 
  var_as <- var(y)
  mean_h <- mean(apply(samples, 1, h))/tau
  var_h <- mean(apply(samples, 1, h)^2)/tau
  return(tau*(var_h-mean_h^2)/var_as)
}
  