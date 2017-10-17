#' get_samples
#' 
#' A function which produces evenly spaced samples from the skeleton of a Zig-zag process
#' 
#' @param skeleton A Zig-zag skeleton produced by "skeleton" function
#' @param k No. of samples to be drawn
#' 
#'

get_samples <- function(skeleton, k){
  
  xi <- skeleton[[1]]
  theta <- skeleton[[2]]
  times <- skeleton[[3]]
  samples <- matrix(0,k,ncol(xi))
  
  t_end <- times[length(times)]
  
  v <- 1:k
  sample_times <- (t_end/k)*v
  zz_integr
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

zz_integrate <- function(f, skeleton, averaging = "discrete", k=0, ...){
  
  if (averaging == "discrete"){
    print("test")
    samples <- get_samples(skeleton, k)
    return(mean(apply(samples,1,f,...)))
    
  }
  
  if (averaging == "continuous") {
    
  xi <- skeleton[[1]]
  theta <- skeleton[[2]]
  times <- skeleton[[3]] 
  
  l <- length(times)
    
  pi_hat <- numeric(l-1)
    
  for (i in 1:(l-1)){
    
    f_i <- function(s){
      f(xi[i, ] + (s - times[i])*theta[i, ],...)
    }
      
    lower_lim <- times[i]
    upper_lim <- times[i+1]
    pi_hat[i] <- integrate(Vectorize(f_i), lower_lim, upper_lim)$value
    }
  
  return(1/times[l]*sum(pi_hat))
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
  