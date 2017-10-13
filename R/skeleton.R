#' next_point_global
#'
#' a function which proposes a new skeleton point for the Zig-zag process given global bounds on the directional
#' derivatives
#'
#'@param xi current skeleton position
#'@param theta current direction
#'@param t_flip most recent flip time
#'@param d number of spatial dimensions
#'@param dervatives function which takes position x as argument and outputs evaluation of directional derivatives at x
#'@param bounds vector of global bounds for directional derivatives
#'
#'

# unclear how much this function should include, see somments below in "skeleton"
next_point_global <- function(xi, theta, t_flip, d, derivatives, bounds){

    flip <- FALSE

    while(flip==FALSE){

      test_taus <- rexp(d, bounds)
      tau <- min(test_taus)
      ind <- which.min(test_taus)

      xi <- xi + tau*theta
      t_flip <- t_flip + tau

      if(runif(1) < max(0, theta[ind]*derivatives(xi)[ind])/bounds[ind]){
        theta[ind] <- -theta[ind]
        flip = TRUE
      }

    }

    return(list(xi = xi, theta = theta, t_flip = t_flip))
}

#' next_point_hess
#'
#'@param xi current skeleton position
#'@param theta current direction
#'@param t_flip most recent flip time
#'@param d number of spatial dimensions
#'@param dervatives function which takes position x as argument and outputs evaluation of directional derivatives at x
#'@param a parameter defining hessian-type bound, calculated elsewhere
#'@param a parameter defining hessian-type bound, calculated elsewhere
#'
#'
next_point_hess <- function(xi, theta, t_flip, d, derivatives, a, b){

  flip <- FALSE

  while(flip==FALSE){

    test_taus <- -a/b + sqrt((a/b)^2 + 2/b*(-log(runif(d))))
    tau <- min(arrivals)
    ind <- which.min(arrivals)

    a <- a + b*tau
    xi <- xi + tau*theta
    t_flip <- t_flip + tau

    a_temp <- theta[ind]*derivatives(xi)[ind]

    if(runif(1) < max(0, a_temp)/max(0, a[ind])){
      theta[ind] = -theta[ind]
      flip = TRUE
    }

    a[ind] = a_temp*(-1)^as.numeric(flip)

  }

  return(list(xi=xi, theta = theta, t_flip <- t_flip, flip=flip, a=a))
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
#'
skeleton <- function(xi, theta, n, derivatives, bounds, bound_type = "global"){

  #set inital time, no. dimensions and record tables for skeleton
  t_flip <- 0
  d <- length(xi)

  xi_rec <- matrix(0,n,d)
  theta_rec <- matrix(0,n,d)
  t_flip_rec <- numeric(n)

  xi_rec[1, ] <- xi; theta_rec[1, ] <- theta

  #start flip counter
  points <- 2

  #order of while & if loops here should be figured out - pointless to check bound_type at every while iteration,
  #but want to do if(test_point$flip == TRUE) part for every bound type, so also pointless to type this many times -
  #maybe arranging next_point functions better would help here

  #In response: Check whether new point is actually flipped in next_point function so that these actually give the next
  #point we want

  if (bound_type == "global"){

    for (i in 2:n){

    next_point <- next_point_global(xi_rec[i-1,], theta_rec[i-1,], t_flip[i-1], d, derivatives, bounds)

        xi_rec[i, ] <- next_point$xi
        theta_rec[i, ] <- next_point$theta
        t_flip_rec[i] <- next_point$t_flip

    }
  }


  if (bound_type == "hessian"){

    a = numeric(d); b = numeric(d)
      for(i in 1:d){
          a[i] = theta[i]*derivatives(xi)[i]
          b[i] = sqrt(d)*sqrt(sum(bounds[,i]^2))
      }

    for (i in 2:n){

      next_point <- next_point_hess(xi, theta, t_flip, d, derivatives, a, b)

        xi_rec[i, ] <- next_point$xi
        theta_rec[i, ] <- next_point$theta
        t_flip_rec[i] <- next_point$t_flip
        a <- next_point$a

      }
    }

  }
  return(list(xi=xi_rec, theta=theta_rec, t_flip=t_flip_rec))
}

#####TESTING
del = function(x){
  return(2*x/(1+x^2))
}# max =1

test=skeleton(c(1,1),c(1,1),10000,del,bounds=c(1,1))
