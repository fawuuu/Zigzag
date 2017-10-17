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



