#' summary.zz
#'
#' Function which outputs summary statistics for an object of class "zz"
#'
#' @param skeleton object of class "zz"
#'
summary.zz = function(skeleton){
  cat("The number of skeleton points is ", length(skeleton$t_flip),".\n", sep="")
  cat("The effective sample size is ", ESS(identity, skeleton, sqrt(length(skeleton$t_flip))),".\n", sep="")
}

#' plot.zz
#'
#' Function which plots traceplots for an object of class "zz" in one and two dimensions
#'
#' @param skeleton object of class "zz"
#' @param contour default false; indicating whether or not to add contour lines to a 2d plot
#' @param f (log)-density of the target distribution used for contour lines
#' @param lev default 15; number of contour levels to plot
#'
plot.zz = function(skeleton, contour = FALSE, f = identity, lev = 15){

  if(ncol(skeleton$xi) == 1){
    plot(skeleton$t_flip, skeleton$xi, type = "l", xlab = "Time units", ylab = expression(xi))
  }

  if(ncol(skeleton$xi) == 2){
    if(contour == FALSE){
      plot(skeleton$xi, type="l", xlab = expression(xi[1]), ylab = expression(xi[2]))
    }
    if(contour == TRUE){
      axe_x = 1.1*seq(range(skeleton$xi[,1])[1], range(skeleton$xi[,1])[2], 0.01*(range(skeleton$xi[,1])[2]-range(skeleton$xi[,1])[1]))
      axe_y = 1.1*seq(range(skeleton$xi[,2])[1], range(skeleton$xi[,2])[2], 0.01*(range(skeleton$xi[,2])[2]-range(skeleton$xi[,2])[1]))
      z = matrix(0, nrow=101, ncol=101)
      for(i in 1:101){
        for(j in 1:101){
          z[i,j] = f(c(axe_x[i], axe_y[j]))
        }
      }
      filled.contour(
        x = axe_x,
        y = axe_y,
        z = z,
        levels = seq(min(z),max(z),(max(z)-min(z))/lev),
        plot.axes={lines(skeleton$xi,cex=0.5); axis(1); axis(2); box()}
      )
    }
  }

  if(ncol(skeleton$xi) > 2){
    cat("Cannot plot more than two dimensions (right now).")
  }

}
