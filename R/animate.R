plot_to_time <- function(skeleton, t){
  
  times <- skeleton$t_flip
  ind <- max(which(times < t))
  
  end_point <- skeleton$xi[ind, ] + (t - times[ind]) * skeleton$theta[ind, ]
  
  plot(skeleton$xi[1:ind, ], type = "l")
  points(end_point[1],end_point[2], col="red", pch = 21, bg ="red")
  segments(skeleton$xi[ind,1], skeleton$xi[ind,2], end_point[1], end_point[2])
  
}

animate <- function(skeleton){
  
  saveHTML({ani.options(interval = 0.05, nmax = 50)
  
    t <- seq(from = 1/200, to = 1, by = 1/200)
  
    for (i in seq_along(t)){
      
      print(t[i])
    
      plot_to_time(skeleton,t[i])
    
      ani.pause()}
    
    }, htmlfile = "Zigzag3.html", ani.width = 600, ani.height = 600)

}