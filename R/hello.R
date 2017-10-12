# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#Initializations for testing
#in this example x=(x1,x2,theta1,theta2)
del1 = function(x){
  return(max(0, x[3]*exp(-(x[1]+x[2])^2)))
}#max 1

del2 = function(x){
  return(x[4]*dcauchy(x[2]))
}#max <0.32

#intensities as list of functions
#global bounds as vector
#state includes both xsi and theta
poispro = function(state, intensities, bounds){
  d = length(intensities)
  flip = FALSE
  arrivals = rexp(d, bounds)
  tau = min(arrivals)
  ind = which.min(arrivals)
  newstate = state
  newstate[1:d] = state[1:d] + tau*state[(d+1):(2*d)]
  if(runif(1)<intensities[[ind]](newstate)/bounds[ind]){
    newstate[d+ind] = -newstate[d+ind]
    flip = TRUE
  }
  return(list(state=newstate, time=tau, index=ind, flip=flip))
}
