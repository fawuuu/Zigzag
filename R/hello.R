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
del = function(x){
  return(2*x/(1+x^2))
}# max =1rep(1,)

#derivatives as function with vector output
#global bounds as vector
#state includes both xsi and theta
poispro_global = function(state, derivatives, bounds){
  d = length(state)/2
  flip = FALSE
  arrivals = rexp(d, bounds)
  print(arrivals)
  tau = min(arrivals)
  ind = which.min(arrivals)
  newstate = state
  newstate[1:d] = state[1:d] + tau*state[(d+1):(2*d)]
  if(runif(1)<max(0, state[d+ind]*derivatives(newstate[1:d])[ind])/bounds[ind]){
    newstate[d+ind] = -newstate[d+ind]
    flip = TRUE
  }
  return(list(state=newstate, time=tau, index=ind, flip=flip))
}

#this should be outside of function call
a = numeric(d); b = numeric(d)
for(i in 1:d){
  a[i] = state[d+i]*derivatives(state[1:d])[i]
  b[i] = sqrt(d)*sqrt(sum(hessian[,i]^2))
}

poispro_hessian = function(state, derivatives, a, b){
  d = length(state)/2
  flip = FALSE
  arrivals =  -a/b + sqrt((a/b)^2 + 2/b*(-log(runif(d))))
  tau = min(arrivals)
  ind = which.min(arrivals)
  a = a + b*tau
  newstate = state
  newstate[1:d] = state[1:d] + tau*state[(d+1):(2*d)]
  a_temp = state[d+ind]*derivatives(newstate[1:d])[ind]
  print(max(0, a_temp)/max(0, a[ind]))
  if(runif(1)<max(0, a_temp)/max(0, a[ind])){
    newstate[d+ind] = -newstate[d+ind]
    flip = TRUE
  }
  a[ind] = a_temp*(-1)^as.numeric(flip)
  return(list(state=newstate, time=tau, index=ind, flip=flip, a=a))
}

poispro_subsamp = function(state, derivatives){

}
