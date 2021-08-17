
#This file contains a function to simplify model selection procedure 
#marginal AIC(mAIC) is used as model selection measure

require("optimParallel")
# Refit the model, possibly with some modifications to map argument using previous
# estimates as starting values. 
# 
# Arguments:
#     model: A list representing a previously fitted model with components
#     estimates:  Estimates from a previous model fit in same format as parameters
#     parameters: Parameter values of restricted models
#     map:        Specifies which parameters were estimated or restricted (input the MakeADFun)
#     args.MakeADFun: Additional arguments to MakeADFun (excluding map and parameters)
#     args.nlminb: Additional argument used to fit the model
#     sdreport:   Standard errors etc...
#   changemap:  A list containing changes we want in the map-argumentn
# 
# Value: 
#    A list with the same format as the model argument containing the new model


subsetpar <- function(x,name)
  subset(x,names(x)==name)

update.tmbmodel <- function(model,optimizer,changemap=NULL,cv=NULL, seed=NULL, bestmodel=FALSE) {
  newmodel <- model
  
  # Compute map for newmodel
  for (parname in names(changemap)) {
    newmodel$map[[parname]] <- changemap[[parname]]
  }
  
  # Compute starting values for newmodel in part using estimates from model
  start <- model$parameters # use parameter values under restricted model as default starting value
  if (!is.null(model$estimates)) { # if previous estimates are available
    for (parname in names(start)) {
      if (parname %in% names(newmodel$map)) {
        estimated.subset <- !is.na(newmodel$map[[parname]])
        # then for elements that we want to estimate, use previous estimate as starting value
        start[[parname]][estimated.subset] <- model$estimates[[parname]][estimated.subset]  
      } else {
        start[[parname]] <- model$estimates[[parname]]
      }
    }
  }
  
  # Make new objective function
  args.MakeADFun <- model$args.MakeADFun
  args.MakeADFun$parameters <- start 
  args.MakeADFun$map <- newmodel$map
  obj <- do.call("MakeADFun", args.MakeADFun)
  
  
  # Add coefficient of variation (cv) of lognormal changes in all starting values
  # seed is to make sure the same random value generated and make the model-fit reproducable
  for (varname in names(newmodel$args.optim))
    obj[[varname]] <- model$args.optim[[varname]]
  if (!is.null(cv)) {
    if (is.null(seed)) {
      oldseed <- .Random.seed
      set.seed(seed)
    }
    obj$par <- obj$par*rlnorm(length(obj$par),sdlog=cv)
    if (is.null(seed)) 
      assign(".Random.seed", oldseed, .GlobalEnv)
  }
  
  # call optimizer
  if (optimizer=="nlminb"){
    opt <- nlminb(obj$par, obj$fn, obj$gr,iter.max=2000)
    newmodel$maic <-  2*length(opt$par) - 2*(-opt$objective)} # 2p - 2l_max
  else if(optimizer=="optim"){
    opt <- optim(obj$par, obj$fn, obj$gr,method="BFGS",control=list(maxit=2000,reltol=1e-12))
    newmodel$maic <-  2*length(opt$par) - 2*(-opt$value)} # 2p - 2l_max
  
  if (bestmodel) {newmodel$sdreport <- sdreport(obj, getReportCovariance = TRUE)}
  else{newmodel$sdreport <- sdreport(obj, getReportCovariance = FALSE)}
  
  newmodel$summary <- summary(newmodel$sdreport,"fixed",p.value = TRUE)
  newmodel$opt <- opt
  #newmodel$maic <-  2*length(opt$par) - 2*(-opt$value) # 2p - 2l_max
  #newmodel$maic <-  2*length(opt$par) - 2*(-opt$objective) # 2p - 2l_max
  newmodel$obj <- obj
  newmodel$convergence  <- opt$convergence
  newmodel$message <- opt$message
  newmodel$mgc <- obj$gr(opt$par)
  #hess <- hessian(obj$fn,newmodel$sdreport$par.fixed)
  #hess[is.na(hess)] <- 0
  #newmodel$eigenvalues <- eigen(hess)$values
  newmodel$noPara<- length(opt$par) 
  
  #Put all estimated and constrained parameters into estimates
  newmodel$estimates <- newmodel$parameters 
  # the initial values for Phi and rho are kept as 0
  for (parname in names(newmodel$parameters)[! names(newmodel$parameters) %in% c("Phi","parCor" )]) {
    if (parname %in% names(newmodel$map)) { 
      estimated.subset <- !is.na(newmodel$map[[parname]])
      newmodel$estimates[[parname]][estimated.subset] <- subsetpar(opt$par,parname)[newmodel$map[[parname]]][estimated.subset]
    } else { # parameter not in map argument 
      if (parname %in% names(opt$par)) { # then if is estimated extract those from opt$par
        newmodel$estimates[[parname]] <- subsetpar(opt$par,parname) 
        names(newmodel$estimates[[parname]]) <- NULL
      } else { # this must be a random effect
        newmodel$estimates[[parname]] <- subsetpar(summary(newmodel$sdreport,"random")[,1], parname) 
        names(newmodel$estimates[[parname]]) <- NULL
        dim(newmodel$estimates[[parname]]) <- dim(newmodel$parameters[[parname]])
      }
    }
  }
  
  newmodel
}



print.tmbmodel <- function(model) {
  #print(summary(model$sdreport,"fixed"))
  #print(summary(model$sdreport,"report"))
  gradient <- model$mgc
  names(gradient) <- names(model$opt$par)
  cat("Gradient\n")
  print(gradient)
  #cat("eigenvalues of hessian matrix\n")
  print(model$eigenvalues)
  cat("convergence = ",model$convergence,"\n")
  cat("maic = ",model$maic,"\n")
  cat("number of paramters = ",model$noPara,"\n")
  #print(model$summary)
}





