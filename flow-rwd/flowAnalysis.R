## Functions for analyzing flowHist datasets

fhAnalyze <- function(fh){
  fh$nls <- fhCall(fh, "nls")
  fh$counts <- fhCount(fh)
  fh$cv <- fhCV(fh)
  fh
}

fhCall <- function(fh, fun){
  on.exit(detach(list(fHModel = fh$model)))
  if (any(grepl("list", search()))) stop("list already on search path")
  attach(list(fHModel = fh$model))
  form1 <- paste("intensity ~ fHModel(")
  args <- as.character(names(formals(fh$model)))
  args <- args[!args %in% c("", "intensity", "xx")]
  args <- paste(args, collapse = ", ")
  form3 <- ", intensity = intensity, xx = x)"
  form <- as.formula(paste(form1, args, form3), env = globalenv())

  eval(call(fun, form, data = fh$data))
}  

fhnls <- function(fh){
  ## the port algorithm with lower = 0 seems to run into infinity values
  ## problems 
  on.exit(detach(list(fHModel = fh$model)))
  if (any(grepl("list", search()))) stop("list already on search path")
  attach(list(fHModel = fh$model))
  form1 <- paste("intensity ~ fHModel(")
  args <- as.character(names(formals(fh$model)))
  args <- args[!args %in% c("", "intensity", "xx")]
  args <- paste(args, collapse = ", ")
  form3 <- ", intensity = intensity, xx = x)"
  form <- as.formula(paste(form1, args, form3), env = globalenv())

  eval(call("nls", form, data = fh$data,
            control = nls.control(minFactor = 1/2048)))  
}  

fhCount <- function(fh, lower = 0, upper = 256, subdivisions = 1000){
  total <-
    do.call(integrate,
            c(substitute(fh$model),
              as.list(coef(fh$nls)),
              intensity = substitute(fh$data$intensity),
              lower = lower, upper = upper,
              subdivisions = subdivisions))
  firstPeak <-
    integrate(fA1, a1 = coef(fh$nls)["a1"],
              Ma = coef(fh$nls)["Ma"],
              Sa = coef(fh$nls)["Sa"],
              lower = lower, upper = upper,
              subdivisions = 1000)
  secondPeak <-
    integrate(fB1, b1 = coef(fh$nls)["b1"],
              Mb = coef(fh$nls)["Mb"],
              Sb = coef(fh$nls)["Sb"],
              lower = lower, upper = upper,
              subdivisions = 1000)
              
  return(list(total = total, firstPeak = firstPeak,
              secondPeak = secondPeak)) 
}  

fhCV <- function(fh){
  CVa <- coef(fh$nls)["Sa"]/coef(fh$nls)["Ma"]
  CVb <- coef(fh$nls)["Sb"]/coef(fh$nls)["Mb"]
  CI <- deltaMethod(fh$nls, "Ma/Mb")
  return(list(CVa = CVa, CVb = CVb, CI = CI))
}

