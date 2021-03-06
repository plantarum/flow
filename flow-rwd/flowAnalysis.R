## Functions for analyzing flowHist datasets

fhAnalyze <- function(fh){
  fh$nls <- fhNLS(fh)
  fh$counts <- fhCount(fh)
  fh$cv <- fhCV(fh)
  fh$RCS <- fhRCS(fh)
  fh
}

fhRCS <- function(fh){
  #########################################################################
  ## This may not be a useful measure of analysis quality. It is heavily ##
  ## influenced by the highest channels, where the expected value is     ##
  ## close to zero. Consequently, observing 2 or 3 stray events in one   ##
  ## of these channels, where the expected value may be < 0.02, produces ##
  ## a higher value of (obs - exp)^2 / exp than larger absolute          ##
  ## differences in the main region of the histogram.                    ##
  ##                                                                     ##
  ## Rabinovitch 1994:                                                   ##
  ##                                                                     ##
  ## The x2 is affected by a large number of variables, not all related  ##
  ## to goodness of the fit; these include the number of cells acquired  ##
  ## in the histogram and the end points of the analysis region used     ##
  ## within the histogram.                                               ##
  #########################################################################

  obs <- fh$data$intensity
  exp <- predict(fh$nls)
  zeros <- zapsmall(exp, digits = 5) == 0
  obs <- obs[!zeros]
  exp <- exp[!zeros]
  chi <- sum(((obs - exp)^2) / exp)

  n <- length(obs)
  m <- length(formals(fh$model)) - 2    # don't count xx or intensity as
                                        # parameters

  chi/summary(fh$nls)$df[2]
}
                  
                  
fhNLS <- function(fh){
  model <- fh$model
  form1 <- paste("intensity ~ model(")
  args <- as.character(names(formals(fh$model)))
  args <- args[!args %in% c("", "intensity", "xx")]
  args <- paste(args, collapse = ", ")
  form3 <- ", intensity = intensity, xx = x)"
  form <- as.formula(paste(form1, args, form3))

  eval(call("nls", form, start = fh$init, data = fh$data)) 
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

