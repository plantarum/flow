library(caTools)                        # for runMax

## Single-cut debris model
## 
## S(x) = a ⱼ₌ₓ₊₁∑ⁿ ³√j YⱼPₛ(j, x)
## Pₛ(j, x) = 2 / (πj √(x/j)(1 - x/j)
##
## a = amplitude parameter
## Yⱼ = intensity in channel j
## Pₛ(j, x) = probability of a nuclei from channel j falling into channel x
## when cut.
##
## By this formula, the debris intensity in a channel/bin is a function of
## the intensity in all the subsequent bins. This recursive relationship is
## tricky to code; in order to take full advantage of all of the R tools
## that support nls, the mean function needs to return one fitted value for
## one predictor value. The following implementation of singleCut therefore
## takes the entire vector of the response vector (intensity), necessary to
## calculate the debris curve, and returns only the value for a single
## predictor value.

singleCutBase <- function(SCa, intensity, xx){
  # Don't extend the model below/beyond the data
  # Modfit appears to cut off the debris slightly above the lowest data,
  # which gives a better fit. Perhaps set first.channel to 2-4? Need to
  # test this and determine best fit. Possibly use an extra parameter to
  # tune this for each data set individually.
  first.channel <- which(fl3.df$intensity > 0)[2]

  res <- 0
  if(xx >= first.channel & xx < length(intensity)){
    channels = (xx + 1):length(intensity)
    for(j in channels){
      res <- res + j^(1/3) * intensity[j] * 2 / (pi * j * sqrt(xx/j * (1 - xx/j)))
    }
  }
  return(SCa * res)
}

## Need to Vectorize this so our comparisons in line 2 make sense. Without
## this, the comparisons try to test all values of the xx vector against
## first.channel and length(intensity). This raises a warning, and the
## results are not what we want. We need them tested one at a time, hence
## the vectorize here: 
singleCut <- Vectorize(singleCutBase, "xx")

findPeaks <- function(dat, window){
  ## extract all peaks from data
  localMax <- runmax(dat, k = window)
  isMax <- localMax == dat
  maxVals <- dat[isMax]
  cbind(mean = (1:length(dat))[isMax], height = maxVals)
}

cleanPeaks <- function(peaks, num, window){
  ## Pick out the top num peaks for use as starting values

  ## Screen out any ties - if two peaks have the same height, and are
  ## within the same 'window', we need to drop one.
  
  ## If a peak has a 'match' at half the size, use the smaller peak (ie.,
  ## take the G1 peak in cases where the G2 peak is higher) 

  ## After the first peak is selected, only consider peaks that are not a
  ## multiple of the size of this peak when selecting the next one.

  peaks <- peaks[order(peaks[,2], decreasing = TRUE), ]

  drop <- numeric()
  for(i in 2: nrow(peaks)){
    if((peaks[i-1, "height"] == peaks[i, "height"]) &
       (abs(peaks[i-1, "mean"] - peaks[i, "mean"]) <= window)){ 
      ## It's a tie!
      drop <- c(drop, i)
    }
  }

  peaks <- peaks[-drop, ]

  out <- matrix(NA, nrow = num, ncol = 2)

  numc <- num
  while(nrow(peaks) > 0 & numc > 0){
    for(i in 1:(nrow(peaks) - 1)){
      if(peaks[1, "mean"] ... 
      
  
}

oneSampleSC <-
  function (xx, a1, Ma, Sa, a2, b1, Mb, Sb, b2, SCa, intensity) { 
    ## a1 == highest G1 peak
    a1/(sqrt(2 * pi) * Sa) * exp(-((xx - Ma)^2)/(2 * Sa^2)) +
      ## a2 == G2 peak for a1
      a2/(sqrt(2 * pi) * Sa * 2) * exp(-((xx - Ma * 2)^2)/(2 * (Sa * 2)^2)) +
        ## b1 == second highest peak, assumed to be the other G1 peak
        b1/(sqrt(2 * pi) * Sb) * exp(-((xx - Mb)^2)/(2 * Sb^2)) +
        ## b2 == G2 peak for b1
        b2/(sqrt(2 * pi) * Sb * 2) * exp(-((xx - Mb * 2)^2)/(2 * (Sb * 2)^2)) +
        ## single cut debris curve
        singleCut(SCa, intensity, xx) 
}

oneSampleSCInit <- function(mCall, LHS, data) {
  ## Not sure we need this fancy stuff, given we have the data already in
  ## hand, and in order:
  
  xy <- sortedXyData(mCall[["xx"]], LHS, data)
  ##xy <- data[, "x", "intensity"]
  window = 20                           # should window be a function
                                        # argument? 
  peaks <- findPeaks(xy[, "y"], window)     

  peaks <- peaks[order(peaks[, "height"], decreasing = TRUE),]

  ## Is the first peak a G2?
  ## - start with the fifth peak and work up
  ## - if there are multiple peaks in the range, we keep the highest one.
  peakOne <- NULL
  for(i in 5:2){                       
    if((peaks[i, "mean"] < 0.55 * peaks[1, "mean"]) &
       (peaks[i, "mean"] > 0.45 * peaks[1, "mean"])){
      peakOne <- i
    }
  }

  if(is.null(peakOne)){
    peakOne <- 1
  }
  
  Ma <- peaks[peakOne, "mean"]  
  Sa <- Ma / 20                         # assume CV = 0.05
  a1 <- peaks[peakOne, "height"] * Sa / 0.4

  ## Is a2 off the chart?
  if((peaks[peakOne, "mean"] * 2) > max(xy[1, ]))
    a2 <- 0
  else
    a2 <- xy[peaks[peakOne, "mean"] * 2, "intensity"] * Sa * 2 / 0.4

  peakTwo <- 
  
  ## Are the two highest peaks the G1 and G2 from the same sample?
  if((peaks[peakTwo, "mean"] < 2.1 * peaks[peakOne, "mean"] &
      peaks[peakTwo, "mean"] > 1.9 * peaks[peakOne, "mean"]) |
     (peaks[peakOne, "mean"] < 2.1 * peaks[peakTwo, "mean"] &
      peaks[peakOne, "mean"] > 1.9 * peaks[peakTwo, "mean"]))
    peaks <- peaks[-peakTwo, ]

  ## Are the two highest peaks tied values from the same peak?
  while((abs(peaks[2, "mean"] - peaks[1, "mean"]) <  window ) &
        (peaks[1, "height"] == peaks[2, "height"]))
    peaks <- peaks[-2, ]
  
  Mb <- peaks[2, "mean"]
  Sb <- Mb / 20
  b1 <- peaks[2, "height"] * Sb / 0.4
  if((peaks[2, "mean"] * 2) > max(xy[1, ]))
    b2 <- 0
  else
    b2 <- as.vector(xy[peaks[2, "mean"] * 2, "intensity"] * Sb * 2 / 0.4)

  SCa <- 0.1                            # just a wild guess for now
  value <- c(a1, Ma, Sa, a2, b1, Mb, Sb, b2, SCa)
  names(value) <-
    mCall[c("a1", "Ma", "Sa", "a2", "b1", "Mb", "Sb", "b2", "SCa")] 
  value
}

SSoneSSC <- selfStart(oneSampleSC, oneSampleSCInit,
                      c("a1", "Ma", "Sa", "a2", "b1", "Mb", "Sb", "b2",
                        "SCa"))  

gaussDipA  <- function(a1, Ma, Sa, a2, xx){
    a1 / (sqrt(2 * pi) * Sa) * exp(-((xx - Ma)^2)/(2 * Sa^2)) +
      a2 / (sqrt(2 * pi) * Sa * 2) * exp(-((xx - Ma * 2)^2)/(2 * (Sa * 2)^2))
}

peakA <- function(a1, Ma, Sa, xx){
    a1 / (sqrt(2 * pi) * Sa) * exp(-((xx - Ma)^2)/(2 * Sa^2)) 
}

peakB <- function(b1, Mb, Sb, xx){
    b1 / (sqrt(2 * pi) * Sb) * exp(-((xx - Mb)^2)/(2 * Sb^2)) 
}

gaussDipB  <- function(b1, Mb, Sb, b2, xx){
  b1 / (sqrt(2 * pi) * Sb) * exp(-((xx - Mb)^2)/(2 * Sb^2)) +
    b2 / (sqrt(2 * pi) * Sb * 2) * exp(-((xx - Mb * 2)^2)/(2 * (Sb * 2)^2))
}

nlsIntegrate <- function(nlsObj, nlsMod, fcsINT, lower = 0, upper =
                         256, subdivisions = 1000){
  do.call(integrate,
          c(substitute(nlsObj),
            as.list(coef(nlsMod)),
            intensity = substitute(fcsINT),
            lower = lower, upper = upper,
            subdivisions = subdivisions))
}

fhIntegrate <- function(fh, lower = 0, upper = 256, subdivisions = 1000){
  total <-
    do.call(integrate,
            c(substitute(fh$model),
              as.list(coef(fh$nls)),
              intensity = substitute(fh$data$intensity),
              lower = lower, upper = upper,
              subdivisions = subdivisions))
  firstPeak <-
    integrate(peakA, a1 = coef(fh$nls)["a1"],
              Ma = coef(fh$nls)["Ma"],
              Sa = coef(fh$nls)["Sa"],
              lower = lower, upper = upper,
              subdivisions = 1000)
  secondPeak <-
    integrate(peakB, b1 = coef(fh$nls)["b1"],
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

nlsIntegrate(SSoneSSC, m3test, fl3.df$intensity)

flowHist <- function(FCS, CHANNEL, MODEL= NULL, bins = 256, maxBins = 1024){
  ## You probably want to subdivide the bins evenly. i.e., if there are
  ## 1024 bins in the data, use 128, 256, 512 bins
  if((1024 %% bins) != 0)
    warning("maxBins is not a multiple of bins!")
  ## Extract the data channel
  chanDat <- exprs(FCS)[, CHANNEL]

  ## remove the top bin - this contains clipped values representing all
  ## out-of-range data, not true values
  chanTrim <- chanDat[chanDat < max(chanDat)]

  ## aggregate bins: combine maxBins into bins via hist
  binAg <- floor(maxBins / bins)

  histBins <- hist(chanTrim, breaks = seq(from = 0, to = 1024, by = binAg),
                   plot = FALSE)

  intensity <- histBins$counts

  res <- list(channel = CHANNEL,
              data = data.frame(x = 1:length(intensity),
                                intensity = intensity),
              model = MODEL, nls = NULL, standard = NULL,
              file = FCS@description$GUID)
  class(res) <- "flowHist"
  
  return(res)
}
              
print.flowHist <- function(self){
  message("flowHist object")
  message("Source file: ", self$file)
  message(paste("Channel: ", self$channel))
  message(paste("Values: ", dim(self$data)[1]))
  message(paste("Total events: ", sum(self$data$intensity)))

  ## if(is.null(self$model)){
  ##   message("No model selected")
  ## } else {
  ##   message(paste("Model: ", substitute(self$model)))
  ## }

  if(is.null(self$nls)){
    message("Not fit")
  } else {
    ## replace this with some measure of goodness-of-fit
    message(paste("Fit!"))
  }

  if(is.null(self$standard)){
    message("Standard not specified")
  } else {
    message(paste("Standard GC value: ", substitute(self$standard)))
  }

  if(!is.null(self$counts)){
    message(paste("Modelled events:", round(self$counts$total$value, 1)))
    message(paste("Peak A:", round(self$counts$firstPeak$value, 1)))
    message(paste("Peak B:", round(self$counts$secondPeak$value, 1)))
  }

  if(!is.null(self$cv)){
    message(paste("CV A:", round(self$cv$CVa, 3)))
    message(paste("CV B:", round(self$cv$CVb, 3)))
  }

}

fHnls <- function(fh){
  nls(intensity ~ SSoneSSC(a1, Ma, Sa, a2, b1, Mb, Sb, b2, SCa,
                           intensity = intensity, xx = x),
      data = fh$data)
}

fHinitial <- function(fh){
  getInitial(intensity ~ SSoneSSC(a1, Ma, Sa, a2, b1, Mb, Sb, b2, SCa,
                                  intensity = intensity, xx = x),
             data = fh$data)
}


plot.flowHist <- function(self, init = FALSE){
  plot(self$data$intensity, type = 'n', main = self$file)
  polygon(x = c(self$data$x, max(self$data$x) + 1), y = c(self$data$intensity, 0),
          col = "lightgray", border = NA)

  if(init){
    iv <- fHinitial(self)
    lines(x = self$data$x,
          y = SSoneSSC(xx = self$data$x, intensity = self$data$intensity,
                       a1 = iv["a1"], Ma = iv["Ma"], Sa = iv["Sa"],
                       a2 = iv["a2"], b1 = iv["b1"], Mb = iv["Mb"],
                       Sb = iv["Sb"], b2 = iv["b2"], SCa = iv["SCa"]),
          col = 1, lwd = 3, lty = 5)
  }
  if(! is.null(self$nls)){
    lines(x = self$data$x, y = predict(self$nls), col = 2)

    ## Debris:
    lines(x = self$data$x,
          y = singleCut(SCa = coef(self$nls)["SCa"],
                        intensity = self$data$intensity, xx = self$data$x),
          col = 3, lty = 2, lwd = 2)
    ##xx <- self$data$x                          # awkward!!

    ## Internal Standard:
    lines(x = self$data$x,
          y = gaussDipA(Ma = coef(self$nls)["Ma"], a1 = coef(self$nls)["a1"],
                        Sa = coef(self$nls)["Sa"], a2 = coef(self$nls)["a2"],
                        xx = self$data$x),
          col = 4)

    ## Sample:
    lines(x = self$data$x,
          y = gaussDipB(Mb = coef(self$nls)["Mb"], b1 = coef(self$nls)["b1"],
                        Sb = coef(self$nls)["Sb"], b2 = coef(self$nls)["b2"],
                        xx = self$data$x),
          col = 5)

  }
}


