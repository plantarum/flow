a2b1sc <-
  function (xx, a1, Ma, Sa, a2, b1, Mb, Sb, SCa, intensity) { 
    ## a1 == highest G1 peak
    a1/(sqrt(2 * pi) * Sa) * exp(-((xx - Ma)^2)/(2 * Sa^2)) +
      ## a2 == G2 peak for a1
      a2/(sqrt(2 * pi) * Sa * 2) * exp(-((xx - Ma * 2)^2)/(2 * (Sa * 2)^2)) +
        ## b1 == second highest peak, assumed to be the other G1 peak
        b1/(sqrt(2 * pi) * Sb) * exp(-((xx - Mb)^2)/(2 * Sb^2)) +
        ## single cut debris curve
        singleCut(SCa, intensity, xx) 
}

a2b1scInit <- function(mCall, LHS, data) {
  ## Not sure we need this fancy stuff, given we have the data already in
  ## hand, and in order:
  
  xy <- sortedXyData(mCall[["xx"]], LHS, data)
  ##xy <- data[, "x", "intensity"]
  window = 20                           # should window be a function
                                        # argument? 
  peaks <- findPeaks(xy[, "y"], window)     

  peaks <- cleanPeaks(peaks, window)

  ## ensure A is the peak with the lower mean:
  peakA <- ifelse(peaks[1, "mean"] < peaks[2, "mean"], 1, 2)
  peakB <- ifelse(peakA == 1, 2, 1)
  
  Ma <- peaks[peakA, "mean"]  
  Sa <- Ma / 20                         # assume CV = 0.05
  a1 <- peaks[peakA, "height"] * Sa / 0.4

  ## Is a2 off the chart?
  if((peaks[peakA, "mean"] * 2) > max(xy[ ,"x"]))
    a2 <- 0
  else
    a2 <- xy[peaks[peakA, "mean"] * 2, "y"] * Sa * 2 / 0.4

  Mb <- peaks[peakB, "mean"]
  Sb <- Mb / 20
  b1 <- peaks[peakB, "height"] * Sb / 0.4

  SCa <- 0.1                            # just a wild guess for now
  value <- c(a1, Ma, Sa, a2, b1, Mb, Sb, SCa)
  names(value) <-
    mCall[c("a1", "Ma", "Sa", "a2", "b1", "Mb", "Sb", "SCa")] 
  value
}

SSa2b1sc <- selfStart(a2b1sc, a2b1scInit,
                      c("a1", "Ma", "Sa", "a2", "b1", "Mb", "Sb", "SCa"))

a2b2sc <-
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

a2b2scInit <- function(mCall, LHS, data) {
  ## Not sure we need this fancy stuff, given we have the data already in
  ## hand, and in order:
  
  xy <- sortedXyData(mCall[["xx"]], LHS, data)
  ##xy <- data[, "x", "intensity"]
  window = 20                           # should window be a function
                                        # argument? 
  peaks <- findPeaks(xy[, "y"], window)     

  peaks <- cleanPeaks(peaks, window)

  Ma <- peaks[1, "mean"]  
  Sa <- Ma / 20                         # assume CV = 0.05
  a1 <- peaks[1, "height"] * Sa / 0.4

  ## Is a2 off the chart?
  if((peaks[1, "mean"] * 2) > max(xy[ ,"x"]))
    a2 <- 0
  else
    a2 <- xy[peaks[1, "mean"] * 2, "y"] * Sa * 2 / 0.4

  Mb <- peaks[2, "mean"]
  Sb <- Mb / 20
  b1 <- peaks[2, "height"] * Sb / 0.4
  if((peaks[2, "mean"] * 2) > max(xy[,"x"]))
    b2 <- 0
  else
    b2 <- as.vector(xy[peaks[2, "mean"] * 2, "y"] * Sb * 2 / 0.4)

  SCa <- 0.1                            # just a wild guess for now
  value <- c(a1, Ma, Sa, a2, b1, Mb, Sb, b2, SCa)
  names(value) <-
    mCall[c("a1", "Ma", "Sa", "a2", "b1", "Mb", "Sb", "b2", "SCa")] 
  value
}

SSa2b2sc <- selfStart(a2b2sc, a2b2scInit,
                      c("a1", "Ma", "Sa", "a2", "b1", "Mb", "Sb", "b2",
                        "SCa"))  

FSa2b2sc <- selfStart(a2b2sc, flowSS,
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

gaussDipB1  <- function(b1, Mb, Sb, xx){
  b1 / (sqrt(2 * pi) * Sb) * exp(-((xx - Mb)^2)/(2 * Sb^2))
}

