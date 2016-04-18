library(caTools)                        # for runmax, runmean

fhPeakPlot <- function(fh, window, smooth = window/2){
  dat <- fh$data[, "intensity"]
  smDat <- runmean(dat, k = floor(smooth), endrule = "mean")
  localMax <- runmax(smDat, k = window)
  isMax <- localMax == smDat
  maxVals <- dat[isMax]                 # use the raw data for heights 
  res <- cbind(mean = (1:length(dat))[isMax], height = maxVals)

  clean <- cleanPeaks(res, window)
  plot(fh$data$intensity, type = 'l')
  points(smDat, type = 'l', col = 'grey')
  points(res, cex = 2, col = 2)
  points(clean, cex = 2, col = 3)
}

findPeaks <- function(fh, window, smooth = window / 2){
  ## Note that this function can be called from the self-starting function,
  ## or from flowHist, each with a different object class for fh. This
  ## needs attention.
  
  ## extract all peaks from data
  ## smoothing removes most of the noisy peaks
  if(class(fh) %in% c("list", "flowHist")){
    dat <- fh$data[, "intensity"]
  } else {
    dat <- fh
  }
  smDat <- runmean(dat, k = floor(smooth), endrule = "mean")
  localMax <- runmax(smDat, k = window)
  isMax <- localMax == smDat
  maxVals <- dat[isMax]                 # use the raw data for heights 
  res <- cbind(mean = (1:length(dat))[isMax], height = maxVals)
  res
}

cleanPeaks <- function(peaks, window){
  ## Remove ties and multiple peaks for histogram analysis

  ## Screen out any ties - if two peaks have the same height, and are
  ## within the same 'window', we need to drop one.
  
  ## If a peak has a 'match' at half the size, use the smaller peak (ie.,
  ## take the G1 peak in cases where the G2 peak is higher) 

  ## After the first peak is selected, only consider peaks that are not a
  ## multiple of the size of this peak when selecting the next one.

  peaks <- peaks[order(peaks[,2], decreasing = TRUE), ]

  ## eliminate the debris field?
  peaks <- peaks[which(peaks[, "mean"] > 40), ]

  drop <- numeric()
  for(i in 2: nrow(peaks)){
    if((peaks[i-1, "height"] == peaks[i, "height"]) &
       (abs(peaks[i-1, "mean"] - peaks[i, "mean"]) <= window)){ 
      ## It's a tie!
      drop <- c(drop, i)
    }
  }

  if(length(drop) > 0){                  # there was at least one tie 
    peaks <- peaks[-drop, ]
  }
  
  out <- matrix(NA, nrow = 0, ncol = 2)

  while(nrow(peaks) > 0){
    ## which peaks are half or double the size of the first peak:
    paircheck <-
      which(((peaks[, "mean"] < 0.53 * peaks[1, "mean"]) &
             (peaks[, "mean"] > 0.47 * peaks[1, "mean"])) |
            ((peaks[, "mean"] < 2.13 * peaks[1, "mean"]) &
             (peaks[, "mean"] > 1.89 * peaks[1, "mean"])))
    ## Add the first peak to that list:
    paircheck <- c(1, paircheck)
    if(length(paircheck) == 1){            # no pairs
      out <- rbind(out, peaks[1, ])
      peaks <- peaks[-1, , drop = FALSE]              # remove peak
    } else if (length(paircheck == 2)) {              # pick the smallest of the pair
      out <- rbind(out, peaks[paircheck[which.min(peaks[paircheck, "mean"])], ])
      peaks <- peaks[-paircheck, , drop = FALSE]      # remove pair
    } else {
      warning("paircheck found more than 2 peaks")
      ## This is a bit convoluted. Hopefully unnecessary now that smoothing
      ## is applied before peak finding?

      ## while (length(paircheck > 2)) {
      ##   if(peaks[paircheck[1], "height"] < peaks[paircheck[2], "height"] &
      ##      peaks[paircheck[1], "height"] < peaks[paircheck[3], "height"]){
      ##     peaks <- peaks[-1, , drop = FALSE]
      ##     paircheck <- paircheck[-1]
      ##   } else if(peaks[paircheck[length(paircheck)], "height"] <
      ##             peaks[paircheck[length(paircheck) - 1], "height"] &
      ##             peaks[paircheck[length(paircheck)], "height"] <
      ##             peaks[paircheck[length(paircheck) - 2], "height"]) {
      ##     peaks <- peaks[-3, , drop = FALSE]
      ##     paircheck <- paircheck[-length(paircheck)]
      ##   } else {
      ##     stop("Funky peaks, I cannot find starting values")
      ##   }
    }

  }

  if(is.vector(peaks))
    out <- rbind(out, peaks)

  rownames(out) <- NULL

  out <- out[1:2, ]
  out <- out[order(out[, "mean"]), ]
  out
}
  
