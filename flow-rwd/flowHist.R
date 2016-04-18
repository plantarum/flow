## Functions for creating and viewing flowHist objects.
flowHist <- function(FCS = NULL, FILE = NULL, CHANNEL,
                     bins = 256, maxBins = 1024){ 
  ## You probably want to subdivide the bins evenly. i.e., if there are
  ## 1024 bins in the data, use 128, 256, 512 bins
  if((1024 %% bins) != 0)
    warning("maxBins is not a multiple of bins!")

  if(sum(c(is.null(FCS), is.null(FILE))) != 1){
    stop("\nOne (and only one) of FCS or FILE must be set.")
  }

  if(!is.null(FILE))
    FCS <- read.FCS(FILE, dataset = 1, alter.names = TRUE)
  
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
              nls = NULL, standard = NULL,
              file = FCS@description$GUID)

  peaks <- cleanPeaks(findPeaks(res, 20, 20), 20)  

  attr(res$data, "peaks") <- peaks
  res$comps <- list(singleCut, fA1, fB1)

  if(peaks[1, "mean"] * 2 <= nrow(res$data))
    res$comps <- c(res$comps, fA2)

  if(peaks[2, "mean"] * 2 <= nrow(res$data))
    res$comps <- c(res$comps, fB2)

  res$model = makeModel(res$comps, env = globalenv())
  
  class(res) <- "flowHist"
  
  return(res)
}

print.flowHist <- function(self){
  message("flowHist object")
  message("Source file: ", self$file)
  message("Channel: ", self$channel)
  message("Values: ", dim(self$data)[1])
  message("Total events: ", sum(self$data$intensity))
  message("Model components: ",
          paste(unlist(lapply(fh1$comps,
                              FUN = function(x) attr(x, "compName"))),
                collapse = ", ") 
          )
  
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
    message(paste("Peak A:", round(self$counts$firstPeak$value, 1), " at ",
                  round(coef(self$nls)["Ma"], 1)))
    message(paste("Peak B:", round(self$counts$secondPeak$value, 1), " at ",
                  round(coef(self$nls)["Mb"], 1)))
  }

  if(!is.null(self$cv)){
    message(paste("CV A:", round(self$cv$CVa, 3)))
    message(paste("CV B:", round(self$cv$CVb, 3)))
    message(paste("Ratio:", round(self$cv$CI[1], 3), "/",
                  round(1/self$cv$CI[1], 3))) 
  }

}

plot.flowHist <- function(self, init = FALSE, nls = TRUE, comps = TRUE){
  plot(self$data$intensity, type = 'n', main = self$file)
  polygon(x = c(self$data$x, max(self$data$x) + 1), y = c(self$data$intensity, 0),
          col = "lightgray", border = NA)

  if(init){
    iv <- fHcall(self, "getInitial")
    yy <- do.call(self$model,
                  args = c(list(intensity = self$data$intensity,
                                xx = self$data$x),
                           as.list(fHcall(self, "getInitial"))))

    lines(x = self$data$x,
          y = yy, 
          col = 1, lwd = 3, lty = 5)
  }
  
  if(nls & (! is.null(self$nls))){
    lines(x = self$data$x, y = predict(self$nls), col = 2)
  }
  
  if(comps & (! is.null(self$nls))){
    for(i in seq_along(self$comps)){
      if("intensity" %in% names(formals(self$comps[[i]]))){
        params <-
          as.list(coef(self$nls)[names(formals(self$comps[[i]]))])
        params <- params[! is.na(names(params))]
        yy <- do.call(self$comps[[i]],
                      args = c(list(intensity = self$data$intensity,
                                    xx = self$data$x),
                               params))
        lines(x = self$data$x, y = yy, col = i + 2)
      } else {
        params <-
          as.list(coef(self$nls)[names(formals(self$comps[[i]]))])
        params <- params[! is.na(names(params))]
        yy <- do.call(self$comps[[i]],
                      args = c(list(xx = self$data$x), params))
        lines(x = self$data$x, y = yy, col = i + 2)
      }
    }
  }
}

