## Functions for building non-linear models for application to flowHist
## objects. 

fA1 <- function(a1, Ma, Sa, xx){
  (a1 / (sqrt(2 * pi) * Sa) * exp(-((xx - Ma)^2)/(2 * Sa^2)))
}
attr(fA1, "compName") <- "fA1"

fA2 <- function(a2, Sa, Ma, xx){
  (a2 / (sqrt(2 * pi) * Sa * 2) * exp(-((xx - Ma * 2)^2)/(2 * (Sa * 2)^2)))
}
attr(fA2, "compName") <- "fA2"

fB1 <- function(b1, Mb, Sb, xx){
  (b1 / (sqrt(2 * pi) * Sb) * exp(-((xx - Mb)^2)/(2 * Sb^2)))
}
attr(fB1, "compName") <- "fB1"

fB2 <- function(b2, Sb, Mb, xx){
  (b2 / (sqrt(2 * pi) * Sb * 2) * exp(-((xx - Mb * 2)^2)/(2 * (Sb * 2)^2)))
}
attr(fB2, "compName") <- "fB2"

derf <- function(x1, x2, s1, s2, x){
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
  (erf((x2 - x) / sqrt(s2)) - erf((x1 -x) / sqrt(s1)))/ 2
}  

broadenedTrapezoid <- function(x1, x2, t1, t2, s1, s2, xx){
  (((t2 - t1)/(x1)) * (xx - x2) + t2) * derf(x1, x2, s1, s2, xx)
}

btA <- function(Ma, t1a, t2a, s1a, s2a, xx){
  broadenedTrapezoid(x1 = Ma, x2 = min(Ma * 2, 255), t1 = t1a, t2 = t2a,
                     s1 = s1a, s2 = s2a, xx)  
}
attr(btA, "compName") <- "btA"

btB <- function(Mb, t1b, t2b, s1b, s2b, xx){
  broadenedTrapezoid(x1 = Mb, x2 = min(Mb * 2, 255), t1 = t1b, t2 = t2b,
                     s1 = s1b, s2 = s2b, xx) 
}
attr(btB, "compName") <- "btB"

broadRect <- function(x1, x2, r, s1, s2, xx){
  r * derf(x1 = x1, x2 = x2, s1 = s1, s2 = s2, x = xx)
}

brA <- function(Ma, ra, s1a, xx){
  broadRect(x1 = Ma, x2 = Ma * 2, r = ra, s1 = s1a, s2 = s1a, xx)
}
attr(brA, "compName") <- "brA"

brB <- function(Mb, rb, s1b, xx){
  broadRect(x1 = Mb, x2 = Mb * 2, r = rb, s1 = s1b, s2 = s1b, xx)
}
attr(brB, "compName") <- "brB"

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
  ## Do not extend the model below/beyond the data
  ## Modfit appears to cut off the debris slightly above the lowest data,
  ## which gives a better fit. Perhaps set first.channel to 2-4? Need to
  ## test this and determine best fit. Possibly use an extra parameter to
  ## tune this for each data set individually.
  first.channel <- which(intensity > 0)[2]

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
singleCutVect <- Vectorize(singleCutBase, "xx")

singleCut <- function(SCa, intensity, xx){
  singleCutVect(SCa, intensity, xx)
}
attr(singleCut, "compName") <- "single cut"

flowInit <- function(fh) {
  xy <- fh$data

  peaks <- fh$peaks

  params <- names(formals(fh$model))
  params <- params[-which(params %in% c("", "xx", "intensity"))]
  value <- c()
  
  if("Ma" %in% params) {
    ## Any model with Ma will require all three of these parameters:
    Ma <- peaks[1, "mean"]  
    Sa <- Ma / 20                         # assume CV = 0.05
    a1 <- peaks[1, "height"] * Sa / 0.4
    ra <- peaks[1, "height"] / 100
    ## t1a <- peaks[1, "height"] / 100
    ## t2a <- t1a
    s1a <- 1
    tmpval <- c(Ma, Sa, a1, ra, s1a)
    names(tmpval) <- c("Ma", "Sa", "a1", "ra", "s1a")
    value <- c(value, tmpval)
  }

  if("a2" %in% params) {
    ## Is a2 off the chart? It shouldn't be! Models with an a2 peak can
    ## break if the a2 peak is beyond the data range.
    if((peaks[1, "mean"] * 2) > max(xy[ ,"x"])){
      warning("a2 peak appears to be out of range")
      a2 <- 0
    } else {
      a2 <- xy[peaks[1, "mean"] * 2, "intensity"] * Sa * 2 / 0.4
    }
    tmpval <- c(a2)
    names(tmpval) <- c("a2")
    value <- c(value, tmpval)
  }

  if("Mb" %in% params){
    Mb <- peaks[2, "mean"]
    Sb <- Mb / 20
    b1 <- peaks[2, "height"] * Sb / 0.4
    rb <- peaks[1, "height"] / 100
    ##t1b <- peaks[2, "height"] / 100
    ##t2b <- t1b
    s1b <- 1
    tmpval <- c(Mb, Sb, b1, rb, s1b)
    names(tmpval) <- c("Mb", "Sb", "b1", "rb", "s1b")
    value <- c(value, tmpval)
  }

  if("b2" %in% params) {
    if((peaks[2, "mean"] * 2) > max(xy[,"x"])){
      warning("b2 peak appears to be out of range")
      b2 <- 0
    } else {
      b2 <- as.vector(xy[peaks[2, "mean"] * 2, "intensity"] * Sb * 2 / 0.4)
    }
    tmpval <- c(b2)
    names(tmpval) <- c("b2")
    value <- c(value, tmpval)
  }

  if("SCa" %in% params){
    SCa <- 0.1                            # just a wild guess for now
    tmpval <- c(SCa)
    names(tmpval) <- c("SCa")
    value <- c(value, tmpval)
  }

  as.list(value)
}

makeModel <- function(components, env = parent.frame()){

  args <- unlist(lapply(components, formals))
  args <- args[unique(names(args))]
  
  bodList <- lapply(components, FUN = body)
  bod <- bodList[[1]]
  bodList <- bodList[-1]

  while(length(bodList) > 0){
    bod <- call("+", bod, bodList[[1]])
    bodList <- bodList[-1]
  }

  eval(call("function", as.pairlist(args), bod), env)

}

