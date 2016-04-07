library(flowCore)
library(car)                            # for deltaMethod

lmd1 <- read.FCS("test-data/LMD good/188-15_+_rad_2015-04-20_647.LMD",
                 dataset = 1, alter.names = TRUE)
fh1 <- flowHist(lmd1, "FL3.INT.LIN")
fh1
fh1$nls <- fHnls(fh1)
fh1$model <- SSoneSSC
plot(fh1)
fhIntegrate(fh1)
fhCV(fh1)

lmd2 <- read.FCS("test-data/LMD good/SM239 + R 2014-12-01 852.LMD",
                 dataset = 1, alter.names = TRUE)
fh2 <- flowHist(lmd2, "FL3.INT.LIN")
fh2
fh2$nls <- fHnls(fh2)
fh2$model <- SSoneSSC
plot(fh2)
fhIntegrate(fh2)
fhCV(fh2)

setwd("~/research/techniques/flow/flow-rwd/test-data/camelina_test")
lmdFiles <- list.files()

out <- list()
for(i in lmdFiles){
  raw <- read.FCS(i, dataset = 1, alter.names = TRUE)
  fhi <- flowHist(raw, "FL3.INT.LIN")
  print(fhi)
  plot(fhi, init = TRUE)
  fhi$nls <- fHnls(fhi)
  fhi$model <- SSoneSSC
  plot(fhi)
  fhi$counts <- fhIntegrate(fhi)
  fhi$cv <- fhCV(fhi)
  out[[i]] <- fhi
}

fl3.trimmed <- exprs(lmd1)[, "FL3.INT.LIN"][exprs(lmd1)[, "FL3.INT.LIN"] <
                                              max(exprs(lmd1)[, "FL3.INT.LIN"])]
## Set one more break than the number of channels we want!
fl3.1024 <- hist(fl3.trimmed, breaks = seq(from = 0, to=1024, by = 1))
fl3.256 <- hist(fl3.trimmed, breaks = seq(from = 0, to = 1024, by = 4),
                plot = FALSE)
fl3.df <- data.frame(x = 1:length(fl3.256$counts),
                     intensity = fl3.256$counts) 

plot(fl3.df$intensity, type = 'n')
polygon(x = c(fl3.df$x, 257), y = c(fl3.df$intensity, 0),
                col = "lightgray", border = NA) 


m3test <- nls(intensity ~ SSoneSSC(a1, Ma, Sa, a2, b1, Mb, Sb, b2, SCa,
                                   intensity = intensity, xx = x),
              data = fl3.df)

## Complete model
lines(x = fl3.df$x, y = predict(m3test), col = 2)

## Debris:
lines(x = fl3.df$x,
      y = singleCut(SCa = coef(m3test)["SCa"],
                    intensity = fl3.df$intensity, xx = fl3.df$x),
      col = 3, lty = 2, lwd = 2)
##xx <- fl3.df$x                          # awkward!!

## Internal Standard:
lines(x = fl3.df$x,
      y = gaussDipA(Ma = coef(m3test)["Ma"], a1 = coef(m3test)["a1"],
                    Sa = coef(m3test)["Sa"], a2 = coef(m3test)["a2"],
                    xx = fl3.df$x),
      col = 4)

## Sample:
lines(x = fl3.df$x,
      y = gaussDipB(Mb = coef(m3test)["Mb"], b1 = coef(m3test)["b1"],
                    Sb = coef(m3test)["Sb"], b2 = coef(m3test)["b2"],
                    xx = fl3.df$x),
      col = 5)


legend(legend = c("Internal Standard", "Sample", "Debris", "Complete Model"),
       col = c(5, 4, 3, 2), x = "topright", lty = c(1, 1, 2, 1))

testFun <- function(x){
  oneSampleSC(xx = x, a1 = coef(m3test)["a1"], Ma = coef(m3test)["Ma"],
              Sa = coef(m3test)["Sa"], a2 = coef(m3test)["a2"],
              b1 = coef(m3test)["b1"], Mb = coef(m3test)["Mb"],
              Sb = coef(m3test)["Sb"], b2 = coef(m3test)["b2"],
              SCa = coef(m3test)["SCa"], intensity = fl3.df$intensity)
}

gaussAmodel <- function(x){
  gaussDipA(xx = x, a1 = coef(m3test)["a1"], Ma = coef(m3test)["Ma"],
              Sa = coef(m3test)["Sa"], a2 = coef(m3test)["a2"])
}



gaussBmodel <- function(x){
  gaussDipB(xx = x, b1 = coef(m3test)["b1"], Mb = coef(m3test)["Mb"],
              Sb = coef(m3test)["Sb"], b2 = coef(m3test)["b2"])
}

modelledEvents <-
  integrate(testFun, 0, 256, subdivisions = 2000)
AEvents <- integrate(gaussAmodel, 0, 256, subdivisions = 2000)
BEvents <- integrate(gaussBmodel, 0, 256, subdivisions = 2000)
CVa <- coef(m3test)["Sa"]/coef(m3test)["Ma"]
CVb <- coef(m3test)["Sb"]/coef(m3test)["Mb"]
CI <- deltaMethod(m3test, "Ma/Mb")

#confint(m3test, "Ma")
#confint(m3test, "Mb")

modelledEvents
AEvents
BEvents
CVa
CVb
c(CI$Estimate - 1.96 * CI$SE, CI$Estimate + 1.96 * CI$SE)

gaussAEvents <-
  integrate(function (x) coef(m3test)["a1"] / (sqrt(2 * pi) * coef(m3test)["Sa"]) *
                         exp(-((x - coef(m3test)["Ma"])^2)/(2 *
                                                        coef(m3test)["Sa"]^2)),
            0, 256)

gaussBEvents <-
  integrate(function (x) coef(m3test)["b1"] / (sqrt(2 * pi) * coef(m3test)["Sb"]) *
                         exp(-((x - coef(m3test)["Mb"])^2)/(2 *
                                                        coef(m3test)["Sb"]^2)),
            0, 256)


