library(flowCore)
library(car)                            # for deltaMethod
source("flowHist.R")
source("flowComponents.R")
source("model-utils.R")
source("flowAnalysis.R")

fh1 <- flowHist(FILE =
                  "test-data/LMD good/188-15_+_rad_2015-04-20_647.LMD",
                CHANNEL = "FL3.INT.LIN", window = 20, smooth = 10)
plot(fh1, init = TRUE)

fh1 <- fhAnalyze(fh1)
plot(fh1, comps = TRUE)
fh1

fh2 <- flowHist(FILE = "test-data/LMD good/SM239 + R 2014-12-01 852.LMD",
                CHANNEL = "FL3.INT.LIN")                
plot(fh2, init = TRUE)

fh2 <- fhAnalyze(fh2)
plot(fh2, comps = TRUE)
fh2

setwd("~/research/techniques/flow/flow-rwd/test-data/camelina_test")
lmdFiles <- list.files()

out <- list()

i <- lmdFiles[1]
i <- lmdFiles[2]
i <- lmdFiles[3]

i <- lmdFiles[6]

## Problem!! intial peaks are not accurate enough after all the smoothing.
## Need a way to manually select peak positions for this.
i <- lmdFiles[4]

## Problem!! only two peaks, one twice the other -- good data?
i <- lmdFiles[5]

i <- lmdFiles[7]
i <- lmdFiles[8]
i <- lmdFiles[9]
i <- lmdFiles[10]
i <- lmdFiles[11]
i <- lmdFiles[12]
i <- lmdFiles[13]


fh1 <- flowHist(FILE = i, CHANNEL = "FL3.INT.LIN")
plot(fh1, init = TRUE)

fh1 <- fhAnalyze(fh1)

fh1 <- pickInit(fh1)

fh1 <- fhAnalyze(fh1)
plot(fh1, comps = TRUE)
fh1


