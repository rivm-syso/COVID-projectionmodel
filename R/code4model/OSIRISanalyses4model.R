##################################################################
### Functions to estimate probabilities and delays from OSIRIS ###
##################################################################

#####################
### Read the data ###
#####################

file.name <- "data/OSIRIS/OSIRISdata_20210106_synthetic.rds"
dataOsiris <- readRDS(file.name)
osirisdata4fit <- dataOsiris %>%
  select(ZIE1eZiekteDt, NCOVdat1ezkhopn, Leeftijdsgroep) %>%
  filter(!is.na(ZIE1eZiekteDt) & !is.na(NCOVdat1ezkhopn)) %>%
  mutate(interval = as.numeric(NCOVdat1ezkhopn - ZIE1eZiekteDt)) %>%
  filter(interval >= 0 & interval < 31) %>% 
  separate(Leeftijdsgroep, c("onder", "boven")) %>% 
  filter(onder != "Niet") %>%
  mutate(onder = as.numeric(onder),
         age_class = floor(onder / 10),
         age_class = if_else(age_class > 8, 8, age_class),
         ti = as.numeric(ZIE1eZiekteDt - as.Date("2020-02-12"))) %>%
  select(ti, age_class, interval)

#############################
### Infection to symptoms ###
#############################
# mean 5, sd 2.5, Backer but a bit lower
get_I2S_delays <- function(maxdelay = 100) {
  cumdist <- pweibull(seq(.5, 0.5 + maxdelay, 1), 2.10, 5.65)
  return(cumdist - c(0, head(cumdist, -1)))
}

#######################################
### Infection to immune (mean only) ###
#######################################
# assume immunity after 19 days
get_I2R_delays <- function() {
  return(19)
}


###############################################
### Symptom to hospital, age group-specific ###
###############################################
meanDT <- function(t, mu1, mu2, rate1, tmid1) {
  mu1 +
    (mu2 - mu1) * exp(rate1 * (t - tmid1)) / (1 + exp(rate1 * (t - tmid1)))
}
minloglikS2A <- function(parms, ints, days, ages) {
  - sum(dnbinom(ints, exp(parms[21]),
                mu = meanDT(days, mu1 = exp(parms[1 + ages]),
                            mu2 = exp(parms[10 + ages]),
                            rate1 = exp(parms[19]),
                            tmid1 = exp(parms[20])), log = T
  ))
}
estresult <- optim(c(log(rep(5,18)),0,5,1), minloglikS2A, 
      ints = osirisdata4fit$interval, days = osirisdata4fit$ti, ages = osirisdata4fit$age_class)$par
bestfitparmsosiris <- exp(estresult)

get_S2A_delays <- function(tmax, maxdelay = 100) {
  allt <- matrix(rep(1:tmax, 9), ncol = 9)
  allmu1 <- matrix(rep(bestfitparmsosiris[1:9], each = tmax), ncol = 9)
  allmu2 <- matrix(rep(bestfitparmsosiris[10:18], each = tmax), ncol = 9)
  
  meandts <- meanDT(allt, allmu1, allmu2, 
                    bestfitparmsosiris[19],
                    bestfitparmsosiris[20])
  
  result <- array(
    sapply(meandts, function(x) dnbinom(0:maxdelay, bestfitparmsosiris[21], mu = x)),
    dim = c(maxdelay + 1, tmax, 9))
  
  return(
    result
  )
}


#############################
### Infection to hospital ###
#############################
get_I2A_delays <- function(tmax, maxdelay = 100) {
  delayi2s <- get_I2S_delays(maxdelay)
  delays2se <- get_S2A_delays(tmax, maxdelay)
  delayi2se <- array(
    # every age class
    apply(delays2se, 3, 
          function(x) 
            # every point in time
            sapply(1:tmax, 
                   function(j) 
                     # probability of delay (0 to naxdelay)
                     sapply(1:(1 + maxdelay), 
                            function(i) 
                              # over all incubation times <= delay
                              sum(delayi2s[1:i] * x[cbind(i:1, pmin(j:(j + i - 1), tmax))])))), 
    dim = c(1 + maxdelay, tmax, 9))
  return(delayi2se)
}




###############################
### Create the results list ###
###############################
PreAdmissionDelays <- list(
  delayI2S = get_I2S_delays(),
  delayI2R = get_I2R_delays(),
  delayS2A = get_S2A_delays(as.numeric(ANALYSISDATE - as.Date("2020-02-12"))),
  delayI2A = get_I2A_delays(as.numeric(ANALYSISDATE - as.Date("2020-02-12")))
  )

rm(file.name, osirisdata4fit,
   get_I2S_delays, get_I2R_delays, get_S2A_delays, get_I2A_delays,
   meanDT, minloglikS2A, estresult, 
   bestfitparmsosiris
   )





