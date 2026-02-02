###########################################
### downstream delays and probabilities ###
###########################################

## The time-dependent delay distributions were calculated up to a maximum delay of 100 days
## For efficiency this maximum can be decreased, but the distributions MUST add up to 1. 
## This normalisation is done here
normalisedelays <- function(delays, maxdelay) {
  toreturn <- delays[1:(1 + maxdelay), ]
  toreturn <- toreturn / rep(colSums(toreturn), each = maxdelay + 1)
  return(toreturn)
}

## The time-dependent delay distributions are defined forward in time, i.e. given infection/admission
## at time t, what is the probability to become symptomatic/be discharged at time t+1, t+2, etc.
## To calculate symptomatic/discharge incidence at time t, each distribution for t-1, t-2 etc. must
## be used once, only at the single time point of the matching delay. For efficiency, it is better
## to tilt the distributions.
tiltdelays <- function(delays, tmax, maxdelay) {
  # add square matrix of columns with final delay distribution
  toreturn <- cbind(delays, matrix(delays[, tmax], nrow = maxdelay + 1, ncol = maxdelay + 1))
  # transpose and remove final maxdelay+1 elements (from the final column, so matrix structure disappears)
  toreturn <- t(toreturn)[1:((tmax + maxdelay) * (maxdelay + 1))]
  # create new matrix with same number of columns, so that number of rows is decreased by 1, and each column
  #   is shifted downwards by one more place than the previous column, then transpose backwards
  toreturn <- t(matrix(toreturn, ncol = maxdelay + 1))
  # make left lower triangle equal to 0
  toreturn[cbind(rep(2:(maxdelay + 1), 1:maxdelay),
                 unlist(sapply(1:maxdelay, function(x) 1:x)))] <- 0
  toreturn <- toreturn[, 1:tmax]
  return(toreturn)
}


############################
### outside the hospital ###
############################
# most defined in OSIRISanalyses4model

## I2Se - probability (hospitalisable)
# defined in ContactsInfectivities
probI2Se <- function(probi2se = "default", nr = 1) {
  if(probi2se[1] == "sample") {
    return(PreAdmissionProbs$probI2Se$samples[nr, ])
  } else if(probi2se[1] == "default") {
    return(PreAdmissionProbs$probI2Se$default)
  } else {
    return(probi2se)
  }
}


## I2S - delay (incubation period)
delayI2S <- function(maxdelay = 100) {
  if(maxdelay < 100) {
    toreturn <- PreAdmissionDelays$delayI2S[1:(1 + maxdelay)]
    toreturn <- toreturn / sum(toreturn)
  } else {
    toreturn <- PreAdmissionDelays$delayI2S
  }
  return(toreturn)
}

## I2R - delay
delayI2R <- function() {
  return(PreAdmissionDelays$delayI2R)
}


## S2Se - delays 
delayS2Se <- function(maxtime = as.numeric(ANALYSISDATE - as.Date("2020-02-12")),
                      maxdelay = 100) {
  if(maxtime > dim(PreAdmissionDelays$delayS2A)[2]) {
    curdims <- dim(PreAdmissionDelays$delayS2A)
    maxtimeinarray <- curdims[2]
    toreturn <- aperm(PreAdmissionDelays$delayS2A, c(1, 3, 2))
    toreturn <- array(c(toreturn, 
                        rep(toreturn[, , maxtimeinarray], maxtime - maxtimeinarray)),
                      dim = c(curdims[1], curdims[3], maxtime))
    toreturn <- aperm(toreturn, c(1, 3, 2))
  } else if(maxtime < dim(PreAdmissionDelays$delayS2A)[2]) {
    toreturn <- PreAdmissionDelays$delayS2A[, 1:maxtime, ]
  } else {
    toreturn <- PreAdmissionDelays$delayS2A
  }
  
  if(maxdelay < 100) {
    toreturn <- apply(toreturn, 3, normalisedelays, maxdelay = maxdelay)
    dim(toreturn) <- c(maxdelay + 1, maxtime, 9)
  }
  
  toreturn <- apply(toreturn, 3, tiltdelays, tmax = maxtime, maxdelay = maxdelay)
  dim(toreturn) <- c(maxdelay + 1, maxtime, 9)
  
  return(toreturn)
}

## I2Se - delays
delayI2Se <- function(maxtime = as.numeric(ANALYSISDATE - as.Date("2020-02-12")),
                      maxdelay = 100) {
  if(maxtime > dim(PreAdmissionDelays$delayI2A)[2]) {
    curdims <- dim(PreAdmissionDelays$delayI2A)
    maxtimeinarray <- curdims[2]
    toreturn <- aperm(PreAdmissionDelays$delayI2A, c(1, 3, 2))
    toreturn <- array(c(toreturn, 
                        rep(toreturn[, , maxtimeinarray], maxtime - maxtimeinarray)),
                      dim = c(curdims[1], curdims[3], maxtime))
    toreturn <- aperm(toreturn, c(1, 3, 2))
  } else if(maxtime < dim(PreAdmissionDelays$delayI2A)[2]) {
    toreturn <- PreAdmissionDelays$delayI2A[, 1:maxtime, ]
  } else {
    toreturn <- PreAdmissionDelays$delayI2A
  }
  
  if(maxdelay < 100) {
    toreturn <- apply(toreturn, 3, normalisedelays, maxdelay = maxdelay)
    dim(toreturn) <- c(maxdelay + 1, maxtime, 9)
  }
  
  toreturn <- apply(toreturn, 3, tiltdelays, tmax = maxtime, maxdelay = maxdelay)
  dim(toreturn) <- c(maxdelay + 1, maxtime, 9)
  
  return(toreturn)
}


### FROM HERE: NICE results

# Se2A - probabilities (which severe are admitted)
# age-dependent (col), time-dependent (row)
probSe2A <- function(maxtime = as.numeric(ANALYSISDATE - as.Date("2020-02-12"))) {
  if(maxtime > nrow(NICEprobabilities$probSe2A)) {
    maxtimeinmatrix<- nrow(NICEprobabilities$probSe2A)
    toreturn <- NICEprobabilities$probSe2A
    toreturn <- rbind(toreturn, 
                      matrix(rep(toreturn[maxtimeinmatrix, ], each = maxtime - maxtimeinmatrix),
                             ncol = ncol(toreturn)))
  } else if(maxtime < nrow(NICEprobabilities$probSe2A)) {
    toreturn <- NICEprobabilities$probSe2A[1:maxtime, ]
  } else {
    toreturn <- NICEprobabilities$probSe2A
  }
  
  return(toreturn)
}

# A2IC - probabilities
# age-dependent (col), time-dependent (row)
probSe2A2IC <- function(maxtime = as.numeric(ANALYSISDATE - as.Date("2020-02-12"))) {
  if(maxtime > nrow(NICEprobabilities$probSe2A2IC)) {
    maxtimeinmatrix<- nrow(NICEprobabilities$probSe2A2IC)
    toreturn <- NICEprobabilities$probSe2A2IC
    toreturn <- rbind(toreturn, 
                      matrix(rep(toreturn[maxtimeinmatrix, ], each = maxtime - maxtimeinmatrix),
                             ncol = ncol(toreturn)))
  } else if(maxtime < nrow(NICEprobabilities$probSe2A2IC)) {
    toreturn <- NICEprobabilities$probSe2A2IC[1:maxtime, ]
  } else {
    toreturn <- NICEprobabilities$probSe2A2IC
  }
  
  return(toreturn)
}

# A2IC - delays
# from NICE
delayA2IC <- function(maxtime = as.numeric(ANALYSISDATE - as.Date("2020-02-12")),
                      maxdelay = 100) {
  if(maxtime > ncol(NICEdelays$delayA2IC)) {
    maxtimeinmatrix<- ncol(NICEdelays$delayA2IC)
    toreturn <- NICEdelays$delayA2IC
    toreturn <- cbind(toreturn, 
                      matrix(rep(toreturn[, maxtimeinmatrix], maxtime - maxtimeinmatrix),
                             nrow = nrow(toreturn)))
  } else if(maxtime < ncol(NICEdelays$delayA2IC)) {
    toreturn <- NICEdelays$delayA2IC[, 1:maxtime]
  } else {
    toreturn <- NICEdelays$delayA2IC
  }
  
  if(maxdelay < 100) {
    toreturn <- normalisedelays(toreturn, maxdelay)
  }
  
  toreturn <- tiltdelays(toreturn, maxtime, maxdelay)
  
  return(toreturn)
}

# A2D - probabilities 
# age-dependent (col), time-dependent (row)
probSe2A2D <- function(maxtime = as.numeric(ANALYSISDATE - as.Date("2020-02-12"))) {
  if(maxtime > nrow(NICEprobabilities$probSe2A2DC)) {
    maxtimeinmatrix<- nrow(NICEprobabilities$probSe2A2DC)
    toreturn <- NICEprobabilities$probSe2A2DC
    toreturn <- rbind(toreturn, 
                      matrix(rep(toreturn[maxtimeinmatrix, ], each = maxtime - maxtimeinmatrix),
                             ncol = ncol(toreturn)))
  } else if(maxtime < nrow(NICEprobabilities$probSe2A2DC)) {
    toreturn <- NICEprobabilities$probSe2A2DC[1:maxtime, ]
  } else {
    toreturn <- NICEprobabilities$probSe2A2DC
  }
  
  return(toreturn)
}

# A2D - delays
# from NICE
delayA2D <- function(maxtime = as.numeric(ANALYSISDATE - as.Date("2020-02-12")),
                     maxdelay = 100) {
  if(maxtime > ncol(NICEdelays$delayA2D)) {
    maxtimeinmatrix<- ncol(NICEdelays$delayA2D)
    toreturn <- NICEdelays$delayA2D
    toreturn <- cbind(toreturn, 
                      matrix(rep(toreturn[, maxtimeinmatrix], maxtime - maxtimeinmatrix),
                             nrow = nrow(toreturn)))
  } else if(maxtime < ncol(NICEdelays$delayA2D)) {
    toreturn <- NICEdelays$delayA2D[, 1:maxtime]
  } else {
    toreturn <- NICEdelays$delayA2D
  }
  
  if(maxdelay < 100) {
    toreturn <- normalisedelays(toreturn, maxdelay)
  }
  
  toreturn <- tiltdelays(toreturn, maxtime, maxdelay)
  
  return(toreturn)
}
  

# IC2D - probabilities
# age-dependent (col), time-dependent (row)
probSe2A2IC2D <- function(maxtime = as.numeric(ANALYSISDATE - as.Date("2020-02-12"))) {
  if(maxtime > nrow(NICEprobabilities$probSe2A2IC2DC)) {
    maxtimeinmatrix<- nrow(NICEprobabilities$probSe2A2IC2DC)
    toreturn <- NICEprobabilities$probSe2A2IC2DC
    toreturn <- rbind(toreturn, 
                      matrix(rep(toreturn[maxtimeinmatrix, ], each = maxtime - maxtimeinmatrix),
                             ncol = ncol(toreturn)))
  } else if(maxtime < nrow(NICEprobabilities$probSe2A2IC2DC)) {
    toreturn <- NICEprobabilities$probSe2A2IC2DC[1:maxtime, ]
  } else {
    toreturn <- NICEprobabilities$probSe2A2IC2DC
  }
  
  return(toreturn)
}
  

# IC2D - delays (or IC2H)
# from NICE
delayA2IC2D <- function(maxtime = as.numeric(ANALYSISDATE - as.Date("2020-02-12")),
                        maxdelay = 100) {
  if(maxtime > ncol(NICEdelays$delayA2IC2D)) {
    maxtimeinmatrix<- ncol(NICEdelays$delayA2IC2D)
    toreturn <- NICEdelays$delayA2IC2D
    toreturn <- cbind(toreturn, 
                      matrix(rep(toreturn[, maxtimeinmatrix], maxtime - maxtimeinmatrix),
                             nrow = nrow(toreturn)))
  } else if(maxtime < ncol(NICEdelays$delayA2IC2D)) {
    toreturn <- NICEdelays$delayA2IC2D[, 1:maxtime]
  } else {
    toreturn <- NICEdelays$delayA2IC2D
  }
  
  if(maxdelay < 100) {
    toreturn <- normalisedelays(toreturn, maxdelay)
  }
  
  toreturn <- tiltdelays(toreturn, maxtime, maxdelay)
  
  return(toreturn)
}
  
  
# IC2H - probabilities
# from NICE
probSe2A2IC2H <- function(maxtime = as.numeric(ANALYSISDATE - as.Date("2020-02-12"))) {
  if(maxtime > nrow(NICEprobabilities$probSe2A2IC2H)) {
    maxtimeinmatrix<- nrow(NICEprobabilities$probSe2A2IC2H)
    toreturn <- NICEprobabilities$probSe2A2IC2H
    toreturn <- rbind(toreturn, 
                      matrix(rep(toreturn[maxtimeinmatrix, ], each = maxtime - maxtimeinmatrix),
                             ncol = ncol(toreturn)))
  } else if(maxtime < nrow(NICEprobabilities$probSe2A2IC2H)) {
    toreturn <- NICEprobabilities$probSe2A2IC2H[1:maxtime, ]
  } else {
    toreturn <- NICEprobabilities$probSe2A2IC2H
  }
  
  return(toreturn)
}  
  


# H22D - delays
# from NICE
delayA2IC2H2D <- function(maxtime = as.numeric(ANALYSISDATE - as.Date("2020-02-12")),
                          maxdelay = 100) {
  if(maxtime > ncol(NICEdelays$delayA2IC2H2D)) {
    maxtimeinmatrix<- ncol(NICEdelays$delayA2IC2H2D)
    toreturn <- NICEdelays$delayA2IC2H2D
    toreturn <- cbind(toreturn, 
                      matrix(rep(toreturn[, maxtimeinmatrix], maxtime - maxtimeinmatrix),
                             nrow = nrow(toreturn)))
  } else if(maxtime < ncol(NICEdelays$delayA2IC2H2D)) {
    toreturn <- NICEdelays$delayA2IC2H2D[, 1:maxtime]
  } else {
    toreturn <- NICEdelays$delayA2IC2H2D
  }
  
  if(maxdelay < 100) {
    toreturn <- normalisedelays(toreturn, maxdelay)
  }
  
  toreturn <- tiltdelays(toreturn, maxtime, maxdelay)
  
  return(toreturn)
}
  


  
customise_delays <- function(tmax, maxdelay, likelihood_only = FALSE) {
  toreturn <- list(
    dI2Se = delayI2Se(tmax, maxdelay),
    dA2IC = delayA2IC(tmax, maxdelay)
  )
  
  if(!likelihood_only) {
    toreturn <- c(toreturn,
                  list(
                    dI2S = delayI2S(maxdelay),
                    dA2D = delayA2D(tmax, maxdelay),
                    dA2IC2D = delayA2IC2D(tmax, maxdelay),
                    dA2IC2H2D = delayA2IC2H2D(tmax, maxdelay)
                  ))
  }
  
  return(toreturn)
}


