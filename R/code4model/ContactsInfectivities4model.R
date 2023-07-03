#######################################################
### contacts, infectivities,                        ###
###   susceptibilities, symptomatic fractions,      ###
###   natural history and control                   ###
#######################################################





# susceptibility/infectivity by age group
LogRelativeSusInf <- function(logrelsusinf = "default", nr = 1) {
  if(logrelsusinf[1] == "sample") {
    return(ContactModelInput$LogRelSusInf$samples[nr, ])
  } else if(logrelsusinf[1] == "default") {
    return(ContactModelInput$LogRelSusInf$default)
  } else {
    return(logrelsusinf)
  }
}



# infectivity profile
InfectivityProfile <- function(infectivity = 9 * 9 * 2.3, ...) {
  return(infectivity * ContactModelInput$InfCurves)
}


# initial number of infections
LogY0 <- function(logy0, nr = 1) {
  if(logy0 == "sample") {
    return(ContactModelInput$LogY0$samples[nr])
  } else {
    return(logy0)
  }
}


# infectivities for different time periods
LogInfectivities <- function(loginfectivities, nr = 1) {
  if(loginfectivities[1] == "sample") {
    return(ContactModelInput$LogInfectivities$samples[nr, ])
  } else {
    return(loginfectivities)
  }
}

# relative infectivity to reflect uncertainty in prognoses
InfectivityDeviation <- function(infectivitydeviation = 1, nr = 1) {
  if(infectivitydeviation == "sample") {
    return(ContactModelInput$InfectivityDeviations$samples[nr])
  } else {
    return(infectivitydeviation)
  }
}


# repeat population functions
popsize <- function(ROAZ = "all", ...) {
  if(ROAZ == "all") {
    return(Populationresults$popsizeNL)
  } else {
    return(Populationresults$popsizeregio[ROAZ])
  }
}

agedist <- function(ROAZ = "all", ngroups = 9, ...) {
  if(ROAZ == "all") {
    if(ngroups == 9) {
      return(Populationresults$agedistNL)
    } else {
      return(Populationresults$agedistNL10)
    }
  } else {
    return(Populationresults$agedistregio[, ROAZ])
  }
}


