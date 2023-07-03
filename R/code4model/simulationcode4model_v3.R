##############################################
### modelcode                              ###
##############################################


################################################################
### functions to create numeric input for simulation engines ###
################################################################
create_contactarray <- function(controlregimes, nr = 1) {
  toreturn <- ContactMatrices[controlregimes]
  whichsample <- sapply(toreturn, length) > 81
  toreturn[whichsample] <- sapply(toreturn[whichsample], function(x) x[nr])
  toreturn <- do.call(c, toreturn)
  dim(toreturn) <- c(9, 9, length(controlregimes))
  return(toreturn)
}

create_infdeviations <- function(controlregimes, nr = 1) {
  toreturn <- c(InfectivityDeviation("sample", nr), 1)[1 + grepl("_mean", controlregimes)]
  return(toreturn)
}

create_dailyfactors <- function(seasonality, Tmax) {
  SeasonalityCurves %>%
    filter(date <= Tmax + as.Date("2020-02-12")) %>%
    pull(eval(seasonality))
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

customise_probs <- function(tmax, probi2se = "default", nr = 1, likelihood_only = FALSE) {
  toreturn <- list(
    pI2Se = probI2Se(probi2se, nr = nr),
    pSe2A2IC = probSe2A2IC(tmax)
  )
  
  if(!likelihood_only) {
    toreturn <- c(toreturn,
                  list(
                    pSe2A = probSe2A(tmax),
                    pSe2A2D = probSe2A2D(tmax),
                    pSe2A2IC2D = probSe2A2IC2D(tmax),
                    pSe2A2IC2H = probSe2A2IC2H(tmax)
                  ))
  }
  
  return(toreturn)
}

#########################################################
### simulation engines: simulate with numerical input ###
#########################################################
## incidence simulator
engine_incidence <- function(Tmax, # length of simulation
                             y0, # length-9 vector of initial incidence (12-day long from 1 Feb to 12 Feb)
                             tperiods, # vector with final days of periods with different contacts and/or infectivities
                             contacts, # 9*9*tperiods array with contact matrices, one for each period
                             infectivities, # vector with infectivities, one for each period
                             dailyfactors, # Tmax-length vector of multiplication factors to reflect seasonality
                             relsusinf) { # normalised length-9 vector with relative susceptibilities/infectivities
  inc <- matrix(ncol = Tmax + 1, nrow = 9)
  inc[, 1] <- y0 / popsize()
  Nsus <- agedist() - y0 / popsize()
  
  whichper <- 1
  for(t in 1:Tmax) {
    if(t > tperiods[whichper]) {
      whichper <- whichper + 1
    }
    foi <- sapply(1:9, function(i) 
      sum(
        inc[, pmax(1, (t - 11):t)] *
          contacts[i, , whichper] *
          relsusinf *
          InfectivityProfile(infectivities[whichper]) *
          dailyfactors[t]) *
        relsusinf[i]
    )
    inc[, t + 1] <- Nsus * (1 - exp(- foi))
    Nsus <- Nsus - inc[, t + 1]
  }
  return(popsize() * inc[, -1])
}

## simulator for likelihood variables
engine_logLikvars <- function(parms, Tmax, tperiods, contacts, infdeviation, dailyfactors,
                              prob_list, delay_list,
                              maxdelay, tsero,
                              reportprob,
                              periodgroups = NULL,
                              estimate_relsusinf = TRUE, # if true, it is part of parms, otherwise provided separately
                              logrelsusinf = c()) {
  # transform relative susceptibility/infectivity
  if(estimate_relsusinf) {
    relsusinf <- c(1, exp(parms[2:9]))
  } else {
    relsusinf <- c(1, exp(logrelsusinf))
  }
  relsusinf <- relsusinf / (sum(relsusinf))
  
  # make infectivities vector
  infectivities <- exp(tail(parms, -1 - 8 * estimate_relsusinf))
  if(!is.null(periodgroups)) {
    infectivities <- infectivities[periodgroups]
  } 
  infectivities <- infectivities * infdeviation
  
  # simulate infection incidence
  incidences <- engine_incidence(Tmax = Tmax, 
                                 y0 = exp(parms[1]),
                                 tperiods = tperiods,
                                 contacts = contacts,
                                 infectivities = infectivities,
                                 dailyfactors = dailyfactors,
                                 relsusinf = relsusinf)

  # calculate incidence of severe cases (ready for hospital admission)
  tobesevere <- t(incidences * prob_list$pI2Se)
 
  incsevere <- sapply(1:Tmax,
                      function(x) colSums(tobesevere[max(1, x - maxdelay):x, , drop = FALSE] *
                                            delay_list$dI2Se[min(1 + maxdelay, x):1, x, ]
                      )) %>% t()
  
  # calculate incidence at IC
  togotoIC <- incsevere * prob_list$pSe2A2IC
  incatIC <- sapply(1:Tmax,
                     function(x) colSums(togotoIC[(max(1, x - maxdelay)):x, , drop = FALSE] *
                                           delay_list$dA2IC[min(1 + maxdelay, x):1, x])) %>% colSums

  # calculate seroprevalences
  seroprevs <- sapply(tsero, function(ts) rowSums(incidences[, 1:ts]) / agedist() / popsize())
  
  # return results in a list
  return(list(incatIC * reportprob, seroprevs))
}

## minus log-Likelihood calculator
engine_minlogLik <- function(parms, dataIC, dataseron, dataserox,
                         Tmax, tperiods, contacts, infdeviation, dailyfactors,
                         prob_list, delay_list,
                         maxdelay, tsero,
                         reportprob,
                         periodgroups = NULL,
                         estimate_relsusinf = TRUE,
                         logrelsusinf = c()) {
  # simulate expected values for observed variables
  expected <- engine_logLikvars(parms, Tmax, tperiods, contacts, infdeviation, dailyfactors,
                                prob_list, delay_list,
                                maxdelay, tsero,
                                reportprob,
                                periodgroups,
                                estimate_relsusinf,
                                logrelsusinf)
  
  # likelihood term for IC incidence
  toreturn <- -sum(dpois(dataIC, expected[[1]], log = TRUE)) 
  
  # likelihood term for seroprevalences
  if(length(tsero) > 0) {
    toreturn <- toreturn -
           sum(dbinom(dataserox, size = dataseron, prob = expected[[2]], log = TRUE))
  }
  
  # return likelihood
  return(toreturn)
}

## simulator for all variables
engine_allvars <- function(logy0, loginfectivities, logrelsusinf, 
                           Tmax, tperiods, contacts, infdeviation, dailyfactors,
                           prob_list, delay_list,
                           maxdelay, periodgroups = NULL) {
  # transform relative susceptibility/infectivity
  relsusinf <- c(1, exp(logrelsusinf))
  relsusinf <- relsusinf / (sum(relsusinf))

  # make infectivities vector
  infectivities <- exp(loginfectivities)
  if(!is.null(periodgroups)) {
    infectivities <- infectivities[periodgroups]
  }
  infectivities <- infectivities * infdeviation
  
  # simulate infection incidence
  incidences <- engine_incidence(Tmax = Tmax, 
                                 y0 = exp(logy0),
                                 tperiods = tperiods,
                                 contacts = contacts,
                                 infectivities = infectivities,
                                 dailyfactors = dailyfactors,
                                 relsusinf = relsusinf)
  
  # calculate incidence of symptom onset (assuming everyone is symptomatic)
  tobesymptomatic <- t(incidences)
  incsymptomatic <- sapply(1:Tmax,
                      function(x) colSums(tobesymptomatic[max(1, x - maxdelay):x, , drop = FALSE] *
                                            delay_list$dI2S[min(1 + maxdelay, x):1]
                      )) %>% t()
  
  # calculate incidence of severe cases (ready for hospital admission)
  tobesevere <- t(incidences * prob_list$pI2Se)
  incsevere <- sapply(1:Tmax,
                      function(x) colSums(tobesevere[max(1, x - maxdelay):x, , drop = FALSE] *
                                            delay_list$dI2Se[min(1 + maxdelay, x):1, x, ]
                      )) %>% t()
  
  # calculate incidence of hospital admissions
  inchospadmissions <- incsevere * prob_list$pSe2A

  # calculate incidence of hospital discharges
  tobedischargedfromhospital <- incsevere * prob_list$pSe2A2D
  inchospdischarges <- sapply(1:Tmax,
                              function(x) colSums(tobedischargedfromhospital[(max(1, x - maxdelay)):x, , drop = FALSE] *
                                                    delay_list$dA2D[min(1 + maxdelay, x):1, x])) %>% t()
  
  # calculate incidence of IC admissions
  togotoIC <- incsevere * prob_list$pSe2A2IC
  incICadmissions <- sapply(1:Tmax,
                     function(x) colSums(togotoIC[(max(1, x - maxdelay)):x, , drop = FALSE] *
                                           delay_list$dA2IC[min(1 + maxdelay, x):1, x])) %>% t()
  
  # calculate incidence of IC discharges
  tobedischargedfromIC <- incsevere * prob_list$pSe2A2IC2D
  incICdischarges <- sapply(1:Tmax,
                            function(x) colSums(tobedischargedfromIC[(max(1, x - maxdelay)):x, , drop = FALSE] *
                                                  delay_list$dA2IC2D[min(1 + maxdelay, x):1, x])) %>% t()
  
  # calculate incidence of re-admission to nursing ward
  tobereadmittedtohospital <- incsevere * prob_list$pSe2A2IC2H
  inchospreadmissions <- sapply(1:Tmax,
                                function(x) colSums(tobereadmittedtohospital[(max(1, x - maxdelay)):x, , drop = FALSE] *
                                                      delay_list$dA2IC2D[min(1 + maxdelay, x):1, x])) %>% t()
  
  # calculate incidence of final discharge from nursing ward (after IC)
  inchospfinaldischarges <- sapply(1:Tmax,
                                function(x) colSums(tobereadmittedtohospital[(max(1, x - maxdelay)):x, , drop = FALSE] *
                                                      delay_list$dA2IC2H2D[min(1 + maxdelay, x):1, x])) %>% t()
  
  # calculate prevalence on IC (adding discharges/readmissions to include discharge day in prevalence)
  prevIC <- apply(incICadmissions, 2, cumsum) -
    apply(incICdischarges, 2, cumsum) + incICdischarges - 
    apply(inchospreadmissions, 2, cumsum) + inchospreadmissions
  
  # calculate prevalence in hospital (adding discharges/readmissions to include discharge day in prevalence)
  prevhosp <- apply(inchospadmissions, 2, cumsum) -
    apply(inchospdischarges, 2, cumsum) + inchospdischarges -
    apply(incICdischarges, 2, cumsum) + incICdischarges - 
    apply(inchospfinaldischarges, 2, cumsum) + inchospfinaldischarges
  
  # calculate prevalence on nursing wards (adding discharges/readmissions to include discharge day in prevalence)
  prevnursing <- apply(inchospadmissions, 2, cumsum) +
    apply(inchospreadmissions, 2, cumsum) -
    apply(inchospdischarges, 2, cumsum) + inchospdischarges -
    apply(incICadmissions, 2, cumsum) + incICadmissions - 
    apply(inchospfinaldischarges, 2, cumsum) + inchospfinaldischarges
  
  # calculate prevalence of recovered individuals
  previmmune <- rbind(matrix(0, ncol = 9, nrow = delayI2R()), 
                      apply(head(t(incidences), -delayI2R()), 2, cumsum))

  # make tibble to return
  # combine and add names
  ageclassnames <- paste0("[", seq(0, 80, 10), ",", c(seq(10, 80, 10), "Inf"), ")")
  toreturn <- as_tibble(t(incidences), .name_repair = ~ ageclassnames) %>%
    mutate(time = 1:Tmax) %>%
    pivot_longer(1:9, names_to = "ageclass", values_to = "incinfection") %>%
    full_join(
      as_tibble(incsymptomatic, .name_repair = ~ ageclassnames) %>%
        mutate(time = 1:Tmax) %>%
        pivot_longer(1:9, names_to = "ageclass", values_to = "incsymp"),
      by = c("time", "ageclass")
    ) %>%
    full_join(
      as_tibble(incsevere, .name_repair = ~ ageclassnames) %>%
        mutate(time = 1:Tmax) %>%
        pivot_longer(1:9, names_to = "ageclass", values_to = "incsevere"),
      by = c("time", "ageclass")
    ) %>%
    full_join(
      as_tibble(inchospadmissions, .name_repair = ~ ageclassnames) %>%
        mutate(time = 1:Tmax) %>%
        pivot_longer(1:9, names_to = "ageclass", values_to = "inchosp"),
      by = c("time", "ageclass")
    ) %>%
    full_join(
      as_tibble(incICadmissions, .name_repair = ~ ageclassnames) %>%
        mutate(time = 1:Tmax) %>%
        pivot_longer(1:9, names_to = "ageclass", values_to = "incICU"),
      by = c("time", "ageclass")
    ) %>%
    full_join(
      as_tibble(prevhosp, .name_repair = ~ ageclassnames) %>%
        mutate(time = 1:Tmax) %>%
        pivot_longer(1:9, names_to = "ageclass", values_to = "prevhosp"),
      by = c("time", "ageclass")
    ) %>%
    full_join(
      as_tibble(prevIC, .name_repair = ~ ageclassnames) %>%
        mutate(time = 1:Tmax) %>%
        pivot_longer(1:9, names_to = "ageclass", values_to = "prevICU"),
      by = c("time", "ageclass")
    ) %>%
    full_join(
      as_tibble(prevnursing, .name_repair = ~ ageclassnames) %>%
        mutate(time = 1:Tmax) %>%
        pivot_longer(1:9, names_to = "ageclass", values_to = "prevnursing"),
      by = c("time", "ageclass")
    ) %>%
    full_join(
      as_tibble(previmmune, .name_repair = ~ ageclassnames) %>%
        mutate(time = 1:Tmax) %>%
        pivot_longer(1:9, names_to = "ageclass", values_to = "previmmune"),
      by = c("time", "ageclass")
    )  
  
  return(toreturn)
}

###########################################
### functions for user-friendlier input ###
###########################################
logLik_simulate <- function(logy0, 
                            loginfectivities,
                            logrelsusinf = "default",
                            probi2se = "default",
                            enddates,
                            contactcontrol,
                            seasonality = "None",
                            reportprob = 1,
                            periodgroups = NULL,
                            Tmax = max(enddates),
                            samplenr = 1,
                            maxdelay = 30,
                            serodates = NULL) {
  parms <- c(LogY0(logy0 = logy0, nr = samplenr), 
             LogRelativeSusInf(logrelsusinf = logrelsusinf, nr = samplenr), 
             LogInfectivities(loginfectivities = loginfectivities, nr = samplenr))
  Tmax <- as.numeric(Tmax - as.Date("2020-02-12"))
  tperiods <- as.numeric(enddates - as.Date("2020-02-12"))
  contacts <- create_contactarray(contactcontrol, samplenr)
  infdeviation <- create_infdeviations(contactcontrol, samplenr)
  dailyfactors <- create_dailyfactors(seasonality, Tmax)
  
  prob_list <- customise_probs(Tmax, probi2se, nr = samplenr, likelihood_only = TRUE)
  delay_list <- customise_delays(Tmax, maxdelay, likelihood_only = TRUE)
  
  if(!is.null(serodates) & length(serodates) > 0) {
    tsero <- as.numeric(serodates - as.Date("2020-02-12") - delayI2R())
  } else {
    tsero <- c()
  }
  
  reportprob <- tail(c(rep(1, Tmax), rev(reportprob)), Tmax)
  
  engine_logLikvars(parms, Tmax, tperiods, contacts, infdeviation, dailyfactors,
                    prob_list, delay_list,
                    maxdelay, tsero, reportprob, periodgroups)
}

logLik_calculate <- function(logy0, 
                             loginfectivities,
                             logrelsusinf = "default",
                             probi2se = "default",
                             dataIC, dataseron = c(), dataserox = c(),
                             enddates,
                             contactcontrol,
                             seasonality = "None",
                             reportprob = 1,
                             periodgroups = NULL,
                             Tmax = max(enddates),
                             samplenr = 1,
                             maxdelay = 30,
                             serodates = NULL) {
  parms <- c(LogY0(logy0 = logy0, nr = samplenr), 
             LogRelativeSusInf(logrelsusinf = logrelsusinf, nr = samplenr), 
             LogInfectivities(loginfectivities = loginfectivities, nr = samplenr))
  Tmax <- as.numeric(Tmax - as.Date("2020-02-12"))
  tperiods <- as.numeric(enddates - as.Date("2020-02-12"))
  contacts <- create_contactarray(contactcontrol, samplenr)
  infdeviation <- create_infdeviations(contactcontrol, samplenr)
  dailyfactors <- create_dailyfactors(seasonality, Tmax)
  
  prob_list <- customise_probs(Tmax, probi2se, nr = samplenr, likelihood_only = TRUE)
  delay_list <- customise_delays(Tmax, maxdelay, likelihood_only = TRUE)
  
  if(!is.null(serodates) & length(serodates) > 0) {
    tsero <- as.numeric(serodates - as.Date("2020-02-12") - delayI2R())
  } else {
    tsero <- c()
  }
  
  reportprob <- tail(c(rep(1, Tmax), rev(reportprob)), Tmax)
  
  engine_minlogLik(parms, dataIC, dataseron, dataserox, 
                   Tmax, tperiods, contacts, infdeviation, dailyfactors,
                   prob_list, delay_list,
                   maxdelay, tsero, reportprob, periodgroups, logrelsusinf = logrelsusinf)
}

logLik_optimise <- function(logy0, 
                            loginfectivities,
                            logrelsusinf = "default",
                            probi2se = "default",
                            dataIC, dataseron = c(), dataserox = c(),
                            enddates,
                            contactcontrol,
                            seasonality = "None",
                            reportprob = 1,
                            periodgroups = NULL,
                            Tmax = max(enddates),
                            samplenr = 1,
                            maxdelay = 30,
                            serodates = NULL,
                            estimate_relsusinf = TRUE, ...) {
  if(estimate_relsusinf) {
    parms <- c(LogY0(logy0 = logy0, nr = samplenr), 
               LogRelativeSusInf(logrelsusinf = logrelsusinf, nr = samplenr), 
               LogInfectivities(loginfectivities = loginfectivities, nr = samplenr))
  } else {
    parms <- c(LogY0(logy0 = logy0, nr = samplenr), 
               LogInfectivities(loginfectivities = loginfectivities, nr = samplenr))
  }
  Tmax <- as.numeric(Tmax - as.Date("2020-02-12"))
  tperiods <- as.numeric(enddates - as.Date("2020-02-12"))
  contacts <- create_contactarray(contactcontrol, samplenr)
  infdeviation <- create_infdeviations(contactcontrol, samplenr)
  dailyfactors <- create_dailyfactors(seasonality, Tmax)
  
  prob_list <- customise_probs(Tmax, probi2se, nr = samplenr, likelihood_only = TRUE)
  delay_list <- customise_delays(Tmax, maxdelay, likelihood_only = TRUE)
  
  if(!is.null(serodates) & length(serodates) > 0) {
    tsero <- as.numeric(serodates - as.Date("2020-02-12") - delayI2R())
  } else {
    tsero <- c()
  }
  
  reportprob <- tail(c(rep(1, Tmax), rev(reportprob)), Tmax)
  
  toreturn <- optim(parms, engine_minlogLik, dataIC = dataIC, dataseron = dataseron, dataserox = dataserox,
        Tmax = Tmax, tperiods = tperiods,
        contacts = contacts,  infdeviation = infdeviation, dailyfactors = dailyfactors, prob_list = prob_list, delay_list = delay_list,
        maxdelay = maxdelay, tsero = tsero, reportprob = reportprob, periodgroups = periodgroups,
        estimate_relsusinf = estimate_relsusinf, 
        logrelsusinf = LogRelativeSusInf(logrelsusinf = logrelsusinf, nr = samplenr), ...)
  
  return(toreturn)
}

simulate_single <- function(logy0, 
                            loginfectivities,
                            logrelsusinf = "default",
                            probi2se = "default",
                            enddates,
                            contactcontrol,
                            seasonality = "None",
                            periodgroups = NULL,
                            Tmax = max(enddates),
                            samplenr = 1,
                            maxdelay = 30) {
  logy0 <- LogY0(logy0 = logy0, nr = samplenr)
  loginfectivities <- LogInfectivities(loginfectivities = loginfectivities, nr = samplenr)
  logrelsusinf <- LogRelativeSusInf(logrelsusinf = logrelsusinf, nr = samplenr)
  cat(logrelsusinf)
  Tmax <- as.numeric(Tmax - as.Date("2020-02-12"))
  tperiods <- as.numeric(enddates - as.Date("2020-02-12"))
  contacts <- create_contactarray(controlregimes = contactcontrol, nr = samplenr)
  infdeviation <- create_infdeviations(controlregimes = contactcontrol, nr = samplenr)
  dailyfactors <- create_dailyfactors(seasonality, Tmax)
  
  prob_list <- customise_probs(Tmax, probi2se, nr = samplenr, likelihood_only = FALSE)
  delay_list <- customise_delays(Tmax, maxdelay, likelihood_only = FALSE)
  
  engine_allvars(logy0, loginfectivities, logrelsusinf, 
                 Tmax, tperiods, contacts, infdeviation, dailyfactors,
                 prob_list, delay_list,
                 maxdelay, periodgroups)
}

# multiple scenarios with the same parameter values: scenarios differ by contactcontrol and/or periodgroups,
# which are provided as lists
simulate_multiple <- function(logy0, 
                              loginfectivities,
                              logrelsusinf = "default",
                              probi2se = "default",
                              enddates,
                              contactcontrol,
                              seasonality = "None",
                              scenarionames = NULL,
                              periodgroups = NULL,
                              Tmax = max(enddates),
                              samplenr = 1,
                              maxdelay = 30) {
  logy0 <- LogY0(logy0 = logy0, nr = samplenr)
  loginfectivities <- LogInfectivities(loginfectivities = loginfectivities, nr = samplenr)
  logrelsusinf <- LogRelativeSusInf(logrelsusinf = logrelsusinf, nr = samplenr)
  Tmax <- as.numeric(Tmax - as.Date("2020-02-12"))
  tperiods <- as.numeric(enddates - as.Date("2020-02-12"))
  dailyfactors <- create_dailyfactors(seasonality, Tmax)
  
  prob_list <- customise_probs(Tmax, probi2se, nr = samplenr, likelihood_only = FALSE)
  delay_list <- customise_delays(Tmax, maxdelay, likelihood_only = FALSE)
  
  if(!is.list(periodgroups) & !is.list(contactcontrol)) {
    periodgroups <- list(periodgroups)
    contactcontrol <- list(contactcontrol)
  }
  if(!is.list(contactcontrol)) {
    contactcontrol <- rep(list(contactcontrol), length(periodgroups))
  }
  if(!is.list(periodgroups)) {
    periodgroups <- rep(list(periodgroups), length(contactcontrol))
  }
  if(is.null(scenarionames)) {
    scenarionames <- paste0("Scenario_", 1:length(contactcontrol))
  }
  
  contacts <- create_contactarray(contactcontrol[[1]], samplenr)
  infdeviation <- create_infdeviations(contactcontrol[[1]], nr = samplenr)
  toreturn <- engine_allvars(logy0, loginfectivities, logrelsusinf, 
                             Tmax, tperiods, contacts,
                             prob_list, delay_list,
                             maxdelay, periodgroups[[1]]) %>%
    mutate(scenario = scenarionames[1])
  
  if(length(contactcontrol) > 1) {
    for(i in 2:length(contactcontrol)) {
      contacts <- create_contactarray(contactcontrol[[i]], samplenr)
      infdeviation <- create_infdeviations(contactcontrol[[i]], nr = samplenr)
      toreturn <- bind_rows(toreturn, 
                            engine_allvars(logy0, loginfectivities, logrelsusinf, 
                                           Tmax, tperiods, contacts, infdeviation, dailyfactors,
                                           prob_list, delay_list,
                                           maxdelay, periodgroups[[i]]) %>%
                              mutate(scenario = scenarionames[i]))
    }
  }
  
  return(toreturn)
}

simulate_replicates <- function(logy0 = "sample", 
                                loginfectivities = "sample",
                                logrelsusinf = "default",
                                probi2se = "default",
                                enddates,
                                contactcontrol,
                                seasonality = "None",
                                scenarionames = NULL,
                                periodgroups = NULL,
                                Tmax = max(enddates),
                                nrof_samples = 200,
                                maxdelay = 30) {
  Tmax <- as.numeric(Tmax - as.Date("2020-02-12"))
  tperiods <- as.numeric(enddates - as.Date("2020-02-12"))
  
  delay_list <- customise_delays(Tmax, maxdelay, likelihood_only = FALSE)
  
  if(!is.list(periodgroups) & !is.list(contactcontrol)) {
    periodgroups <- list(periodgroups)
    contactcontrol <- list(contactcontrol)
  }
  if(!is.list(contactcontrol)) {
    contactcontrol <- rep(list(contactcontrol), length(periodgroups))
  }
  if(!is.list(periodgroups)) {
    periodgroups <- rep(list(periodgroups), length(contactcontrol))
  }
  if(is.null(scenarionames)) {
    scenarionames <- paste0("Scenario_", 1:length(contactcontrol))
  }
  if(!is.list(seasonality)) {
    seasonality <- rep(list(seasonality), length(contactcontrol))
  }
  
  toreturn <- NULL
  for(repl in 1:nrof_samples) {
    sampledlogy0 <- LogY0(logy0 = logy0, nr = repl)
    sampledloginfectivities <- LogInfectivities(loginfectivities = loginfectivities, nr = repl)
    sampledlogrelsusinf <- LogRelativeSusInf(logrelsusinf = logrelsusinf, nr = repl)
    prob_list <- customise_probs(Tmax, probi2se = probi2se, nr = repl, likelihood_only = FALSE)
    
    for(i in 1:length(contactcontrol)) {
      contacts <- create_contactarray(contactcontrol[[i]], nr = repl)
      infdeviation <- create_infdeviations(contactcontrol[[i]], nr = repl)
      dailyfactors <- create_dailyfactors(seasonality[[i]], Tmax)
      toreturn <- bind_rows(toreturn, 
                            engine_allvars(sampledlogy0, sampledloginfectivities, sampledlogrelsusinf, 
                                           Tmax, tperiods, contacts, infdeviation, dailyfactors,
                                           prob_list, delay_list,
                                           maxdelay, periodgroups[[i]]) %>%
                              mutate(scenario = scenarionames[i], 
                                     repl = repl))
    }
    
  }
  
  return(toreturn)
}

