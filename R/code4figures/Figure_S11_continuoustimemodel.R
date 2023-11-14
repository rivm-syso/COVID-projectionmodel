###########################################################
### Supplementary code with Klinkenberg et al (2023)
### Additional analysis: fit ODE model to data from 6 January 2021
###########################################################

################################################################
### Load all results until model fitting (from Masterscript) ###
################################################################
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(lubridate)
library(deSolve)
library(cowplot)
options(dplyr.summarise.inform = FALSE)

ANALYSISDATE <- as.Date("2021-01-06")
load(paste0("data/processedmodelinput/modeldatainput", ANALYSISDATE, ".RData"))
source("R/code4model/Replace_syntheticresults.R")
llo_discretetime <- readRDS("results/maxlikelihood_20210106.rds")

##############################################
### Define functions to simulate ODE model ###
##############################################
# compartments: S, E, I, and cumulative incidence
ODEcompartments <- function() return(c("S", "EA", "EB", "IA", "IB", "cumI"))

# natural history parameters (SI = serial interval)
parmsNatHist <- function(SI = 4) {
  return(list(beta = 3.5 / SI / 2,
              gamma = 3.5 / SI))
}

### rate of new infections into different age classes 
# y: current state of the ODE model (named vector)
# contacts: contact matrix
# susceptibility, infectiousness: vectors with relative susc/inf by age
# parms: list with parameter values
covidincidence <- function(y, contacts, susceptibility, infectiousness, parms) {
  res <- rep(0, 9)
  for(i in ODEcompartments()) {
    if(substr(i, 1, 1) == "I") {
      res <- res + parms[["beta"]] * 
        susceptibility *
        colSums(y[paste0(i, 1:9)] * infectiousness * contacts) * y[paste0("S", 1:9)]
    }
  }
  names(res) <- paste0("la", 1:9)
  return(res)
}

### model equations of SEEIIR model: S, I, and cumulative incidence
# arguments as required by ode() function of deSolve
covidODEmodel <- function(t, y, parms) {
  incid <- covidincidence(y, parms$contacts, parms$susceptibility, parms$infectiousness, parms)
  
  res <- rep(0, length(y))
  names(res) <- names(y)
  
  # dS/dt
  for(i in 1:9) {
    res[paste0("S", i)] <- -incid[paste0("la", i)]
  }
  
  # dEA/dt
  for(i in 1:9) {
    res[paste0("EA", i)] <- incid[paste0("la", i)] - parms[["gamma"]] * y[paste0("EA", i)]
  }
  
  # dEB/dt
  for(i in 1:9) {
    res[paste0("EB", i)] <- parms[["gamma"]] * y[paste0("EA", i)] - parms[["gamma"]] * y[paste0("EB", i)]
  }
  
  # dIA/dt
  for(i in 1:9) {
    res[paste0("IA", i)] <- parms[["gamma"]] * y[paste0("EB", i)] - parms[["gamma"]] * y[paste0("IA", i)]
  }
  
  # dIB/dt
  for(i in 1:9) {
    res[paste0("IB", i)] <- parms[["gamma"]] * y[paste0("IA", i)] - parms[["gamma"]] * y[paste0("IB", i)]
  }
  
  # dcumI/dt
  for(i in 1:9) {
    res[paste0("cumI", i)] <- incid[paste0("la", i)]
  }
  
  # finished
  return(list(res))
}


### incidence simulator: this function replaces the discrete-time equivalent originally used on 6 January 2021
# first keep the original function in memory
engine_incidence_discretetime <- engine_incidence
engine_incidence_continuoustime <- function(Tmax, # length of simulation
                                            y0, # length-9 vector of initial incidence (12-day long from 1 Feb to 12 Feb)
                                            tperiods, # vector with final days of periods with different contacts and/or infectivities
                                            contacts, # 9*9*tperiods array with contact matrices, one for each period
                                            infectivities, # vector with infectivities, one for each period
                                            dailyfactors, # Tmax-length vector of multiplication factors to reflect seasonality
                                            relsusinf) { # normalised length-9 vector with relative susceptibilities/infectivities
  ### define all compartments
  compartments <- t(outer(ODEcompartments(), 1:9, paste0))
  
  ### initialize compartments
  initialstate <- rep(0, length(compartments))
  names(initialstate) <- compartments
  imports <- y0 / popsize()
  initialstate[paste0("S", 1:9)] <- agedist() - imports
  initialstate[paste0("EA", 1:9)] <- imports / 4
  initialstate[paste0("EB", 1:9)] <- imports / 4
  initialstate[paste0("IA", 1:9)] <- imports / 4
  initialstate[paste0("IA", 1:9)] <- imports / 4
  
  ### simulate no control
  nathistparms2use <- parmsNatHist()
  nathistparms2use$beta <- nathistparms2use$beta * infectivities[1]
  ressimul <- ode(
    func = covidODEmodel,
    y = initialstate,
    times = seq(0, tperiods[1], 1),
    parms = c(list(contacts = contacts[,,1],
                   susceptibility = relsusinf,
                   infectiousness = relsusinf),
              nathistparms2use),
    method = "lsoda"
  )
  
  ### simulate controls
  for(i in 2:length(tperiods)) {
    nathistparms2use <- parmsNatHist()
    nathistparms2use$beta <- nathistparms2use$beta * infectivities[i]
    ressimul <- rbind(
      head(ressimul, -1),
      ode(
        func = covidODEmodel,
        y = tail(ressimul, 1)[, compartments],
        times = seq(tperiods[i - 1], tperiods[i], 1),
        parms = c(list(contacts = contacts[,,i],
                       susceptibility = relsusinf,
                       infectiousness = relsusinf),
                  nathistparms2use),
        method = "lsoda"
      )
    )
  }
  
  return(t(popsize() * diff(ressimul[, paste0("cumI", 1:9)])))
}

  
#################################
### fit continuous-time model ###
#################################
engine_incidence <- engine_incidence_continuoustime
llo_continuoustime <- logLik_optimise(logy0 = 3.4,
                           loginfectivities = c(5.23,4.71,5.20,5.38,5.27,5.51),
                           logrelsusinf = "default",
                           probi2se = "default",
                           dataIC = inc_nice_ic(),
                           enddates = c(as.Date(c("2020-03-18", "2020-03-28", "2020-05-10", "2020-06-01",
                                                  "2020-07-05", "2020-08-30", "2020-09-28", "2020-10-14",
                                                  "2020-10-25", "2020-11-04", "2020-11-18", "2020-12-14",
                                                  "2020-12-20")), ANALYSISDATE),
                           contactcontrol = c("precontrol_mean", "intellockdown_mean", "intellockdown_mean", "batch1_mean",
                                              "batch2_mean", "batch3_mean", "sep2020_mean", "okt2020_mean",
                                              "okt2020holiday_mean", "partlockdown_mean", "nov2020_mean", "partlockdown_mean",
                                              "winterlockdown_mean", "winterlockdownChristmas2020_mean"),
                           reportprob = ReportingDelays$delayIC,
                           periodgroups = c(1,2,3,3, 3,3,3,4, 3,5,5,6, 6,6),
                           estimate_relsusinf = F,
                           hessian = T, maxdelay = 30)
saveRDS(llo_continuoustime, "results/additionalresults/maxlikelihoodcontinuoustime_20210106.rds")

####################
### Make 200 samples for parameters, from point estimates with Hessian matrix
####################
samplesforsimulation(llo_continuoustime)

### run all scenarios and prognoses: as an example the possibility of relaxing control up to the level of November (partial lockdown)
simrepls_continuoustime <- simulate_replicates(logrelsusinf = "default",
                                enddates = c(as.Date(c("2020-03-18", "2020-03-28", "2020-05-10", "2020-06-01",
                                                       "2020-07-05", "2020-08-30", "2020-09-28", "2020-10-14",
                                                       "2020-10-25", "2020-11-04", "2020-11-18", "2020-12-14", 
                                                       "2020-12-20")), ANALYSISDATE - 14, 
                                             as.Date(c("2021-01-03", "2021-01-18", "2021-08-01"))),
                                contactcontrol = 
                                  list(c("precontrol_mean", "intellockdown_mean", "intellockdown_mean", "batch1_mean",
                                         "batch2_mean", "batch3_mean", "sep2020_mean", "okt2020_mean",
                                         "okt2020holiday_mean", "partlockdown_mean", "nov2020_mean", "partlockdown_mean", 
                                         "winterlockdown_mean", "winterlockdownChristmas2020_mean", 
                                         "winterlockdownChristmas2020", "winterlockdown", "winterlockdown"), # default
                                       c("precontrol_mean", "intellockdown_mean", "intellockdown_mean", "batch1_mean",
                                         "batch2_mean", "batch3_mean", "sep2020_mean", "okt2020_mean",
                                         "okt2020holiday_mean", "partlockdown_mean", "nov2020_mean", "partlockdown_mean", 
                                         "winterlockdown_mean", "winterlockdownChristmas2020_mean", 
                                         "winterlockdownChristmas2020", "winterlockdown", "partlockdown") # opheffen 19 januari
                                  ), 
                                periodgroups = list(c(1,2,3,3, 3,3,3,4, 3,5,5,6, 6,6, 6,6,6),
                                                    c(1,2,3,3, 3,3,3,4, 3,5,5,6, 6,6, 6,6,6)), 
                                seasonality = "None",
                                nrof_samples = 200,
                                maxdelay = 100)

saveRDS(simrepls_continuoustime, "results/additionalresults/simulations_continuoustime_20210106.rds")
# save only 10% to reduce file size for github
saveRDS(simrepls_continuoustime %>% filter(repl <= 20), "results/additionalresults/simulations_continuoustime_10pct_20210106.rds")

simrepls_discretetime <- readRDS("results/simulations_20210106.rds")

############ PROGNOSES ##################
### plot results IC and HOSP, with original Dutch captions
#########################################
pICinc_d <- plotICinc(simrepls_discretetime, from = as.Date("2020-02-13"), to = as.Date("2021-05-01"), scenarios = paste0("Scenario_", c(1, 2)),
                      scenarionames_in_plot = c(Scenario_1 = "control measures\n unchanged",
                                                Scenario_2 = "control measures back\n to November levels\n on 19 January")) + 
  labs(title = "Discrete-time model: number of ICU admissions per day", x = "Day", y = "Daily ICU admissions") +
  scale_shape(labels = "NICE ICU admissions", name = "data points")

pICinc_c <- plotICinc(simrepls_continuoustime, from = as.Date("2020-02-13"), to = as.Date("2021-05-01"), scenarios = paste0("Scenario_", c(1, 2)),
                      scenarionames_in_plot = c(Scenario_1 = "control measures\n unchanged",
                                                Scenario_2 = "control measures back\n to November levels\n on 19 January")) + 
  labs(title = "Continuous-time model: number of ICU admissions per day", x = "Day", y = "Daily ICU admissions") +
  scale_shape(labels = "NICE ICU admissions", name = "data points")

pICinc_diff <- 
  simrepls_continuoustime %>%
  select(time, ageclass, incICU, scenario, repl) %>%
  mutate(model = "continuous") %>%
  bind_rows(simrepls_discretetime %>%
              select(time, ageclass, incICU, scenario, repl) %>%
              mutate(model = "discrete")) %>%
  group_by(scenario, repl, time, model) %>%
  reframe(incICU = sum(incICU)) %>%
  group_by(time, scenario, repl) %>%
  reframe(difference = incICU[model == "continuous"] - incICU[model == "discrete"]) %>%
  group_by(time, scenario) %>%
  reframe(diffmean = mean(difference)) %>%
  ungroup() %>%
  mutate(time = time + as.Date("2020-02-12")) %>%
  filter(time <= as.Date("2021-05-01")) %>%
  ggplot(aes(x = time, y = diffmean, color = scenario)) +
  geom_line() + 
  theme_light() +
  scale_x_date(date_breaks = "1 month", date_labels = "1 %b") +
  scale_color_discrete(labels = c(Scenario_1 = "control measures\n unchanged",
                                  Scenario_2 = "control measures back\n to November levels\n on 19 January")) +
  labs(title = "Difference between continuous-time and discrete-time model projections",
       x = "Day", y = "Difference in mean number of ICU admissions")

plot_grid(pICinc_c, pICinc_d, pICinc_diff, nrow = 3)
ggsave("results/additionalresults/FigS11.pdf", bg = "white", width = 20, height = 25, units = "cm")
ggsave("results/additionalresults/FigS11.jpg", bg = "white", width = 20, height = 25, units = "cm", dpi = 600)


pICinc_d +
  labs(title = "Six-week model projection of 6 January 2021", x = "Day", y = "Daily ICU admissions") +
  scale_x_date(limits = as.Date(c("2020-02-01", "2021-02-21"))) +
  scale_shape(labels = "ICU admission data", name = "")
ggsave("results/additionalresults/FigPoster2.jpg", bg = "white", width = 20, height = 8, units = "cm", dpi = 600)

####################################
### Results Supplementary Tables ###
####################################
### confidence intervals
apply(MASS::mvrnorm(10000, mu = llo_continuoustime$par, Sigma = solve(llo_continuoustime$hessian)), 2, quantile, probs = c(.025,.5,.975))
apply(MASS::mvrnorm(10000, mu = llo_discretetime$par, Sigma = solve(llo_discretetime$hessian)), 2, quantile, probs = c(.025,.5,.975))

### correlation matrices
cov2cor(solve(llo_discretetime$hessian))
cov2cor(solve(llo_continuoustime$hessian))
