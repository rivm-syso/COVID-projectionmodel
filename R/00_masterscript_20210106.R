###########################################################
### Supplementary code with Klinkenberg et al (2023)
### Analyse data and simulate projections as on January 6 2021 
###########################################################

options(tidyverse.quiet = TRUE)
library(tidyverse)
library(lubridate)
library(deSolve)
options(dplyr.summarise.inform = FALSE)

#############################
### Data analyses: step 1 ###
#############################
ANALYSISDATE <- as.Date("2021-01-06")
DELAYSTARTDATE <- as.Date("2020-09-01")
source("R/code4model/Populationdata4model.R")
source("R/code4model/EPI_NICEanalyses4model.R") # warnings can be ignored
source("R/code4model/EPI_Reportingdelays4model.R")
source("R/code4model/OSIRISanalyses4model.R")
source("R/code4model/SEROanalysesContactinput4model.R")
source("R/code4model/readmatrices4model.R") # is in fact step 2
source("R/code4model/Seasonality4model.R") # not yet implemented on 6 Jan 2021
# image saved so that these estimates can be used later, with more recent data in step 3
rm(dataOsiris) # reduce file size
save(list = ls(), file = paste0("data/processedmodelinput/modelinput", ANALYSISDATE, ".RData"))
### After downloading the repo, use the line below to unzip RData files
# unzip("data/processedmodelinput.zip")

#####################################
### Functions for the simulations ###
#####################################
source("R/code4model/ContactsInfectivities4model.R")
source("R/code4model/DelaysProbabilities4model.R")
source("R/code4model/simulationcode4model_v3.R")
source("R/code4model/SimulatePlot4model.R")

############################################
### incidence data for estimation step 3 ###
############################################
# reload data: possible to use with older results from step 1
############################################
ANALYSISDATE <- as.Date("2021-01-06")
source("R/code4model/EPI_NICEdata4fit.R") 
# save complete data and parameter values prior to step 3
save(list = ls(), file = paste0("data/processedmodelinput/modeldatainput", ANALYSISDATE, ".RData"))
### create zip files to reduce size
zip("data/processedmodelinput.zip",
    c(paste0("data/processedmodelinput/modelinput", ANALYSISDATE, ".RData"),
      paste0("data/processedmodelinput/modeldatainput", ANALYSISDATE, ".RData")))

###########
### replace results from synthetic data by actual results from 6 Jan 2021 (datasets not replaced)
###########
source("R/code4model/Replace_syntheticresults.R")



###############
# optimise likelihood (step 3 of data analysis)
###############
### NB: many options were developed and tested and optimised on days that projections were not made
### Also many change point options for the relative infectivity were tried every week.
### This is only the estimation in the final step
llo <- logLik_optimise(logy0 = 1.68,
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

### Some alternative change point options, potentially needed in coming weeks
# periodgroups = c(1,2,3,3, 3,3,3,4, 5,6,7,8, 9,9) ; minloglik = 822.99
# periodgroups = c(1,2,3,3, 3,3,3,4, 3,5,5,6, 6,6) ; minloglik = 824.84  => final model
# periodgroups = c(1,2,3,3, 3,3,3,4, 3,5,5,6, 7,7) ; minloglik = 822.97  

saveRDS(llo, "results/maxlikelihood_20210106.rds")

####################
### Make 200 samples for parameters, from point estimates with Hessian matrix
####################
samplesforsimulation(llo)

### run all scenarios and prognoses: as an example the possibility of relaxing control up to the level of November (partial lockdown)
simrepls <- simulate_replicates(logrelsusinf = "default",
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

saveRDS(simrepls, "results/simulations_20210106.rds")
# save only 10% to reduce file size for github
saveRDS(simrepls %>% filter(repl <= 20), "results/simulations_10pct_20210106.rds")



############ PROGNOSES ##################
### plot results IC and HOSP, with original Dutch captions
#########################################
plotICinc(simrepls, from = as.Date("2020-12-01"), to = as.Date("2021-05-01"), scenarios = paste0("Scenario_", c(1, 2)),
          scenarionames_in_plot = c(Scenario_1 = "met gevolgde\n maatregelen",
                                    Scenario_2 = "afschalen per\n 19 januari"))

plotICprev(simrepls, from = as.Date("2020-12-01"), to = as.Date("2021-05-01"), scenarios = paste0("Scenario_", c(1, 2)),
          scenarionames_in_plot = c(Scenario_1 = "met gevolgde\n maatregelen",
                                    Scenario_2 = "afschalen per\n 19 januari"))

plothospinc(simrepls, from = as.Date("2020-12-01"), to = as.Date("2021-05-01"), scenarios = paste0("Scenario_", c(1, 2)),
          scenarionames_in_plot = c(Scenario_1 = "met gevolgde\n maatregelen",
                                    Scenario_2 = "afschalen per\n 19 januari"))

plothospprev(simrepls, from = as.Date("2020-12-01"), to = as.Date("2021-05-01"), scenarios = paste0("Scenario_", c(1, 2)),
          scenarionames_in_plot = c(Scenario_1 = "met gevolgde\n maatregelen",
                                    Scenario_2 = "afschalen per\n 19 januari"))



