##################
### replace synthetic results by actual results
#################

#################
NICEtimeseries <- readRDS("data/originalresults/NICEtimeseries_20210106.rds")
NICEdelays <- readRDS("data/originalresults/NICEdelays_20210106.rds")
NICEprobabilities <- readRDS("data/originalresults/NICEprobabilities_20210106.rds")
PreAdmissionDelays <- readRDS("data/originalresults/PreAdmissionDelays_20210106.rds")
PreAdmissionProbs <- readRDS("data/originalresults/PreAdmissionProbs_20210106.rds")
ReportingDelays <- readRDS("data/originalresults/ReportingDelays_20210106.rds")


#################
inc_nice_ic_orig <- inc_nice_ic
prev_nice_ic_orig <- prev_nice_ic
inc_nice_hosp_orig <- inc_nice_hosp
prev_nice_hosp_orig <- prev_nice_hosp

inc_nice_ic <- function(ages = FALSE, ROAZ = "all", includeblacklisted = TRUE, original = TRUE) {
  if(original) return(NICEtimeseries$inc_ic)
  
  return(inc_nice_ic_orig(ages, ROAZ, includeblacklisted))
}

prev_nice_ic <- function(ages = FALSE, ROAZ = "all", original = TRUE) {
  if(original) return(NICEtimeseries$prev_ic)
  
  return(prev_nice_ic_orig(ages, ROAZ))
}

inc_nice_hosp <- function(ages = FALSE, ROAZ = "all", excludeIC = FALSE, original = TRUE) {
  if(original) return(NICEtimeseries$inc_hosp)
  
  return(inc_nice_hosp_orig(ages, ROAZ, excludeIC))
}

prev_nice_hosp <- function(ages = FALSE, ROAZ = "all", maxstay = Inf, original = TRUE) {
  if(original) return(NICEtimeseries$prev_hosp)
  
  return(prev_nice_hosp_orig(ages, ROAZ, maxstay))
}

