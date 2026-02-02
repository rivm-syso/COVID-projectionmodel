############################################################
# Function to determine reporting delays in NICE
############################################################


get_reportingdelay_NICE <- function(start_date = as.Date("2020-03-28"), end_date = today(), IC = TRUE) {
  if(exists("DELAYSTARTDATE")) {
    start_date <- DELAYSTARTDATE
  }
  if(exists("ANALYSISDATE")) {
    end_date <- ANALYSISDATE
  }
  
  temp <- data_nice %>%
    filter(IsIc == as.numeric(IC)) %>%
    group_by(Pid) %>%
    filter(min(AdmissionDate) >= start_date & min(AdmissionDate) <= end_date) %>%
    summarise(admission = min(AdmissionDate),
              reporting = min(as.Date(CreationDateNice)),
              delay = 1 + as.numeric(reporting - admission)) %>%
    mutate(admissionday = weekdays(admission))
  temp <- sapply(weekdays(7:1 + end_date),
                 function(x) temp %>%
                   filter(admissionday == x) %>%
                   pull(delay) %>%
                   tabulate(nbins = 101))
  temp <- apply(temp,2, function(x) cumsum(x)/sum(x))
  toreturn <- temp[cbind(1:101, rep(1:7, 15)[1:101])]
  
  
  return(toreturn)
}


ReportingDelays <- 
  list(
    delayIC = get_reportingdelay_NICE(),
    delayHOSP = get_reportingdelay_NICE(IC = FALSE)
  )

rm(get_reportingdelay_NICE)
