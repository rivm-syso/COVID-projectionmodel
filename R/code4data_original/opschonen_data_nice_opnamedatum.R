
# Make a data table to define First admission date; last discharge; first report date for...
  # Hospital admission
    # Episode 
    # Pid
  # IC admission
    # Episode 
    # Pid

# Episode dates are handy for calculating the number of patients in the hospital at a certain moment
# Pid date are handy for calculating incidence of hospitalization. One patient doesn't count twice

# Hospital admission...

# Hospital admission - EPISODE...
# Toon voortgang
message("Opschonen data nice opnamedatum")


# Make last discharge date per Pid per Episode (separately for hospital and IC)
data_nice <- data_nice %>% 
  group_by(Pid, Episode, IsIc) %>% 
  # calculate last discharge date, separately for hospital and IC
  mutate(DischargeDate_Episode = max(DischargeDate_temp, na.rm = TRUE)) %>%
  group_by(Pid, Episode) %>%
  # separate hospital and IC discharge date per episode; discharge from IC cannot be later than from hospital
  mutate(DischargeDate_Episode = 
           if_else(IsIc == 1, min(DischargeDate_Episode), DischargeDate_Episode),
         DischargeDate_Episode_Ic = 
           if_else(any(IsIc == 1), min(DischargeDate_Episode), as.Date(NA)),
         DischargeDate_Episode = max(DischargeDate_Episode)) %>%
  # bring discharge date at ANALYSISDATE + 1 back to NA
  mutate(DischargeDate_Episode = if_else(DischargeDate_Episode == ANALYSISDATE + 1, as.Date(NA), DischargeDate_Episode),
         DischargeDate_Episode_Ic = if_else(DischargeDate_Episode_Ic == ANALYSISDATE + 1, as.Date(NA), DischargeDate_Episode_Ic)) 

# Make first admission date (do not use negative admissions to calcualate the first admission date)
data_nice <- data_nice %>% 
  # Remove negative results to find first positive admission date
  filter(Covid19Status %in% c("lab", "ct")) %>% 
  group_by(Pid, Episode) %>% 
  # calculate first admission date
  mutate(AdmissionDate_Episode = min(AdmissionDate, na.rm = TRUE)) %>% 
  ungroup %>% 
  select(Pid, Episode, AdmissionDate_Episode) %>% 
  distinct %>% 
  # merge this info back to the data
  right_join(data_nice, by = c("Pid", "Episode"))

# remove records of first admissions with negative resuls %>% 
data_nice <- data_nice %>% 
  group_by(Pid, Episode) %>% 
  filter( !(AdmissionDate <= AdmissionDate_Episode & !Covid19Status %in% c("lab", "ct")) )   %>% 
  ungroup

# Hospital admission - Pid

# Make first admission date and first report date per Pid
data_nice <- data_nice %>%
  group_by(Pid) %>%
  mutate(
    AdmissionDate_Pid    = min(AdmissionDate_Episode, na.rm = TRUE),
    CreationDateNice_Pid = min(CreationDateNice, na.rm = TRUE) %>% as.POSIXct(tz = "UTC"))

# IC admission...

# Find IC admission day and report day 
tmp <- data_nice %>% 
  select(Pid, Episode, AdmissionDate, CreationDateNice, IsIc) %>% 
  filter(IsIc == 1) %>%
  # IC admission day and discharge day per episode
  group_by(Pid, Episode)  %>% 
  mutate(
    AdmissionDate_Episode_Ic = min(AdmissionDate, na.rm = TRUE)) %>% 
  ungroup %>% 
  # IC admission day and report day  per Pid
  group_by(Pid) %>% 
  mutate(
    AdmissionDate_Pid_Ic    = min(AdmissionDate, na.rm = TRUE),
    CreationDateNice_Pid_Ic = min(CreationDateNice, na.rm = TRUE)) %>% 
  ungroup %>% 
  select(!c(AdmissionDate, IsIc, CreationDateNice)) %>%  
  distinct

# Join it back to the data
data_nice <- tmp  %>% 
  right_join(data_nice, by = c("Pid", "Episode")) %>%
  select(!DischargeDate_temp)

rm(tmp)
