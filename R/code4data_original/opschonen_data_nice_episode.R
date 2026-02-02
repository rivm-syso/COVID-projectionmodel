# Toon voortgang
message("Opschonen data nice episode")

# Admissions with 3 or less day difference in between are the same admission episode
data_nice <- data_nice %>% 
  group_by(Pid, IsIc) %>%
  arrange(AdmissionDate, DischargeDate) %>%
  mutate(
      # impute discharge dates if missing, so that they can be used in defining episodes
    DischargeDate_temp = if_else(
      # first, if there is a later admission, use that as discharge date
      condition = is.na(DischargeDate), 
      true = lead(AdmissionDate), false = DischargeDate),
    DischargeDate_temp = if_else(
      # second, if there is a later discharge in an earlier record, use the admission date as discharge date (transfer on same day)
      condition = is.na(DischargeDate_temp) & lag(cummax(as.numeric(DischargeDate))) > as.numeric(AdmissionDate), 
      true = AdmissionDate, false = DischargeDate_temp),
    DischargeDate_temp = if_else(
      # third, assume that the patient has not been discharged yet, by using the current day + 1
      condition = is.na(DischargeDate_temp), 
      true = ANALYSISDATE + 1, false = DischargeDate_temp)) %>%
  group_by(Pid) %>% 
  # Organize the records inside a Pid per AdmissionDate, DischargeDate, IsIc
  arrange(AdmissionDate, DischargeDate, IsIc) %>% 
  mutate(
    Episode = if_else(
      # admission greater than 3 days after other admission/discharge of this patient is seen as a next event
      # condition = AdmissionDate > min(AdmissionDate) + 3 | AdmissionDate > max(DischargeDate) + 3,
      condition = 
        AdmissionDate >= lag(AdmissionDate) + 3 & as.numeric(AdmissionDate) >= lag(cummax(as.numeric(DischargeDate_temp))) + 3,
      # every episode starts with 1
      # the records equal to zero (which are not greater than 3 days) follows the previous "1".
      # therefore, with cumulative sum, every 1 counts as a new/extra episode, 
      # while 0 just takes the last episode number
      # missing marks the start of every Pid. so it is always 1
      true = 1, false = 0, missing = 1),
    # Make a temporarily Episode_tmp, just to be able to aggregate data per episode
    # But the real episode sequential number is recalculated in the next step
    Episode_tmp = cumsum(Episode)) %>% 
  ungroup 

# Remove episode per Pid with only non-positive results
data_nice <- data_nice %>% 
  group_by(Pid, Episode_tmp) %>% 
  filter(any(Covid19Status %in% c("lab", "ct"))) %>%  
  ungroup 

# Recalculate cumsum (now only with Pid episode with at least one positive result)
data_nice <- data_nice %>% 
  group_by(Pid) %>% 
  mutate(
    Episode = cumsum(Episode)) %>% 
  ungroup %>% 
  select(!Episode_tmp)

