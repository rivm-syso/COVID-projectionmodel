#
# Opschonen data nice
#

# Toon voortgang
message("Opschonen data nice filter")

# Elke patient heeft meerdere records
# ZKH opname
# IC opname
# Verplatsing van afdeling
# Verplatsing naar ander zikenhuis

# Filter records of today
# Data should contain only records until yesterday 23:59
data_nice <- data_nice %>% 
  filter(!as.Date(CreationDateNice) >= ANALYSISDATE + 1)

# # Make a Pid of every record from Pid 13428 using PidSequenceNumber
# # Pid 13428 contains multiple patients
# # This Pid number is used to assing patients with an unkown BSN
# data_nice_new <- data_nice_new %>% 
#   mutate(Pid = if_else(condition = Pid == 13428, true = -PidSequenceNumber, false = Pid))

# Remove records with only negative results/heropname per Pid 
data_nice <- data_nice %>% 
  group_by(Pid) %>% 
  filter(any(Covid19Status %in% c("lab", "ct"))) %>% 
  ungroup()

# Clean Age
# Some Pid have different ages per record
# Choose the age in the first notification created
# If this age is equal to zero, than choose the age of the next notification
data_nice <- data_nice %>% 
  group_by(Pid) %>% 
  # Arrange per CreationDateNice to later on choose the Age in the first notification
  arrange(ModifiedDateNice) %>% 
  mutate(
    # If Age is zero than exchange it for the next notification date
    # Therefore, exchange  first all Age equal zero (if other Age is avilable) to NA
    Age = if_else(condition = Age == 0 & any(Age > 0), true = NA_integer_, false = Age),
    # Choose the Age in the first non-NA notification
    Age = first(Age %>% na.omit)) %>% 
  ungroup

# Remove duplicates...

# # Tranform POSIXct data in Date
# data_nice_new <- data_nice_new %>% 
#   mutate(
#     AdmissionDate = as.Date(AdmissionDate),
#     DischargeDate = as.Date(DischargeDate))

# Remove duplicate records in Pid, AdmissionDate, IsIc, DischargeDate and Covid19Status
data_nice <- data_nice %>% 
  distinct(Pid, AdmissionDate, IsIc, DischargeDate, Covid19Status, .keep_all = TRUE) 

# Remove duplicates based on conditions...

# Identify duplicates
data_nice <- data_nice %>% 
  # find duplicated admission record per Pid
  count(Pid, AdmissionDate, IsIc, Hospital, name = "Duplicates") %>% 
  # join this info back to the data
  right_join(data_nice, by = c("Pid", "AdmissionDate", "IsIc", "Hospital"))

# Among the duplicates choose the records to be kept

# Later bind it back to the data_nice
tmp <- data_nice %>% 
  # Among duplicated admission record per Pid, AdmissionDate, IsIc...
  filter(Duplicates > 1) %>% 
  # make choices based on unique Pid, AdmissionDate, IsIc records
  group_by(Pid, AdmissionDate, IsIc, Hospital) %>% 
  # Keep the records with a discharge datum
  # Or the records with the latest discharge datum
  # Still keep the records with no discharge date at all
  filter(all(is.na(DischargeDate)) | DischargeDate == max(DischargeDate, na.rm = TRUE)   )  %>% 
  # Among the records with no discharge date at all or with duplicate dates, 
  # prefer the ones with covid status (if existent)
  mutate(Duplicated_DischargeDate = sum(n())) %>% 
  filter(!(
    # Among the records which are still duplicated 
    (Duplicated_DischargeDate > 1) & 
      # And with any covid status
      any(Covid19Status %in% c("lab", "ct")) & 
      # Remove the ones which are not covid
      !Covid19Status %in% c("lab", "ct"))) %>% 
  # Among the records which are still duplicated 
  # Keep the last updated record (PidSequenceNumber)
  filter(! (Duplicated_DischargeDate > 1 & !PidSequenceNumber == max(PidSequenceNumber, na.rm = TRUE) )) %>% 
  ungroup

# Rebind duplicates in tmp to non duplicates in data_nice
data_nice <- data_nice %>% 
  filter(Duplicates <= 1) %>% 
  bind_rows(tmp) %>% 
  select(!c(Duplicates, Duplicated_DischargeDate))

rm(tmp)


