#
# Opschonen data nice
#

# Toon voortgang
message("Opschonen data nice algemeen")

# Elke patient heeft meerdere records
# ZKH opname
# IC opname
# Verplatsing van afdeling
# Verplatsing naar ander zikenhuis


# Make a Pid of every record from Pid 13428 using PidSequenceNumber
# Pid 13428 contains multiple patients
# This Pid number is used to assing patients with an unkown BSN
data_nice <- data_nice_org %>% 
  mutate(Pid = if_else(condition = Pid == 13428, true = -PidSequenceNumber, false = Pid))


# Tranform POSIXct data in Date
data_nice <- data_nice %>% 
  mutate(
    AdmissionDate = as.Date(AdmissionDate),
    DischargeDate = as.Date(DischargeDate))

