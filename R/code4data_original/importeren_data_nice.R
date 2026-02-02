#
# importing the NICE data
#

message("Importing data nice")
files <- list.files("data_original/NICE", full.names = TRUE)
files <- files[grepl(files, pattern = "rds")]
files <- files[grepl(files, pattern = "ruw")]
filedates <- as.Date(
  substr(files, regexpr("[[:digit:]]{8}", files), regexpr("[[:digit:]]{8}", files) + 7),
  "%Y%m%d")
if(exists("ANALYSISDATE")) {
  file.name <- files[filedates == ANALYSISDATE + 1]
} else {
  file.name <- files[filedates %>% which.max]
  ANALYSISDATE <- max(filedates) - 1
}
data_nice_org <- readRDS(file.name)

rm(files, filedates, file.name)


