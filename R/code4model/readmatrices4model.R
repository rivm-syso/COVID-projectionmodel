###################################
### collect all contactmatrices ###
###################################

# possibility to add alternative matrices
files2use <- c(
  "ContactmatricesD3praktijk_midpoint_24mrt2020.rds", # pre-control matrix `99`
  intellockdown = "Contactmatrices_intelligentlockdown_march_2020-11-26.rds", # during first lockdown
  batch1 = "Contactmatrices_batch1_11may_2020-11-26.rds", # starting 11 May 2020
  batch2 = "Contactmatrices_batch2_1june_2020-11-26.rds", # starting 02 June 2020
  batch3 = "Contactmatrices_batch3_summerholiday_2020-11-26.rds", # starting 06 July 2020
  sep2020 = "Contactmatrices_start-schoolyear-september_2020-11-26.rds", # starting 31 August 2020
  okt2020 = "Contactmatrices_restrictions28sep_2020-11-26.rds", # starting 30 September 2020
  okt2020holiday = "Contactmatrices_partiallockdown-october-holiday_2020-11-26.rds", # starting 15 October 2020
  partlockdown = "Contactmatrices_partiallockdown-october_2020-11-26.rds", # starting 26 October & 19 November & 6 January
  nov2020 = "Contactmatrices_2week-lockdown-november_2020-11-26.rds", # starting 5 November
  winterlockdown = "Contactmatrices_winter-lockdown_2020-12-16.rds", # starting 15 december, until (including) 24 Jan (except holiday)
  winterlockdownChristmas2020 = "Contactmatrices_winter-lockdown-christmas_2020-12-16.rds"  # between Mo 21 Dec and Su 3 Jan
)

ContactMatrices <- list(readRDS(paste0("data/contacts/", files2use[1]))[[1]])
names(ContactMatrices)[1] <- "precontrol_mean"
for(i in 2:length(files2use)) {
  ContactMatrices <- c(ContactMatrices, list(readRDS(paste0("data/contacts/", files2use[i]))))
  names(ContactMatrices)[i] <- names(files2use)[i]
}
for(i in 2:length(ContactMatrices)) {
  toadd <- list(Reduce(`+`, ContactMatrices[[i]])/200)
  names(toadd)[1] <- paste0(names(ContactMatrices)[i], "_mean")
  ContactMatrices <- c(ContactMatrices, toadd)
}

rm(files2use, i, toadd)
