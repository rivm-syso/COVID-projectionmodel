datum_vandaag <- lubridate::today()
source("R/code4data_original/importeren_data_nice.R")
source("R/code4data_original/opschonen_data_nice.R")

# Toon voortgang
message("Transformeer data nice voor model")


# choice: admission: count everyone once, add durations of separate admissions (HOSP), 
# but treat IC duration as first admission-to-last discharge (close look at data shows that in most cases admissions are missing from the data)
NICEIDdataprocessed <- data_nice %>%
  arrange(Pid, AdmissionDate, DischargeDate, IsIc) %>%
  mutate(DischargedTo = case_when(
           IsIc == 0 ~ "hospitalrecord",
           is.na(DischargedTo) ~ "censored",
           DischargedTo %in% as.character(c(1, 3, 4, 6)) ~ "HOSP",
           DischargedTo == "7" ~ "death",
           DischargedTo == "8" ~ "home",
           DischargedTo == "9" ~ "unknown",
           DischargedTo == "10" ~ "buitenland",
           DischargedTo %in% as.character(c(2, 5)) ~ "otherICU",
           TRUE ~ "hospitalrecord"
         )) %>%
  group_by(Pid) %>%
  mutate(
    admissionIC = any(IsIc == 1)
  ) %>%
  group_by(Pid) %>%
  # opnameZH and laatsteZH are meant to (1) get regional data and (2) to remove records for badly registering hospitals
  mutate(
    # first hospital
    opnameZH = Hospital[1],
    # last hospital
    laatsteZH = tail(Hospital, 1)
  ) %>% 
  # per admissionseq
  # join Episodes with IC admission
  group_by(Pid, Episode) %>%
  mutate(
    everIC = any(IsIc == 1)
  ) %>%
  group_by(Pid) %>%
  mutate(
    Episode = if_else(everIC, min(c(100, Episode[IsIc == 1])), Episode)
  ) %>%
  group_by(Pid, Episode) %>%
  mutate(
    IC2D = any(DischargedTo == "death"),
    IC2C = any(DischargedTo == "home"),
    IC2G = any(DischargedTo == "buitenland"),
    HOSP2D = all(!is.na(DischargeDate_Episode)) & any(DiedInHospital == 1, na.rm = T) & !IC2D,
    HOSP2C = all(!is.na(DischargeDate_Episode)) & !IC2D & !IC2C & !HOSP2D,
    onIC = any(is.na(DischargeDate_Episode_Ic) & IsIc),
    inHOSP = !IC2D & !IC2C & !HOSP2D & !HOSP2C & !onIC,
    inHOSPalt = any(is.na(DischargeDate_Episode) & !onIC),
    opname = min(AdmissionDate_Episode),
    opnameIC = if_else(everIC, 
                       min(c(ANALYSISDATE, AdmissionDate_Episode_Ic), na.rm = T), 
                       as.Date(NA)),
    einde = if_else(!inHOSP,
                    max(c(as.Date("2020-01-01"), DischargeDate_Episode), na.rm = T), 
                    as.Date(NA)),
    eindeIC = if_else(everIC & !onIC, 
                      max(c(as.Date("2020-01-01"), DischargeDate_Episode_Ic), na.rm = T), 
                      as.Date(NA))
  ) %>%
  ungroup() %>%
  select(Pid, Episode, Age, opnameZH, laatsteZH, everIC, 
         IC2D, IC2C, IC2G, HOSP2D, HOSP2C, onIC, inHOSP,
         opname, opnameIC, einde, eindeIC) %>%
  distinct() 


### Add ROAZ regions
# read data
ROAZdata <- read_csv("data_original/population/ROAZregios.csv", col_types = cols(
  PC4_CODE = col_character(),
  Categorie_5 = col_character(),
  Categorie_model = col_character(),
  Gemeente = col_character(),
  Regio_ROAZ = col_character(),
  Regio_Ziekenhuisregio = col_character(),
  Ziekenhuis_loc_nr = col_character(),
  Ziekenhuis_org_nr = col_character(),
  Bevolking_2019 = col_double()
))     # PC4 + hospital codes + ROAZ regions
ICUzhcodes <- read_csv("data_original/population/NICE_zhcodes.csv", col_types = cols(
  name = col_character(),
  Ziekenhuis_loc_nr = col_double(),
  Ziekenhuis_org_nr = col_double()
))   # hospital codes + hospital names as in NICE data

# match hospital to region
ZHROAZmatched <- ROAZdata %>% 
  group_by(Regio_ROAZ, Ziekenhuis_loc_nr) %>%
  # population in shared ROAZ-hospital region (some hospitals have more than one ROAZ)
  summarise(bevgrootte = sum(Bevolking_2019, na.rm = T)) %>%
  group_by(Ziekenhuis_loc_nr) %>%
  # attribute hospital to ROAZ with largest shared population
  mutate(ZH_ROAZ = Regio_ROAZ[bevgrootte == max(bevgrootte)]) %>%
  select(ZH_ROAZ, Ziekenhuis_loc_nr) %>%
  distinct()

# match ICU patients to region
NICEIDdataprocessed <- NICEIDdataprocessed %>%
  mutate(name = opnameZH) %>%
  # give hospitals their codes
  left_join(ICUzhcodes, by = "name") %>%
  mutate(Ziekenhuis_loc_nr = as.character(Ziekenhuis_loc_nr),
         Ziekenhuis_org_nr = as.character(Ziekenhuis_org_nr)) %>%
  # match codes to ROAZ
  left_join(ZHROAZmatched, 
            by = c("Ziekenhuis_loc_nr")) 

# Categorise hospitals for reporting policy, separately by wave
BlackListedOpname <-
  NICEIDdataprocessed %>%
  select(opnameZH, everIC, opname) %>%
  mutate(maand = match(months(opname), month.name),
         wavenr = if_else(maand <= 6, 1, 2)) %>%
  group_by(wavenr, opnameZH) %>%
  summarise(pctIC = sum(everIC)/n(),
            Underreported = pctIC > 0.6) %>%
  group_by(opnameZH) %>%
  filter(any(Underreported))


########################
### incidence curves ###
########################

inc_nice_hosp <- function(ages = FALSE, ROAZ = "all", excludeIC = FALSE) {
  if(is.numeric(ROAZ)) {
    ROAZ <- sort(unique(NICEIDdataprocessed$ZH_ROAZ))[ROAZ]
  }
  
  HOSPdata <- NICEIDdataprocessed %>%
    filter(!everIC | !excludeIC) %>%
    # select(Pid, Age, opname, ZH_ROAZ) %>%
    select(Pid, Age, opname) %>%
    group_by(Pid) %>%
    mutate(opname = min(opname)) %>%
    ungroup() %>%
    distinct() %>%
    filter(opname <= ANALYSISDATE) %>%
    filter(!is.na(Age)) %>%
    mutate(
      ageclass = floor(Age / 10),
      ageclass = if_else(ageclass > 8, 8, ageclass)) %>%
    mutate(HOSPdag = as.numeric(opname - as.Date("2020-02-12"))) 
  
  if(ROAZ != "all") {
    HOSPdata <- HOSPdata %>%
      filter(ZH_ROAZ == ROAZ)
  }
  
  if(ages) {
    toreturn <- lapply(0:8, function(x) tabulate(HOSPdata %>% 
                                                   filter(ageclass == x) %>% 
                                                   pull(HOSPdag), nbins = as.numeric(ANALYSISDATE - as.Date("2020-02-12"))))
    names(toreturn) <- c("[0,10)", "[10,20)", "[20,30)", "[30,40)",
                         "[40,50)", "[50,60)", "[60,70)", "[70,80)", "[80,Inf]")
    return(toreturn)
  } else {
    return(HOSPdata %>%
             pull(HOSPdag) %>% 
             tabulate(nbins = as.numeric(ANALYSISDATE - as.Date("2020-02-12"))))
  }
}


inc_nice_ic <- function(ages = FALSE, ROAZ = "all", includeblacklisted = TRUE) {
  if(is.numeric(ROAZ)) {
    ROAZ <- sort(unique(NICEIDdataprocessed$ZH_ROAZ))[ROAZ]
  }
  
  if(includeblacklisted) {
    ICdata <- NICEIDdataprocessed
  } else {
    ICdata <- NICEIDdataprocessed %>%
      mutate(
        toremove = opname < as.Date("2020-07-01") & 
          opnameZH %in% (BlackListedOpname %>% filter(wavenr == 1 & Underreported) %>% pull(opnameZH)),
        toremove = toremove | 
          (opname >= as.Date("2020-07-01") &
             opnameZH %in% (BlackListedOpname %>% filter(wavenr == 2 & Underreported)  %>% pull(opnameZH)))) %>%
      filter(!toremove) %>%
      select(-toremove)
  }
  
  ICdata <- ICdata %>%
    # select(Age, opnameIC, ZH_ROAZ) %>%
    select(Age, opnameIC) %>%
    filter(opnameIC <= ANALYSISDATE) %>%
    filter(!is.na(Age)) %>%
    mutate(
      ageclass = floor(Age / 10),
      ageclass = if_else(ageclass > 8, 8, ageclass)) %>%
    mutate(ICdag = as.numeric(opnameIC - as.Date("2020-02-12")))
  
  if(ROAZ != "all") {
    ICdata <- ICdata %>%
      filter(ZH_ROAZ == ROAZ)
  }
  
  if(ages) {
    toreturn <- lapply(0:8, function(x) tabulate(ICdata %>% 
                                                   filter(ageclass == x) %>% 
                                                   pull(ICdag), nbins = as.numeric(ANALYSISDATE - as.Date("2020-02-12"))))
    names(toreturn) <- c("[0,10)", "[10,20)", "[20,30)", "[30,40)",
                         "[40,50)", "[50,60)", "[60,70)", "[70,80)", "[80,Inf]")
    return(toreturn)
  } else {
    return(ICdata %>%
             pull(ICdag) %>% 
             tabulate(nbins = as.numeric(ANALYSISDATE - as.Date("2020-02-12"))))
  }
}

#########################
### prevalence curves ###
#########################

prev_nice_hosp <- function(ages = FALSE, ROAZ = "all", maxstay = Inf) {
  if(is.numeric(ROAZ)) {
    ROAZ <- sort(unique(NICEIDdataprocessed$ZH_ROAZ))[ROAZ]
  }
  
  mutations <- NICEIDdataprocessed %>%
    # select(opname, einde, Age, ZH_ROAZ) %>%
    filter(opname <= ANALYSISDATE & (is.na(einde) | einde <= ANALYSISDATE)) %>%
    mutate(
      ageclass = floor(Age / 10),
      ageclass = if_else(ageclass > 8, 8, ageclass),
      einde = if_else(is.na(einde), ANALYSISDATE + 1, einde)) %>%
    mutate(hospdag = as.numeric(opname - as.Date("2020-02-12")),
           einddag = 1 + as.numeric(einde - as.Date("2020-02-12")))
  
  if(ROAZ != "all") {
    mutations <- mutations %>%
      filter(ZH_ROAZ == ROAZ)
  }
  
  if(ages) {
    opnames <- lapply(0:8, function(x) tabulate(mutations %>% 
                                                  filter(ageclass == x) %>% 
                                                  pull(hospdag), nbins = as.numeric(ANALYSISDATE - as.Date("2020-02-12"))))
    ontslagen <- lapply(0:8, function(x) tabulate(mutations %>% 
                                                    filter(ageclass == x) %>% 
                                                    pull(einddag), nbins = as.numeric(ANALYSISDATE - as.Date("2020-02-12"))))
    toreturn <- lapply(1:9, function(x) cumsum(opnames[[x]]) - cumsum(ontslagen[[x]]))
    names(toreturn) <- c("[0,10)", "[10,20)", "[20,30)", "[30,40)",
                         "[40,50)", "[50,60)", "[60,70)", "[70,80)", "[80,Inf]")
    return(toreturn)
  } else {
    opnames <- tabulate(mutations$hospdag, nbins = as.numeric(ANALYSISDATE - as.Date("2020-02-12")))
    ontslagen <- tabulate(mutations$einddag, nbins = as.numeric(ANALYSISDATE - as.Date("2020-02-12")))
    return(cumsum(opnames) - cumsum(ontslagen))
  }
}


prev_nice_ic <- function(ages = FALSE, ROAZ = "all") {
  if(is.numeric(ROAZ)) {
    ROAZ <- sort(unique(NICEIDdataprocessed$ZH_ROAZ))[ROAZ]
  }
  
  mutations <- NICEIDdataprocessed %>%
    # select(opnameIC, eindeIC, Age, ZH_ROAZ) %>%
    filter(opnameIC <= ANALYSISDATE & (is.na(eindeIC) | eindeIC <= ANALYSISDATE)) %>%
    mutate(
      ageclass = floor(Age / 10),
      ageclass = if_else(ageclass > 8, 8, ageclass),
      eindeIC = if_else(is.na(eindeIC), ANALYSISDATE + 1, eindeIC)) %>%
    mutate(ICdag = as.numeric(opnameIC - as.Date("2020-02-12")),
           einddag = 1 + as.numeric(eindeIC - as.Date("2020-02-12")))
  
  if(ROAZ != "all") {
    mutations <- mutations %>%
      filter(ZH_ROAZ == ROAZ)
  }
  
  if(ages) {
    opnames <- lapply(0:8, function(x) tabulate(mutations %>% 
                                                  filter(ageclass == x) %>% 
                                                  pull(ICdag), nbins = as.numeric(ANALYSISDATE - as.Date("2020-02-12"))))
    ontslagen <- lapply(0:8, function(x) tabulate(mutations %>% 
                                                    filter(ageclass == x) %>% 
                                                    pull(einddag), nbins = as.numeric(ANALYSISDATE - as.Date("2020-02-12"))))
    toreturn <- lapply(1:9, function(x) cumsum(opnames[[x]]) - cumsum(ontslagen[[x]]))
    names(toreturn) <- c("[0,10)", "[10,20)", "[20,30)", "[30,40)",
                         "[40,50)", "[50,60)", "[60,70)", "[70,80)", "[80,Inf]")
    return(toreturn)
  } else {
    opnames <- tabulate(mutations$ICdag, nbins = as.numeric(ANALYSISDATE - as.Date("2020-02-12")))
    ontslagen <- tabulate(mutations$einddag, nbins = as.numeric(ANALYSISDATE - as.Date("2020-02-12")))
    return(cumsum(opnames) - cumsum(ontslagen))
  }
}



################################
### Remove temporary objects ###
################################
rm(ROAZdata, ICUzhcodes, ZHROAZmatched
)

