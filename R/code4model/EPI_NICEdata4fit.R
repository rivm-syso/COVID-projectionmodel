### cleaning up to this point from raw NICE synthetic data is not included, as 
### creating reliable synthetic data for that part to work and create a 
### dataset similar to the original is too complicated.
### The code is provided in the folders 'code4model_original' and 'code4data_original'
NICEIDdataprocessed <- readRDS("data/NICE/NICEIDdata_20210106_synthetic.rds")
data_nice <- readRDS("data/NICE/data_nicedelay_20210106_synthetic.rds")


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



