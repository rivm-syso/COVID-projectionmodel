##########################################################
### Functions for population size and age distribution ###
##########################################################


# read data
ROAZdata <- read_csv("data/population/ROAZregios_synthetic.csv", col_types = cols(
  Categorie_5 = col_character(),
  Regio_ROAZ = col_character(),
  Bevolking_2019 = col_double()
))     # PC4 + hospital codes + ROAZ regions

agedist <- function(ROAZ = "all", ngroups = 9, ...) {
  if(is.numeric(ROAZ)) {
    ROAZ <- sort(unique(ROAZdata$Regio_ROAZ))[-9][ROAZ]
  }
  if(ROAZ == "all") {
      ROAZdata %>%
        group_by(Categorie_5) %>%
        summarise(bevgrootte = sum(Bevolking_2019, na.rm = T)) %>%
        separate(Categorie_5, into = c("lower", "upper")) %>%
        mutate(lower = as.numeric(lower),
               lower = if_else(is.na(lower), 95, lower),
               lower = floor(lower/10),
               lower = if_else(lower == 9, ngroups - 1, lower)) %>%
        group_by(lower) %>%
        summarise(bevgrootte = sum(bevgrootte)) %>%
        mutate(bevpct = bevgrootte/sum(bevgrootte)) %>% pull(bevpct)
    } else {
    ROAZdata %>%
      filter(Regio_ROAZ == ROAZ) %>%
      group_by(Categorie_5) %>%
      summarise(bevgrootte = sum(Bevolking_2019, na.rm = T)) %>%
      separate(Categorie_5, into = c("lower", "upper")) %>%
      mutate(lower = as.numeric(lower),
             lower = if_else(is.na(lower), 95, lower),
             lower = floor(lower/10),
             lower = if_else(lower == 9, ngroups - 1, lower)) %>%
      group_by(lower) %>%
      summarise(bevgrootte = sum(bevgrootte)) %>%
      mutate(bevpct = bevgrootte/sum(bevgrootte)) %>% pull(bevpct)
  }
} 
popsize <- function(ROAZ = "all", ...) {
  if(is.numeric(ROAZ)) {
    ROAZ <- sort(unique(ROAZdata$Regio_ROAZ))[-9][ROAZ]
  }
  if(ROAZ == "all") {
    return(ROAZdata %>%
             summarise(bevgrootte = sum(Bevolking_2019, na.rm = T)) %>% pull(bevgrootte))
  } else {
    return(ROAZdata %>%
             filter(Regio_ROAZ == ROAZ) %>%
             summarise(bevgrootte = sum(Bevolking_2019, na.rm = T))  %>% pull(bevgrootte))
  }
}

Populationresults <- list(
  popsizeregio = sapply(sort(unique(ROAZdata$Regio_ROAZ))[-9], popsize),
  popsizeNL = popsize(),
  agedistregio = sapply(sort(unique(ROAZdata$Regio_ROAZ))[-9], agedist),
  agedistNL = agedist(),
  agedistNL10 = agedist(ngroups = 10)
)

rm(popsize, agedist, ROAZdata)


popsize <- function(ROAZ = "all", ...) {
  if(ROAZ == "all") {
    return(Populationresults$popsizeNL)
  } else {
    return(Populationresults$popsizeregio[ROAZ])
  }
}

agedist <- function(ROAZ = "all", ngroups = 9, ...) {
  if(ROAZ == "all") {
    if(ngroups == 9) {
      return(Populationresults$agedistNL)
    } else {
      return(Populationresults$agedistNL10)
    }
  } else {
    return(Populationresults$agedistregio[, ROAZ])
  }
}
