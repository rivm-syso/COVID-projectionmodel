################################################################
### Functions to estimate probabilities and delays from NICE ###
################################################################

#####################
### Read the data ###
#####################

source("R/code4model/EPI_NICEdata4fit.R")

message("analyseer nice data: kansen en delays")

################
### ANALYSES ###
################
### further clean data ###
pddata <- NICEIDdataprocessed %>%
  
  # remove hospitals with underreporting in non-ICU wards
  filter(!(opnameZH %in% BlackListedOpname$opnameZH)) %>%
  
  # numeric dates (admission, admissionIC, end, endIC)
  mutate(
    opname = as.numeric(opname - as.Date("2020-02-12")),
    opnameIC = as.numeric(opnameIC - as.Date("2020-02-12")),
    einde = as.numeric(einde - as.Date("2020-02-12")),
    einde = if_else(inHOSP | onIC, as.numeric(ANALYSISDATE - as.Date("2020-02-12")), einde),
    eindeIC = as.numeric(eindeIC - as.Date("2020-02-12")),
    eindeIC = if_else(onIC, as.numeric(ANALYSISDATE - as.Date("2020-02-12")), eindeIC)
    ) %>%
  
  # turn times to durations, summing Episodes per patient, redefine state of each patient (IC2G becomes IC2C)
  group_by(Pid) %>%
  mutate(
    everIC = any(everIC),
    opnameIC = max(-100, opnameIC, na.rm = T),
    eindeIC = max(-100, eindeIC, na.rm = T),
    durnoIC = sum(einde - opname + 1) - 1,
    durpreIC = if_else(everIC, max(opnameIC - opname), NA_real_),
    duronIC = if_else(everIC, max(eindeIC - opnameIC), NA_real_),
    durpostIC = if_else(everIC, min((einde - eindeIC)[einde - eindeIC >= 0]) +
                          sum((einde - opname + 1)[opname > opnameIC]), NA_real_),
    HOSP2D = any(HOSP2D),
    inHOSP = any(inHOSP),
    IC2D = any(IC2D),
    onIC = any(onIC),
    HOSP2C = !(HOSP2D | inHOSP | IC2D | onIC) & (!everIC | eindeIC < max(einde)),
    IC2C = !(HOSP2D | inHOSP | IC2D | onIC) & everIC & (eindeIC == max(einde))
  ) %>% 
  
  # classify durations by outcome, define admission of first Episode as 'the' admission
  mutate(ti = min(opname),
         durA2IC = durpreIC,
         durA2D = if_else(!everIC & HOSP2D, durnoIC, NA_real_),
         durA2C = if_else(!everIC & HOSP2C, durnoIC, NA_real_),
         durA2U = if_else(!everIC & inHOSP, durnoIC, NA_real_),
         durIC2H = if_else(everIC & !(IC2D | IC2C | onIC), duronIC, NA_real_),
         durIC2D = if_else(IC2D, duronIC, NA_real_),
         durIC2C = if_else(IC2C, duronIC, NA_real_),
         durIC2U = if_else(onIC, duronIC, NA_real_),
         durH2D = if_else(everIC & HOSP2D, durpostIC, NA_real_),
         durH2C = if_else(everIC & !(HOSP2D | inHOSP | IC2D | onIC | IC2C), durpostIC, NA_real_),
         durH2U = if_else(everIC & inHOSP, durpostIC, NA_real_)) %>%
  ungroup() %>%
  select(-Episode) %>%
  
  # censor admission time and age for fitting
  mutate(ti = pmax(ti,0),
         age = pmin(90, Age)) %>%
  select(Pid, ti, age, HOSP2D, HOSP2C, inHOSP, IC2D, IC2C, onIC, everIC, 
         durA2IC, durA2D, durA2C, durA2U, durIC2H, durIC2D, durIC2C, durIC2U, durH2D, durH2C, durH2U) %>%
  filter(!is.na(ti) & !is.na(age)) %>%
  distinct()
pddata4fit <- pddata %>%
  mutate(age = floor(age/10),
         age = pmin(8, age))
pddata4fitdurs <- pddata4fit %>%
  
  # remove earliest cases, often earlier admitted for non-COVID reasons, and therefore having unusual lengths of stay
  filter(ti > 20) %>%
  
  # distinguish censored from uncensored observations
  mutate(
    intervalA2IC = if_else(!is.na(durA2IC), pmin(30, durA2IC), durA2IC),
    intervalA2D_u = pmin(durA2C, durA2D, na.rm = T),
    intervalA2D_c = if_else(!is.na(durA2U), pmin(30, durA2U), durA2U),
    intervalIC2D_u = pmin(durIC2D, durIC2C, durIC2H, na.rm = T),
    intervalIC2D_c = if_else(!is.na(durIC2U), pmin(100, durIC2U), durIC2U),
    intervalH2D_u = pmin(durH2D, durH2C, na.rm = T),
    intervalH2D_c = if_else(!is.na(durH2U), pmin(30, durH2U), durH2U)
  )  %>% 
  filter(
    (is.na(intervalA2IC) | intervalA2IC >= 0) & 
      (is.na(intervalA2D_u) | intervalA2D_u >= 0) & 
      (is.na(intervalA2D_c) | intervalA2D_c >= 0) & 
      (is.na(intervalIC2D_u) | intervalIC2D_u >= 0) & 
      (is.na(intervalIC2D_c) | intervalIC2D_c >= 0) & 
      (is.na(intervalH2D_u) | intervalH2D_u >= 0) & 
      (is.na(intervalH2D_c) | intervalH2D_c >= 0)
  )

################################
### Define all probabilities ###
################################
### The four basic probabilities
# probability to die untreated (age, first wave)
prob_death_untreated <- function(a, ad0, bd0) {
  exp(ad0 + bd0 * a) / (1 + exp(ad0 + bd0 * a)) 
}
# probability to be hospitalised among those that will die (changes in first wave)
prob_admission_inD <- function(a, t, pAinDmin, rate1, tmid1) {
  (1 / (1 + exp(rate1 * (t - tmid1))))  +
    pAinDmin[a] * exp(rate1 * (t - tmid1)) / (1 + exp(rate1 * (t - tmid1)))
}
# probability to die among admitted
prob_death_inA <- function(a, t, ad0, bd0, pAinDmin, rate1, tmid1) {
  prob_admission_inD(a, t, pAinDmin, rate1, tmid1) *
    prob_death_untreated(a, ad0, bd0) /
    (1 - prob_death_untreated(a, ad0, bd0) *
       (1 - prob_admission_inD(a, t, pAinDmin, rate1, tmid1)))
}
# probability to go to ICU (one age group)
prob_ic_age <- function(t, L1, L2, rate1, tmid1) {
  L1 +
    (L2 - L1) * exp(rate1 * (t - tmid1)) / (1 + exp(rate1 * (t - tmid1)))  
}
# mean delay (A2IC, A2D, IC2D, H2D)
mean_DT <- function(t, mu1, mu2, mu3, mu4, rate1, rate2, rate3, tmid1, tmid2, tmid3) {
  mu1 +
    (mu2 - mu1) * exp(rate1 * (t - tmid1)) / (1 + exp(rate1 * (t - tmid1))) +
    (mu3 - mu2) * exp(rate2 * (t - tmid2)) / (1 + exp(rate2 * (t - tmid2))) +
    (mu4 - mu3) * exp(rate3 * (t - tmid3)) / (1 + exp(rate3 * (t - tmid3)))
}

##################################################################
### likelihood 1: probability of admission, probability to die ###
##################################################################
minloglik1 <- function(parms) {
  - sum(log(prob_death_inA(1 + pddata4fit$age[pddata4fit$IC2D | pddata4fit$HOSP2D], 
                           pddata4fit$ti[pddata4fit$IC2D | pddata4fit$HOSP2D],
                           ad0 = parms[1], bd0 = parms[2],
                           pAinDmin = c(rep(1, 4), exp(parms[3:7])/ (1 + exp(parms[3:7]))),
                           rate1 = exp(parms[8]), tmid1 = 10*parms[9]))) -
    sum(log(1 - prob_death_inA(1 + pddata4fit$age[pddata4fit$IC2C | pddata4fit$HOSP2C], 
                           pddata4fit$ti[pddata4fit$IC2C | pddata4fit$HOSP2C],
                           ad0 = parms[1], bd0 = parms[2],
                           pAinDmin = c(rep(1, 4), exp(parms[3:7])/ (1 + exp(parms[3:7]))),
                           rate1 = exp(parms[8]), tmid1 = 10*parms[9])))
}
res1 <- optim(c(-7.1, 0.8, c(.1,-.5,-0.1,0.1,-.6), -1.8, 4.1),
              minloglik1, method = "CG")
bestfitparmsdeath <- c(as.list(res1$par[1:2]), 
                       list(c(rep(1, 4), exp(res1$par[3:7]) / (1 + exp(res1$par[3:7])))),
                       as.list(c(exp(res1$par[8]),10*res1$par[9])))

######################################################################
### likelihood 2: probability of going to the ICU, given admission ###
######################################################################
minloglik2 <- function(parms, age) {
  - sum(log(prob_ic_age(pddata4fit$ti[pddata4fit$everIC & pddata4fit$age == age],
                        exp(parms[1])/(1 + exp(parms[1])), exp(parms[2])/(1 + exp(parms[2])),
                        bestfitparmsdeath[[4]], bestfitparmsdeath[[5]]))) -
    sum(log(1 - prob_ic_age(pddata4fit$ti[!pddata4fit$everIC & pddata4fit$age == age],
                            exp(parms[1])/(1 + exp(parms[1])), exp(parms[2])/(1 + exp(parms[2])),
                            bestfitparmsdeath[[4]], bestfitparmsdeath[[5]]))) 
}
res2 <- sapply(0:8, function(x) optim(c(rep(-1,2)), minloglik2, age = x)$par)
bestfitparmsic <- exp(res2) / (1 + exp(res2))

#########################################################################
### likelihood 2a: probability of going back to HOSP, given admission ###
#########################################################################
minloglik2a <- function(parms, age) {
  - sum(log(prob_ic_age(pddata4fit$ti[pddata4fit$everIC & 
                                        (pddata4fit$HOSP2D | pddata4fit$HOSP2C | pddata4fit$inHOSP) & 
                                        pddata4fit$age %in% age],
                        exp(parms[1])/(1 + exp(parms[1])), exp(parms[2])/(1 + exp(parms[2])),
                        bestfitparmsdeath[[4]], bestfitparmsdeath[[5]]))) -
    sum(log(1 - prob_ic_age(pddata4fit$ti[pddata4fit$everIC & 
                                            (pddata4fit$IC2D | pddata4fit$IC2C) & pddata4fit$age %in% age],
                            exp(parms[1])/(1 + exp(parms[1])), exp(parms[2])/(1 + exp(parms[2])),
                            bestfitparmsdeath[[4]], bestfitparmsdeath[[5]]))) 
}
res2a <- 
  cbind(matrix(rep(optim(c(rep(-1,2)), minloglik2a, age = 0:4)$par, 5), nrow = 2),
        sapply(5:8, function(x) optim(c(rep(-1,2)), minloglik2a, age = x)$par))
bestfitparmshosp <- exp(res2a) / (1 + exp(res2a))


#######################################
### likelihood 3: durations of stay ###
#######################################
minlogLik3 <- function(pars, DTuncensored, tuncensored, DTcensored, tcensored) {
  # voltooide ligtijden: dnbinom(DT, ...) is de kans op de waargenomen ligtijd
  - sum(dnbinom(DTuncensored, exp(pars[11]),
                mu = mean_DT(tuncensored, exp(pars[1]), exp(pars[2]), exp(pars[3]), exp(pars[4]),
                             exp(pars[5]), exp(pars[6]), exp(pars[7]),
                             pars[8], pars[9], pars[10]), log = T)) -
  # onvoltooide ligtijden: pnbinom(DT - 1, ..., lower.tail = FALSE) is de kans op een ligtijd groter dan DT - 1
    sum(pnbinom(DTcensored - 1, exp(pars[11]),
                mu = mean_DT(tcensored, exp(pars[1]), exp(pars[2]), exp(pars[3]), exp(pars[4]),
                             exp(pars[5]), exp(pars[6]), exp(pars[7]),
                             pars[8], pars[9], pars[10]), lower.tail = FALSE, log.p = T))
}
res3 <- optim(c(log(c(2,2,1,2,.1,.1,.1)),20,40,160,1), minlogLik3, 
              DTuncensored = pddata4fitdurs$intervalA2IC[!is.na(pddata4fitdurs$intervalA2IC)],
              tuncensored = pddata4fitdurs$ti[!is.na(pddata4fitdurs$intervalA2IC)],
              DTcensored = c(),
              tcensored = c())$par
res3 <- cbind(res3, optim(c(log(c(10,10,10,10,.1,.1,.1)),20,110,180,1), minlogLik3, 
                          DTuncensored = pddata4fitdurs$intervalA2D_u[!is.na(pddata4fitdurs$intervalA2D_u)],
                          tuncensored = pddata4fitdurs$ti[!is.na(pddata4fitdurs$intervalA2D_u)],
                          DTcensored = pddata4fitdurs$intervalA2D_c[!is.na(pddata4fitdurs$intervalA2D_c)],
                          tcensored = pddata4fitdurs$ti[!is.na(pddata4fitdurs$intervalA2D_c)])$par)
res3 <- cbind(res3, optim(c(log(c(20,20,20,20,.1,.1,.1)),10,30,80,1), minlogLik3, 
                          DTuncensored = pddata4fitdurs$intervalIC2D_u[!is.na(pddata4fitdurs$intervalIC2D_u)],
                          tuncensored = pddata4fitdurs$ti[!is.na(pddata4fitdurs$intervalIC2D_u)],
                          DTcensored = pddata4fitdurs$intervalIC2D_c[!is.na(pddata4fitdurs$intervalIC2D_c)],
                          tcensored = pddata4fitdurs$ti[!is.na(pddata4fitdurs$intervalIC2D_c)])$par)
res3 <- cbind(res3, optim(c(log(c(10,10,10,10,.1,.1,.1)),20,110,180,1), minlogLik3, 
                          DTuncensored = pddata4fitdurs$intervalH2D_u[!is.na(pddata4fitdurs$intervalH2D_u)],
                          tuncensored = pddata4fitdurs$ti[!is.na(pddata4fitdurs$intervalH2D_u)],
                          DTcensored = pddata4fitdurs$intervalH2D_c[!is.na(pddata4fitdurs$intervalH2D_c)],
                          tcensored = pddata4fitdurs$ti[!is.na(pddata4fitdurs$intervalH2D_c)])$par)
bestfitparmsdur <- res3
bestfitparmsdur[c(1:7,11),] <- exp(res3[c(1:7,11),])


##############################
## possibility to check results with plots
##############################
# lowfunc <- function(nobs, xobs) {
#   if(xobs == 0) return(0)
#   uniroot(function(x) pbinom(xobs - 1, nobs, x) - 0.975, c(0,1))$root
# }
# uppfunc <- function(nobs, xobs) {
#   if(xobs == nobs) return(1)
#   uniroot(function(x) pbinom(xobs, nobs, x) - 0.025, c(0,1))$root
# }
# # function for estimated probability of death given admission
# probdeath_est <- function(time, age) {
#   prob_death_inA(age + 1, time, bestfitparmsdeath[[1]], bestfitparmsdeath[[2]],
#                  bestfitparmsdeath[[3]], bestfitparmsdeath[[4]],
#                  bestfitparmsdeath[[5]])
# }
# pddata4fit %>%
#   mutate(
#     tiperiode = ceiling(ti/20)
#   ) %>%
#   group_by(tiperiode, age) %>%
#     summarise(
#       nD = sum(HOSP2D | IC2D),
#       ntotal = sum(HOSP2D | IC2D | HOSP2C | IC2C)
#     ) %>%
#     ungroup() %>%
#     mutate(
#       probdie = nD/ntotal,
#       lolim = sapply(1:n(), function(x) lowfunc(ntotal[x], nD[x])),
#       uplim = sapply(1:n(), function(x) uppfunc(ntotal[x], nD[x]))
#     ) %>%
#   mutate(
#     estpdie = probdeath_est(20 * tiperiode - 10, age)
#   ) %>%
#   ggplot(aes(x = tiperiode, y = probdie)) +
#   geom_pointrange(aes(ymin = lolim, ymax = uplim)) +
#   geom_line(aes(y = estpdie)) +
#   facet_wrap(~age)
# 
# 
# # function for estimated probability of ic given admission
# prob_ic_est <- function(a, t) {
#   prob_ic_age(t, bestfitparmsic[1, a + 1], bestfitparmsic[2, a + 1],
#               bestfitparmsdeath[[4]], bestfitparmsdeath[[5]])
# }
# pddata4fit %>%
#   mutate(
#     tiperiode = ceiling(ti/20),
#     lft = pmax(6, age)
#   ) %>%
#   group_by(tiperiode, age) %>%
#   summarise(
#     nIC = sum(everIC),
#     ntotal = n()
#   ) %>%
#   ungroup() %>%
#   mutate(
#     probic = nIC/ntotal,
#     lolim = sapply(1:n(), function(x) lowfunc(ntotal[x], nIC[x])),
#     uplim = sapply(1:n(), function(x) uppfunc(ntotal[x], nIC[x]))
#   ) %>%
#   mutate(
#     estpic = prob_ic_est(age, 20 * tiperiode - 10)
#   ) %>%
#   ggplot(aes(x = tiperiode, y = probic)) +
#   geom_pointrange(aes(ymin = lolim, ymax = uplim)) +
#   geom_line(aes(y = estpic)) +
#   facet_wrap(~age, scales = "free_y")
# 
# pddata4fit %>%
#   filter(everIC) %>%
#   mutate(tiperiode = ceiling(ti/10),
#          ICdur = durA2IC,
#          sddur = sd(ICdur, na.rm = T)) %>%
#   # filter(ICdur <= 100) %>%
#   group_by(tiperiode) %>%
#   summarise(
#     meandur = mean(ICdur, na.rm = T),
#     sedur = mean(sddur/sqrt(n()))
#   ) %>%
#   ungroup() %>%
#   mutate(
#     lowdur = meandur - 1.96 * sedur,
#     uppdur = meandur + 1.96 * sedur
#   ) %>%
#   mutate(
#     estdur = mean_DT(10 * tiperiode - 5, bestfitparmsdur[1, 1], bestfitparmsdur[2, 1],
#                      bestfitparmsdur[3, 1], bestfitparmsdur[4, 1], bestfitparmsdur[5, 1],
#                      bestfitparmsdur[6, 1], bestfitparmsdur[7, 1], bestfitparmsdur[8, 1],
#                      bestfitparmsdur[9, 1], bestfitparmsdur[10, 1])) %>%
#   ggplot(aes(x = tiperiode, y = meandur)) +
#   geom_line(aes(y = estdur)) +
#   geom_pointrange(aes(ymin = lowdur, ymax = uppdur), size = 0.3) +
#   coord_cartesian(ylim = c(0,5)) +
#   theme_light() +
#   labs(title = "Gemiddelde ligduur op ZKH", y = "ZKH-ligduur")
# pddata4fit %>%
#   filter(!everIC) %>%
#   mutate(tiperiode = ceiling(ti/10),
#          ICdur = pmax(durA2D, durA2C, durA2U, na.rm = T),
#          sddur = sd(ICdur, na.rm = T)) %>%
#   # filter(ICdur <= 100) %>%
#   group_by(tiperiode) %>%
#   summarise(
#     meandur = mean(ICdur, na.rm = T),
#     sedur = mean(sddur/sqrt(n()))
#   ) %>%
#   ungroup() %>%
#   mutate(
#     lowdur = meandur - 1.96 * sedur,
#     uppdur = meandur + 1.96 * sedur
#   ) %>%
#   mutate(
#     estdur = mean_DT(10 * tiperiode - 5, bestfitparmsdur[1, 2], bestfitparmsdur[2, 2],
#                      bestfitparmsdur[3, 2], bestfitparmsdur[4, 2], bestfitparmsdur[5, 2],
#                      bestfitparmsdur[6, 2], bestfitparmsdur[7, 2], bestfitparmsdur[8, 2],
#                      bestfitparmsdur[9, 2], bestfitparmsdur[10, 2])) %>%
#   ggplot(aes(x = tiperiode, y = meandur)) +
#   geom_line(aes(y = estdur)) +
#   geom_pointrange(aes(ymin = lowdur, ymax = uppdur), size = 0.3) +
#   coord_cartesian(ylim = c(0,30)) +
#   theme_light() +
#   labs(title = "Gemiddelde ligduur op ZKH", y = "ZKH-ligduur")
# pddata4fit %>%
#   filter(everIC) %>%
#   mutate(tiperiode = ceiling(ti/10),
#          ICdur = pmin(durIC2H, durIC2D, durIC2C, durIC2U, na.rm = T),
#          sddur = sd(ICdur, na.rm = T)) %>%
#   filter(ICdur <= 100) %>%
#   group_by(tiperiode) %>%
#   summarise(
#     meandur = mean(ICdur, na.rm = T),
#     sedur = mean(sddur/sqrt(n()))
#   ) %>%
#   ungroup() %>%
#   mutate(
#     lowdur = meandur - 1.96 * sedur,
#     uppdur = meandur + 1.96 * sedur
#   ) %>%
#   mutate(
#     estdur = mean_DT(10 * tiperiode - 5, bestfitparmsdur[1,3], bestfitparmsdur[2, 3],
#                      bestfitparmsdur[3, 3], bestfitparmsdur[4, 3], bestfitparmsdur[5, 3],
#                      bestfitparmsdur[6, 3], bestfitparmsdur[7, 3], bestfitparmsdur[8, 3],
#                      bestfitparmsdur[9, 3], bestfitparmsdur[10, 3])) %>%
#   ggplot(aes(x = tiperiode, y = meandur)) +
#   geom_line(aes(y = estdur)) +
#   geom_pointrange(aes(ymin = lowdur, ymax = uppdur), size = 0.3) +
#   coord_cartesian(ylim = c(0,30)) +
#   theme_light() +
#   labs(title = "Gemiddelde ligduur op IC", y = "IC-ligduur")
# pddata4fit %>%
#   filter(!is.na(durH2D) | !is.na(durH2C) | !is.na(durH2U)) %>%
#   mutate(tiperiode = ceiling(ti/10),
#          ICdur = pmin(durH2D, durH2C, durH2U, na.rm = T),
#          sddur = sd(ICdur, na.rm = T)) %>%
#   # filter(ICdur <= 100) %>%
#   group_by(tiperiode) %>%
#   summarise(
#     meandur = mean(ICdur, na.rm = T),
#     sedur = mean(sddur/sqrt(n()))
#   ) %>%
#   ungroup() %>%
#   mutate(
#     lowdur = meandur - 1.96 * sedur,
#     uppdur = meandur + 1.96 * sedur
#   ) %>%
#   mutate(
#     estdur = mean_DT(10 * tiperiode - 5, bestfitparmsdur[1,4], bestfitparmsdur[2, 4],
#                      bestfitparmsdur[3, 4], bestfitparmsdur[4, 4], bestfitparmsdur[5, 4],
#                      bestfitparmsdur[6, 4], bestfitparmsdur[7, 4], bestfitparmsdur[8, 4],
#                      bestfitparmsdur[9, 4], bestfitparmsdur[10, 4])) %>%
#   ggplot(aes(x = tiperiode, y = meandur)) +
#   geom_line(aes(y = estdur)) +
#   geom_pointrange(aes(ymin = lowdur, ymax = uppdur), size = 0.3) +
#   coord_cartesian(ylim = c(0,30)) +
#   theme_light() +
#   labs(title = "Gemiddelde ligduur op IC", y = "IC-ligduur")


### check with latest results that required adaptation (data upto 04 nov 2020)
if(any(abs(res1$par - c(-7.1,0.82,0.02,-.56,-.18,.01,-.68,-1.9,4.0)) > 0.5)) {
  cat("Please check NICE results 1 carefully: \nExpected: ")
  cat(c(-7.1,0.82,0.02,-.56,-.18,.01,-.68,-1.9,4.0))
  cat("\nObserved: ", res1$par)
}
if(any(abs(bestfitparmsic[2,] - c(0,.06,.10,.11,.16,.20,.24,.20,.03)) > 0.05)) {
  cat("Please check NICE results 2 carefully: \nExpected: ")
  cat(c(0,.06,.10,.11,.16,.20,.24,.20,.03))
  cat("\nObserved: ", bestfitparmsic[2,])
}
if(any(abs(bestfitparmshosp[2,] - c(.87,.87,.87,.87,.87,.83,.72,.57,.55)) > 0.05)) {
  cat("Please check NICE results 2a carefully: \nExpected: ")
  cat(c(.87,.87,.87,.87,.87,.83,.72,.57,.55))
  cat("\nObserved: ", bestfitparmshosp[2,])
}
# aangepast data 6/11
if(any(abs(bestfitparmsdur[4,1:4] - c(2.3,7.9,15.4,9.9)) > 1.0)) {
  cat("Please check NICE results 3 carefully: \nExpected: ")
  cat(c(2.3,7.9,15.4,9.9))
  cat("\nObserved: ", bestfitparmsdur[4,1:4])
}


########################################################
### Functions to obtain all probabilities and delays ###
########################################################
### function to calculate the probabilities
get_allprobs <- function(tmax) {
  a <- matrix(rep(1:9, each = tmax), ncol = 9)
  t <- matrix(rep(1:tmax, 9), ncol = 9)
  
  pD <- prob_death_untreated(a, bestfitparmsdeath[[1]], bestfitparmsdeath[[2]])
  pAinD <- prob_admission_inD(a, t, bestfitparmsdeath[[3]], bestfitparmsdeath[[4]],
                              bestfitparmsdeath[[5]])
  pIC <- sapply(1:9, function(x) 
    prob_ic_age(t[, x],
                bestfitparmsic[1, x], bestfitparmsic[2, x], 
                bestfitparmsdeath[[4]], bestfitparmsdeath[[5]]))
  pHOSP <- sapply(1:9, function(x) 
    prob_ic_age(t[, x],
                bestfitparmshosp[1, x], bestfitparmshosp[2, x], 
                bestfitparmsdeath[[4]], bestfitparmsdeath[[5]]))

  probSe2A <- pD * pAinD + 1 - pD
  probSe2D <- pD * (1 - pAinD)
  
  probA2IC <- pIC
  probA2D <- pD * pAinD / probSe2A
  probA2C <- (1 - pD) / probSe2A
  probIC2H <- pHOSP

  probSe2A2IC <- probSe2A * probA2IC
  probSe2A2DC <- probSe2A * (1 - probA2IC)

  probSe2A2IC2H <- probSe2A2IC * probIC2H
  probSe2A2IC2DC <- probSe2A2IC * (1 - probIC2H)

  return(list(
    probSe2A = probSe2A,
    probSe2D = probSe2D,
    probA2IC = probA2IC,
    probA2D = probA2D,
    probA2C = probA2C,
    probIC2H = probIC2H,
    probSe2A2IC = probSe2A2IC,
    probSe2A2DC = probSe2A2DC,
    probSe2A2IC2H = probSe2A2IC2H,
    probSe2A2IC2DC = probSe2A2IC2DC
  ))
}

### function to calculate the delays
normalisedelays <- function(delays, maxdelay) {
  toreturn <- delays[1:(1 + maxdelay), ]
  toreturn <- toreturn / rep(colSums(toreturn), each = maxdelay + 1)
  return(toreturn)
}
get_alldelays <- function(tmax, maxdelay = 100) {
  delayA2IC <- sapply(mean_DT(1:tmax, bestfitparmsdur[1, 1], bestfitparmsdur[2, 1],
                              bestfitparmsdur[3, 1], bestfitparmsdur[4, 1], bestfitparmsdur[5, 1], 
                              bestfitparmsdur[6, 1], bestfitparmsdur[7, 1], bestfitparmsdur[8, 1], 
                              bestfitparmsdur[9, 1], bestfitparmsdur[10, 1]), 
                      function(x) dnbinom(0:(2 * maxdelay), bestfitparmsdur[11, 1], mu = x))
  delayA2D <- sapply(mean_DT(1:tmax, bestfitparmsdur[1, 2], bestfitparmsdur[2, 2],
                             bestfitparmsdur[3, 2], bestfitparmsdur[4, 2], bestfitparmsdur[5, 2], 
                             bestfitparmsdur[6, 2], bestfitparmsdur[7, 2], bestfitparmsdur[8, 2], 
                             bestfitparmsdur[9, 2], bestfitparmsdur[10, 2]), 
                     function(x) dnbinom(0:(2 * maxdelay), bestfitparmsdur[11, 2], mu = x))
  delayIC2D <- sapply(mean_DT(1:tmax, bestfitparmsdur[1, 3], bestfitparmsdur[2, 3],
                              bestfitparmsdur[3, 3], bestfitparmsdur[4, 3], bestfitparmsdur[5, 3], 
                              bestfitparmsdur[6, 3], bestfitparmsdur[7, 3], bestfitparmsdur[8, 3], 
                              bestfitparmsdur[9, 3], bestfitparmsdur[10, 3]), 
                      function(x) dnbinom(0:(2 * maxdelay), bestfitparmsdur[11, 3], mu = x))
  delayH2D <- sapply(mean_DT(1:tmax, bestfitparmsdur[1, 4], bestfitparmsdur[2, 4],
                             bestfitparmsdur[3, 4], bestfitparmsdur[4, 4], bestfitparmsdur[5, 4], 
                             bestfitparmsdur[6, 4], bestfitparmsdur[7, 4], bestfitparmsdur[8, 4], 
                             bestfitparmsdur[9, 4], bestfitparmsdur[10, 4]), 
                     function(x) dnbinom(0:(2 * maxdelay), bestfitparmsdur[11, 4], mu = x))
  
  delayA2IC2D <- t(sapply(1:(1 + 2 * maxdelay), function(x) colSums(delayA2IC[1:x,,drop=F] * delayIC2D[x:1,,drop = F])))
  delayA2IC2H2D <- t(sapply(1:(1 + 2 * maxdelay), function(x) colSums(delayA2IC2D[1:x,,drop=F] * delayH2D[x:1,,drop=F])))
  
  # normalise delay distributions
  delayA2IC <- normalisedelays(delayA2IC, maxdelay)
  delayA2D <- normalisedelays(delayA2D, maxdelay)
  delayIC2D <- normalisedelays(delayIC2D, maxdelay)
  delayH2D <- normalisedelays(delayH2D, maxdelay)
  delayA2IC2D <- normalisedelays(delayA2IC2D, maxdelay)
  delayA2IC2H2D <- normalisedelays(delayA2IC2H2D, maxdelay)
 
  return(list(
    delayA2IC = delayA2IC,
    delayA2D = delayA2D,
    delayIC2D = delayIC2D,
    delayH2D = delayH2D,
    delayA2IC2D = delayA2IC2D,
    delayA2IC2H2D = delayA2IC2H2D
  ))
}


##################################################
### Create lists with probabilities and delays ###
##################################################
NICEprobabilities <- get_allprobs(as.numeric(ANALYSISDATE - as.Date("2020-02-12")))
NICEdelays <- get_alldelays(as.numeric(ANALYSISDATE - as.Date("2020-02-12")))

################################
### Remove temporary objects ###
################################
rm(pddata, pddata4fitdurs,
   prob_admission_inD, 
   prob_death_inA, prob_death_untreated,
   prob_ic_age, mean_DT,
   normalisedelays,
   get_alldelays, get_allprobs,
   minloglik1, minloglik2, minloglik2a, minlogLik3,
   res1, res2, res2a, res3, 
   bestfitparmsdeath, bestfitparmsic,
   bestfitparmshosp, bestfitparmsdur
   )
