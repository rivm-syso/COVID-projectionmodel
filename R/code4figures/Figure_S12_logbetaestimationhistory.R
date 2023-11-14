###########################################################
### Supplementary code with Klinkenberg et al (2023)
### Additional analysis: reconstruct fit to 40 datasets up to 6 January 2021
###########################################################

################################################################
### Load all results until model fitting (from Masterscript) ###
################################################################
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(lubridate)
library(deSolve)
library(cowplot)
options(dplyr.summarise.inform = FALSE)
library(RColorBrewer)
color_palette <- scales::hue_pal()(6)

ANALYSISDATE <- as.Date("2021-01-06")
load(paste0("data/processedmodelinput/modeldatainput", ANALYSISDATE, ".RData"))
source("R/code4model/Replace_syntheticresults.R")

#################
### function to reconstruct approximate IC admission data set of earlier analysis day, 
###   by multiplying incidence with reporting probability, and rounding
### Function uses global ANALYSISDATE
#################
inc_nice_ic_hist <- function() {
  reportingprobs <- c(rep(1, 228), rev(ReportingDelays$delayIC))
  inflated_inc <- inc_nice_ic() / reportingprobs
  shorter_inc <- head(inflated_inc, as.numeric(ANALYSISDATE - ymd("2020-02-12")))
  toreturn <- round(shorter_inc * tail(reportingprobs, as.numeric(ANALYSISDATE - ymd("2020-02-12"))))
  return(toreturn)
}

######################################
### Input data for all 40 analyses ###
######################################
### dates of the 40 analyses
ANALYSISDATES <- as.Date(c(paste0("2020-",
                                       c("03-28","03-30","04-03","04-10","04-16","04-24",
                                         "05-01","05-08","05-15","05-22","05-29","06-05",
                                         "06-12","06-19","06-27","07-03","07-10","07-17",
                                         "07-24","08-06","08-21","08-28","09-04","09-10",
                                         "09-18","09-25","10-02","10-08","10-16","10-23",
                                         "10-30","11-05","11-12","11-19","11-25","12-02",
                                         "12-09","12-16","12-28")), "2021-01-06"))
### contact matrices up to 6 January 2021
allcontactcontrols <- c("precontrol_mean", "intellockdown_mean", "intellockdown_mean", "batch1_mean",
                        "batch2_mean", "batch3_mean", "sep2020_mean", "okt2020_mean",
                        "okt2020holiday_mean", "partlockdown_mean", "nov2020_mean", "partlockdown_mean",
                        "winterlockdown_mean", "winterlockdownChristmas2020_mean")
### end dates for the contact matrices, as used on 6 January 2021 (transition days u_i)
allenddates<- c("2020-03-18", "2020-03-28", "2020-05-10", "2020-06-01",
                "2020-07-05", "2020-08-30", "2020-09-28", "2020-10-14",
                "2020-10-25", "2020-11-04", "2020-11-18", "2020-12-14",
                "2020-12-20")
### for each of the 40 analyses: by how many days are u_1 and u_2 different from 6 January 2021
enddateshifts <- sapply(c(1,1,1,2,2,2, 2,2,2,2,3,3, 3,4,5,5,5,5, 5,5,5,5,5,5, 5,5,5,5,5,6, 6,6,6,6,6,7, 7,7,7,7),
                        function(x) list(c(-5,-4,0,0,0,0,0,0,0,0,0,0,0),
                                         c(-3,-5,0,0,0,0,0,0,0,0,0,0,0),
                                         c(-3,-1,0,0,0,0,0,0,0,0,0,0,0),
                                         c(-3,7,0,0,0,0,0,0,0,0,0,0,0),
                                         c(-3,5,0,0,0,0,0,0,0,0,0,0,0),
                                         c(-1,0,0,0,0,0,0,0,0,0,0,0,0),
                                         c(0,0,0,0,0,0,0,0,0,0,0,0,0))[[x]])
### for each of the 40 analyses: what is the pattern of equal log(beta)s among the control periods (changepoints)
periodgroupsets <- sapply(c(1,1,1,1,2,2, 2,2,2,3,3,3, 4,4,4,4,4,5, 5,6,6,7,7,8, 8,9,10,10,10,11, 11,12,13,13,14,14, 15,15,16,16),
                          function(x) list(c(1,2,2),
                                           c(1,2,3),
                                           c(1,2,3,3),
                                           c(1,2,3,3,3),
                                           c(1,2,3,3,4,4),
                                           c(1,2,3,3,4,5),
                                           c(1,2,3,3,3,3),
                                           c(1,2,3,3,3,3,3),
                                           c(1,2,3,3,3,3,4),
                                           c(1,2,3,3,3,3,3,3),
                                           c(1,2,3,3,3,3,3,4,3),
                                           c(1,2,3,3,3,3,3,4,3,3),
                                           c(1,2,3,3,3,3,3,4,3,3,3),
                                           c(1,2,3,3,3,3,3,4,3,5,5,5),
                                           c(1,2,3,3,3,3,3,4,3,5,5,6),
                                           c(1,2,3,3,3,3,3,4,3,5,5,6,6,6))[[x]])

#######################
### all 40 analyses ###
#######################
### run the 40 analyses
llosfinal <- list()
for(i in 1:40) {
  # show progress
  cat(i, "\n")
  
  # change global ANALYSISDATE
  ANALYSISDATE <<- ANALYSISDATES[i]
  
  # define changepoints, enddates of control periods, matrices for control periods, and initial values for 'optim'
  periodgroups2use <- periodgroupsets[[i]]
  enddates2use <- c.Date(head(as.Date(allenddates) + enddateshifts[, i], length(periodgroups2use) - 1), ANALYSISDATE) 
  print(enddates2use) # check end dates
  contactcontrol2use <- head(allcontactcontrols, length(periodgroups2use))
  loginfectivities2use <- c(5.23,4.71,5.20,5.38,5.27,5.51,5.51)[1:max(periodgroups2use)]
  
  # fit the model as in 00_masterscript_20210106.R, but without hessian (variance matrix)
  llo <- logLik_optimise(logy0 = 1.68,
                         loginfectivities = loginfectivities2use,
                         logrelsusinf = "default",
                         probi2se = "default",
                         dataIC = inc_nice_ic_hist(),
                         enddates = enddates2use,
                         contactcontrol = contactcontrol2use,
                         reportprob = ReportingDelays$delayIC,
                         periodgroups = periodgroups2use,
                         estimate_relsusinf = F,
                         hessian = F, maxdelay = 30)
  
  # add to resultslist
  llosfinal <- c(llosfinal, list(llo))
}

saveRDS(llosfinal, "results/additionalresults/maxlikelihoodsestimationhistory_20210106.rds")


### put the results of the 40 analyses in a single tibble
resulttibble <- tibble()
for(x in 1:40) {
  resulttibble <- bind_rows(resulttibble, 
                            tibble(
                              analysis = x,
                              period = 1:(length(periodgroupsets[[x]]) + 1),
                              logbeta = llosfinal[[x]]$par[1 + c(periodgroupsets[[x]], tail(periodgroupsets[[x]], 1))]))
}

### load the data for the upper panel of Figure S2 (same as Figure 2b)
prognosisdata <- readRDS("data/figuredata/PrognosisData20210106.rds")
observations <- readRDS("data/figuredata/Observations20220321.rds")


#######################
### make Figure S12 ###
#######################
progplotdata <- prognosisdata %>%
  group_by(tijd, prognosisday) %>%
  summarise(mininc = quantile(incICU, .025),
            maxinc = quantile(incICU, .975),
            obsinc = median(observed)) %>%
  ungroup() %>%
  mutate(correct = obsinc > mininc & obsinc < maxinc) %>%
  arrange(prognosisday, tijd) %>%
  group_by(prognosisday) %>%
  mutate(correct = cumsum(correct) < max(cumsum(correct)),
         correct = lag(correct, default = T)) %>%
  ungroup() %>% 
  mutate(correct = if_else(correct, "Yes", "No"),
         correct = factor(correct, levels = c("Yes", "No")),
         analysis = match(prognosisday, ANALYSISDATES)) %>%
  rename(`Projection within \n95% interval` = correct,
         epidate = tijd)

FigS12a <- prognosisdata %>%
  filter(tijd <= prognosisday + 21, tijd >= prognosisday) %>%
  filter(prognosisday <= ANALYSISDATE) %>%
  group_by(tijd, prognosisday) %>%
  summarise(medinc = median(incICU),
            mininc = quantile(incICU, .025),
            maxinc = quantile(incICU, .975),
            obsinc = median(observed)) %>%
  ungroup() %>%
  mutate(correct = obsinc > mininc & obsinc < maxinc) %>%
  arrange(prognosisday, tijd) %>%
  group_by(prognosisday) %>%
  mutate(correct = cumsum(correct) < max(cumsum(correct)),
         correct = lag(correct, default = T)) %>%
  ungroup() %>% 
  mutate(correct = if_else(correct, "Yes", "No"),
         correct = factor(correct, levels = c("Yes", "No"))) %>%
  rename(`Projection within \n95% interval` = correct) %>%
  mutate(pointsize = T,
         progstart = if_else(tijd == prognosisday, medinc, NA_real_),
         punten = "Start") %>%
  ggplot(aes(x = tijd)) +
  # geom_ribbon(aes(ymin = mininc, ymax = maxinc, group = prognosisday), alpha = 0.2) +
  geom_point(data = observations %>% filter(tijd < ANALYSISDATE + 22) %>% mutate(punten = "Data"),
             aes(y = observed, size = punten)) +
             # aes(y = observed), size = 1) +
  geom_line(aes(y = medinc, colour = `Projection within \n95% interval`, group = prognosisday), 
            size = 1) +
  geom_point(aes(y = progstart, shape = punten), size = 3, colour = "black") +
  scale_color_manual(values = color_palette[c(4, 6)]) +
  scale_size_manual(values = c(1), labels = c("Data"), name = "") +
  scale_shape_manual(values = c(3), labels = c("Start of\nprojection"), name = "") +
  scale_x_date(limits = c(ymd("2020-02-12"), ANALYSISDATE + 22),
               expand = c(0, 0), date_breaks = "1 month", date_labels = "1 %b") +
  scale_y_continuous(limits = c(-1, NA),
                     expand = expansion(mult = c(0, 0.02))) +
  theme_light(base_size = 11) +
  theme(legend.position = c(0.4, 0.9), 
        legend.justification = c(0, 1),
        legend.box = "horizontal") +
  #theme(text = element_text(size = 6)) +
  labs(#title = "Daily ICU admissions: forecasts and data",
    x = "Date", 
    y = "ICU admissions") +
  guides(size = guide_legend(order = 2),color = guide_legend(order = 1)) 

FigS12b <- resulttibble %>%
  mutate(epidate = 1 + as.Date(allenddates[period]) + enddateshifts[cbind(pmin(13, period), analysis)]) %>%
  full_join(expand_grid(analysis = 1:40,
                        epidate = seq(ymd("2020-02-13"), ymd("2021-01-06"), by = "day"))) %>%
  mutate(analysisdate = ANALYSISDATES[analysis],
         period = if_else(is.na(period), 0L, period + 1L),
         period = if_else(epidate == as.Date("2020-02-13"), 1L, period)) %>%
  arrange(analysisdate, epidate) %>%
  group_by(analysisdate) %>%
  mutate(period = cummax(period)) %>%
  ungroup() %>%
  select(-logbeta) %>% 
  left_join(resulttibble) %>% 
  filter(epidate <= analysisdate) %>% 
  ggplot(aes(x = epidate, y = analysis)) +
  geom_tile(aes(fill = logbeta)) +
  scale_fill_gradientn(colours = c("yellow", "orange", "red"), values = scales::rescale(c(4.7, 5.20, 5.52))) +
  geom_line(data = progplotdata, aes(colour = `Projection within \n95% interval`, group = prognosisday), size = 1) +
  scale_x_date(limits = c(ymd("2020-02-12"), ANALYSISDATE + 22),
               expand = c(0, 0), date_breaks = "1 month", date_labels = "1 %b") +
  scale_color_manual(values = color_palette[c(4, 6)]) +
  theme_light() +
  theme(legend.position = c(.9,.3)) + 
  labs(x = "Day", y = "Model fit (order)", fill = "log(beta)")

plot_grid(FigS12a, FigS12b, nrow = 2, rel_heights = c(.3,.7))
ggsave("results/additionalresults/FigS12.pdf", bg = "white", width = 20, height = 25, units = "cm")
ggsave("results/additionalresults/FigS12.jpg", bg = "white", width = 20, height = 25, units = "cm", dpi = 600)

