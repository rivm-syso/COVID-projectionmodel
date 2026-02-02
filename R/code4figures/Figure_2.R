#########################
### Code for Figure 2 ###
#########################
library(tidyverse)
library(lubridate)
library(deSolve)
library(cowplot)
library(ggforce)
library(scales)
library(RColorBrewer)
color_palette <- hue_pal()(6)

### From 00_masterscript_20210106.R
ANALYSISDATE <- ymd("2021-01-06")
load(paste0("data/processedmodelinput/modeldatainput", ANALYSISDATE, ".RData"))
source("R/code4model/Replace_syntheticresults.R")
llo <- readRDS("results/maxlikelihood_20210106.rds")
prognosisdata <- readRDS("data/figuredata/PrognosisData20210106.rds")
observations <- readRDS("data/figuredata/Observations20220321.rds")

# dates
enddates <- c(as.Date(c("2020-03-18", "2020-03-28", "2020-05-10", "2020-06-01",
                        "2020-07-05", "2020-08-30", "2020-09-28", "2020-10-14",
                        "2020-10-25", "2020-11-04", "2020-11-18", "2020-12-14",
                        "2020-12-20")), ANALYSISDATE)
startdates <- c(ymd("2020-02-12"), head(enddates, -1) + 1)

# values for infectivity (estimates in llo$par)
priormeans <- sapply(c("precontrol_mean", "intellockdown_mean", "intellockdown_mean", "batch1_mean",
                       "batch2_mean", "batch3_mean", "sep2020_mean", "okt2020_mean",
                       "okt2020holiday_mean", "partlockdown_mean", "nov2020_mean", "partlockdown_mean",
                       "winterlockdown_mean", "winterlockdownChristmas2020_mean"), 
                     function(x) eigen(t(ContactMatrices[[x]] * exp(c(0, ContactModelInput$LogRelSusInf$default))) * 
                                         exp(c(0, ContactModelInput$LogRelSusInf$default)))$values[1]) / 
  eigen(t(ContactMatrices$precontrol_mean * exp(c(0, ContactModelInput$LogRelSusInf$default))) *
          exp(c(0, ContactModelInput$LogRelSusInf$default)))$values[1]

posteriormeans <- sapply(c("precontrol_mean", "intellockdown_mean", "intellockdown_mean", "batch1_mean",
                           "batch2_mean", "batch3_mean", "sep2020_mean", "okt2020_mean",
                           "okt2020holiday_mean", "partlockdown_mean", "nov2020_mean", "partlockdown_mean",
                           "winterlockdown_mean", "winterlockdownChristmas2020_mean"), 
                         function(x) eigen(t(ContactMatrices[[x]] * exp(c(0, ContactModelInput$LogRelSusInf$default))) *
                                             exp(c(0, ContactModelInput$LogRelSusInf$default)))$values[1]) * 
  exp(llo$par[2:7][
    c(1,2,3,3, 3,3,3,4, 3,5,5,6, 6,6)
  ]) / exp(llo$par[2]) / eigen(t(ContactMatrices$precontrol_mean * exp(c(0, ContactModelInput$LogRelSusInf$default))) *
                            exp(c(0, ContactModelInput$LogRelSusInf$default)))$values[1]


# dates and infectivities in one dataframe
infectivitytibble <- tibble(
  tdate = c(startdates, enddates),
  start = c(rep(T, 14), rep(F, 14)),
  priormean = rep(priormeans, 2),
  posteriormean = rep(posteriormeans, 2)
) %>% 
  arrange(tdate) 


change_points <- infectivitytibble %>% 
  mutate(ratio = priormean/posteriormean) %>% 
  filter(abs(c(diff(ratio), 0)) > 0.001)


p2A <- infectivitytibble %>%
  pivot_longer(priormean:posteriormean, names_to = "Estimate", values_to = "mean") %>%
  mutate(Estimate = factor(Estimate, levels = c("priormean", "posteriormean"), labels = c("before fitting", "after fitting"))) %>% 
  filter(tdate < ymd("2021-01-08")) %>%
  ggplot(aes(x = tdate, y = mean, color = Estimate)) + 
  geom_vline(data = tibble(tdate = change_points$tdate),
             aes(xintercept = tdate),
             lty = 2,
             color = "darkgrey") +
  geom_line() +
  geom_segment(aes(y = y0, yend = y1, xend = tdate), 
               data = tibble(tdate = enddates, 
                             y0 = 0.05, 
                             y1 = 0), 
               color = "black",
               arrow = arrow(length = unit(.2, "cm"))) +
  geom_text(aes(y = y0, label = lab), 
            data = tibble(tdate = enddates, 
                          y0 = 0.1, 
                          lab = c("u[1]","u[2]","u[3]","u[4]",
                                  "u[5]","u[6]","u[7]","u[8]",
                                  "u[9]","u[10]","u[11]","u[12]","u[13]","u[14]")), 
            color = "black", parse = T) +
  scale_color_manual(values = color_palette[c(1, 3)]) +
  scale_x_date(limits = c(ymd("2020-02-12"), ANALYSISDATE + 22),
               expand = c(0, 0), date_breaks = "1 month", date_labels = "1 %b") +
  scale_y_continuous(limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.02))) +
  # theme(legend.position = c(.42,.75), legend.direction = "horizontal") +
  #theme(legend.position = c(.9,.8), legend.spacing.y = unit(.3, "cm"), legend.key.height = unit(.2, "cm")) +
  labs(y = "Relative transmission rate", x = "Date") +
  theme_light(base_size = 11) +
  theme(legend.position = c(0.14, 0.95),
        legend.justification = c(0, 1))



p2B <- prognosisdata %>%
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
  geom_line(aes(y = medinc, colour = `Projection within \n95% interval`, group = prognosisday), 
            size = 1) +
  geom_point(aes(y = progstart, shape = punten), size = 3, colour = "black") +
  scale_color_manual(values = color_palette[c(4, 6)]) +
  scale_size_manual(values = c(1), labels = c("Data"), name = "") +
  scale_shape_manual(values = c(3), labels = c("Starts of\nprojections"), name = "") +
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


plot_grid(p2A, p2B,
          nrow = 2,
          rel_heights = c(0.8, 1),
          labels = c("a", "b"),
          align = "h")

ggsave("results/Fig2.pdf", bg = "white", width = 20, height = 16, units = "cm")
ggsave("results/Fig2.jpg", bg = "white", width = 20, height = 16, units = "cm", dpi = 600)


p2A
ggsave("results/additionalresults/FigPoster3.jpg", bg = "white", width = 20, height = 8, units = "cm", dpi = 600)
p2B
ggsave("results/additionalresults/FigPoster4.jpg", bg = "white", width = 20, height = 8, units = "cm", dpi = 600)
