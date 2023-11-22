#########################
### Code for Figure 1 ###
#########################
library(tidyverse)
library(lubridate)
library(deSolve)
library(cowplot)
library(ggforce)
library(scales)
library(RColorBrewer)

### From 00_masterscript_20210106.R
ANALYSISDATE <- ymd("2021-01-06")
load(paste0("data/processedmodelinput/modeldatainput", ANALYSISDATE, ".RData"))
source("R/code4model/Replace_syntheticresults.R")

color_palette <- hue_pal()(6)

################## Fig 1: model structure and delays ##################

# time point 6 Jan 2021, and age group [60-70)
time <- as.integer(ymd("2021-01-06") - ymd("2020-02-12"))
age_group <- 7


# Structure transmission model

rect_left <- c(3, 8, 12, 17, 21, 26)
rect_size <- 3
y_mid <- 3


p1A <- ggplot() +
  geom_rect(data = tibble(xmin = rect_left) %>% 
              mutate(xmax = xmin + rect_size,
                     ymin = rep(y_mid - 0.5*rect_size, n()),
                     ymax = ymin + rect_size),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            col = 1, 
            fill = "white") +
  geom_segment(data = tibble(x = rect_left[-which.max(rect_left)] + rect_size,
                             xend = rect_left[-which.min(rect_left)],
                             y = y_mid,
                             yend = y_mid,
                             type = factor(c(0, rep(0, length(rect_left)-2)))),
               aes(x = x, xend = xend, y = y, yend = yend, color = type),
               arrow = arrow(length = unit(0.15, "npc")),
               linewidth = 1) +
  geom_segment(data = tibble(x = c(mean(rect_left[1:2] + c(rect_size, 0)),
                                   rect_left[4] + rect_size/2,
                                   rect_left[5] + rect_size/2),
                             xend = c(rect_left[5] + rect_size/2,
                                      rect_left[4] + rect_size/2,
                                      rect_left[5] + rect_size/2),
                             y = y_mid + rect_size,
                             yend = y_mid + rect_size / c(1, 2, 2)),
               aes(x = x, xend = xend, y = y, yend = yend),
               linetype = 2, color = color_palette[1],
               linewidth = 0.7) +
  geom_segment(data = tibble(x = c(mean(rect_left[1:2] + c(rect_size, 0))),
                             xend = c(mean(rect_left[1:2] + c(rect_size, 0))),
                             y = y_mid + rect_size,
                             yend = y_mid + rect_size/5),
               aes(x = x, xend = xend, y = y, yend = yend),
               linetype = 2, color = color_palette[1],
               arrow = arrow(length = unit(0.15, "npc")),
               linewidth = 0.7) +
  geom_text(data = tibble(x = c(rect_left + 0.5*rect_size, mean(rect_left[c(1,2)] + 0.5*rect_size)),
                          y = c(rep(y_mid, length(rect_left)), y_mid - 0.5*rect_size),
                          label = c("S", "~E^{1}", "~E^{2}", "~I^{1}", "~I^{2}", "R", "y")),
            aes(x = x, y = y, label = label),
            size = 7,
            parse = TRUE) +
  scale_color_manual(values = c("1" = color_palette[1], "0" = 1)) +
  scale_y_continuous(limits = c(y_mid - 0.75*rect_size, NA)) +
  scale_x_continuous(limits = c(-1, 30)) +
  coord_equal() +
  labs(x = NULL,
       y = NULL) +
  guides(color = guide_none()) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank())


# Structure clinical outcome model

circle_mid <- c(3, 6, 9)
circle_rad <- 1


p1C <- ggplot() +
  geom_circle(data = tibble(x = circle_mid,
                            y = y_mid,
                            r = circle_rad),
              aes(x0 = x, y0 = y, r = r)) +
  geom_segment(data = tibble(x = c(1, circle_mid[-which.max(circle_mid)] + circle_rad, circle_mid),
                             xend = c(circle_mid - circle_rad, circle_mid),
                             y = c(rep(y_mid, 3), rep(y_mid - circle_rad, 3)),
                             yend = c(rep(y_mid, 3), rep(y_mid - 2*circle_rad, 3)),
                             type = factor(c(2, 3, 4, 5, 4, 6))),
               aes(x = x, xend = xend, y = y, yend = yend, color = type),
               arrow = arrow(length = unit(0.1, "npc")),
               size = 1) +
  geom_text(data = tibble(x = c(0.5, circle_mid),
                          y = y_mid,
                          label = c("y", "A", "ICU", "H")),
            # label = c("y", "H[1]", "ICU", "H[2]")),
            aes(x = x, y = y, label = label),
            size = 7,
            parse = TRUE) +
  scale_x_continuous(limits = c(-2, 12)) +
  scale_color_manual(values = c("2" = color_palette[4], 
                                "3" = color_palette[2],
                                "4" = color_palette[6],
                                "5" = color_palette[3],
                                "6" = color_palette[5])) +
  coord_equal() +
  labs(x = NULL,
       y = NULL) +
  guides(color = guide_none()) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank())



# Generation interval distribution
ressim1 <- ode(y = c(1,0,0,0),
               times = seq(0,50,.01),
               func = function(t, y, parms) {list(c(0, parms[1] * y[1:3]) - parms[1] * y)},
               parms = c(0.875))

p1B <- tibble(genint = 1:12, 
              probability = ContactModelInput$InfCurves[1, ] %>% rev) %>%
  ggplot(aes(x = genint, y = probability)) +
  geom_bar(stat = "identity",
           width = 1, 
           fill = color_palette[1]) +
  geom_line(data=tibble(genint = ressim1[1:1201,1], 
                        probability = 100*rowSums(ressim1[1:1201,4:5])/sum(ressim1[,4:5]))) +
  scale_x_continuous(breaks = 1:12) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_light() +
  labs(x = "Generation interval (days)", y = "Probability")


# Infection to admission distribution (for age group 60-69)

p1D <- ggplot(data = tibble(x = 0:100,
                            y = PreAdmissionDelays$delayI2A[, time, age_group]),
              aes(x = x, y = y)) +
  geom_bar(stat = "identity",
           width = 1, 
           fill = color_palette[4]) +
  coord_cartesian(xlim = c(0, 31)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(subtitle = "Infection to hospital admission (age 60-69)",
       x = "Delay (days)", 
       y = "Probability") +
  theme_light()


# Admission to ICU distribution

p1F <- ggplot(data = tibble(x = 0:100,
                            y = NICEdelays$delayA2IC[, time + 1]),
              aes(x = x, y = y)) +
  geom_bar(stat = "identity",
           width = 1, 
           fill = color_palette[2]) +
  coord_cartesian(xlim = c(-0.5, 21)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(subtitle = "Hospital to ICU admission (all ages)",
       x = "Delay (days)", 
       y = "Probability") +
  theme_light()


# Lenghts of stay
# delayA2D: length of stay H1 (not to ICU)
# delayH2D: length of stay H2
# delayIC2D: length of stay ICU


p1E <- bind_rows(
  "A" = tibble(probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
               d = sapply(probs, function(q) which(cumsum(NICEdelays$delayA2D[, time + 1]) > q)[1] - 1),
               mean = sum((0:100)*NICEdelays$delayA2D[, time + 1]),
               type = "5"),
  "ICU" = tibble(probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                 d = sapply(probs, function(q) which(cumsum(NICEdelays$delayIC2D[, time + 1]) > q)[1] - 1),
                 mean = sum((0:100)*NICEdelays$delayIC2D[, time + 1]),
                 type = "4"),
  "H" = tibble(probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
               d = sapply(probs, function(q) which(cumsum(NICEdelays$delayH2D[, time + 1]) > q)[1] - 1),
               mean = sum((0:100)*NICEdelays$delayH2D[, time + 1]),
               type = "6"),
  .id = "delay") %>% 
  mutate(delay = factor(delay, levels = c("H", "ICU", "A"), labels = c("H", "ICU", "A\n(not to ICU)"))) %>% 
  pivot_wider(names_from = probs, values_from = d, names_prefix = "p_") %>% 
  ggplot(aes(y = delay, yend = delay, col = type)) +
  geom_segment(aes(x = p_0.05, xend = p_0.95),
               alpha = 0.5,
               lwd = 2) +
  geom_segment(aes(x = p_0.25, xend = p_0.75),
               lwd = 2) +
  geom_point(aes(x = mean),
             size = 4) +
  geom_segment(aes(x = p_0.5 - 0.1, xend = p_0.5 + 0.1),
               col = "white", 
               lwd = 4) +
  scale_color_manual(values = c("2" = color_palette[4], 
                                "3" = color_palette[2],
                                "4" = color_palette[6],
                                "5" = color_palette[3],
                                "6" = color_palette[5])) +
  labs(x = "Length of stay (days)",
       y = NULL) +
  guides(col = guide_none()) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))


plot_grid(p1A, p1B, p1C, p1D, p1E, p1F,
          nrow = 3,
          rel_widths = c(1, 0.8),
          labels = c("a", "b", "c", "d", "e", "f"))

ggsave("results/Fig1.pdf", bg = "white", width = 20, height = 16, units = "cm")
ggsave("results/Fig1.jpg", bg = "white", width = 20, height = 16, units = "cm", dpi = 600)


