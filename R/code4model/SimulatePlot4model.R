samplesforsimulation <- function(optimresult, nrof_samples = 200, dev_futureinfectivity = 0.02,
                                 estimate_relsusinf = FALSE, estimate_logy0 = TRUE, logy0 = NULL) {
  if(any(diag(solve(optimresult$hessian)) < 0)) {
    stop("cannot sample: negative variances")
  }
  
  set.seed(ANALYSISDATE)
  
  randompars <- mvtnorm::rmvnorm(nrof_samples, optimresult$par, solve(optimresult$hessian))
  ContactModelInput$InfectivityDeviations$samples <<- runif(nrof_samples, min = 1 - dev_futureinfectivity, max = 1 + dev_futureinfectivity)

  if(!estimate_logy0) {
    randompars <- cbind(logy0, randompars)
  }
  
  ContactModelInput$LogY0$samples <<- randompars[, 1]
  
  if(estimate_relsusinf) {
    ContactModelInput$LogRelSusInf$samples <<- randompars[, 2:9]
    ContactModelInput$LogInfectivities$samples <<- randompars[, -(1:9), drop = FALSE]
  } else {
    ContactModelInput$LogInfectivities$samples <<- randompars[, -1, drop = FALSE]
  }
}




# provide scenarionames_in_plot as named vector, where names are original names, and the vector
# contains the new names. In short this is done like
# c(Scenario_1 = "new name of first scenario", Scenario_2 = "new name of second scenario")
# provide them in the order to be shown in the legend
plotICinc <- function(simulationoutput, 
                      from = as.Date("2020-03-08"), to,
                      y_axis_max = 150,
                      scenarios = "all", 
                      scenarionames_in_plot = "as_is") {
  if(scenarios[1] != "all") {
    simulationoutput <- simulationoutput %>%
      filter(scenario %in% scenarios)
  }
  
  if(scenarionames_in_plot[1] != "as_is") {
    simulationoutput <- simulationoutput %>%
      mutate(scenario = scenarionames_in_plot[scenario],
             scenario = factor(scenario, scenarionames_in_plot))
  }
  
  simulationoutput %>% 
    mutate(tijd = time + as.Date("2020-02-12")) %>%
    filter(tijd >= from & tijd <= to) %>%
    group_by(tijd, scenario, repl) %>%
    summarise(
      totsimulated = sum(incICU)
    ) %>%
    group_by(tijd, scenario) %>%
    summarise(
      lowerbound = qpois(0.025, quantile(totsimulated, .025)),
      upperbound = qpois(0.975, quantile(totsimulated, .975)),
      midbound = median(totsimulated)
    ) %>%
    ungroup() %>%
    left_join(tibble(tijd = seq(as.Date("2020-02-13"), ANALYSISDATE, 1),
                     datapoints = inc_nice_ic()), 
              by = c("tijd")) %>%
    mutate(datapunten = "NICE IC-opnames") %>%
    ggplot(aes(x = tijd, y = midbound)) +
    geom_line(aes(color = scenario), size = 2) +
    geom_ribbon(aes(ymin = lowerbound, ymax = upperbound, fill = scenario, color = NULL), alpha = 0.3) +
    geom_point(mapping = aes(y = datapoints, shape = datapunten), size = 2) +
    theme_light() + 
    labs(x = "Dag", y = "Dagelijks aantal IC-opnames") +
    labs(title = "Aantal IC-opnames per dag") +
    scale_color_discrete(name = "scenario") +
    scale_fill_discrete(name = "scenario") +
    scale_x_date(date_breaks = "1 month", date_labels = "1 %b") +
    coord_cartesian(ylim = c(0, y_axis_max)) +
    guides(color = guide_legend(order = 1), 
           fill = guide_legend(order = 1),
           shape = guide_legend(order = 2))
}





plotICprev <- function(simulationoutput, 
                      from = as.Date("2020-03-08"), to,
                      y_axis_max = 1300,
                      scenarios = "all", 
                      scenarionames_in_plot = "as_is") {
  if(scenarios[1] != "all") {
    simulationoutput <- simulationoutput %>%
      filter(scenario %in% scenarios)
  }
  
  if(scenarionames_in_plot[1] != "as_is") {
    simulationoutput <- simulationoutput %>%
      mutate(scenario = scenarionames_in_plot[scenario],
             scenario = factor(scenario, scenarionames_in_plot))
  }
  
  simulationoutput %>% 
    mutate(tijd = time + as.Date("2020-02-12")) %>%
    filter(tijd >= from & tijd <= to) %>%
    group_by(tijd, scenario, repl) %>%
    summarise(
      totsimulated = sum(prevICU)
    ) %>%
    group_by(tijd, scenario) %>%
    summarise(
      lowerbound = qpois(0.025, quantile(totsimulated, .025)),
      upperbound = qpois(0.975, quantile(totsimulated, .975)),
      midbound = median(totsimulated)
    ) %>%
    ungroup() %>%
    left_join(tibble(tijd = seq(as.Date("2020-02-13"), ANALYSISDATE, 1),
                     datapoints = prev_nice_ic()), 
              by = c("tijd")) %>%
    mutate(datapunten = "NICE IC-bezetting") %>%
    ggplot(aes(x = tijd, y = midbound)) +
    geom_line(aes(color = scenario), size = 2) +
    geom_ribbon(aes(ymin = lowerbound, ymax = upperbound, fill = scenario, color = NULL), alpha = 0.3) +
    geom_point(mapping = aes(y = datapoints, shape = datapunten), size = 2) +
    theme_light() + 
    labs(x = "Dag", y = "IC-bezetting") +
    labs(title = "Aantal bezette IC-bedden") +
    scale_color_discrete(name = "scenario") +
    scale_fill_discrete(name = "scenario") +
    scale_x_date(date_breaks = "1 month", date_labels = "1 %b") +
    coord_cartesian(ylim = c(0, y_axis_max)) +
    guides(color = guide_legend(order = 1), 
           fill = guide_legend(order = 1),
           shape = guide_legend(order = 2))
}


plothospinc <- function(simulationoutput, 
                       from = as.Date("2020-03-08"), to,
                       y_axis_max = 600,
                       scenarios = "all", 
                       scenarionames_in_plot = "as_is") {
  if(scenarios[1] != "all") {
    simulationoutput <- simulationoutput %>%
      filter(scenario %in% scenarios)
  }
  
  if(scenarionames_in_plot[1] != "as_is") {
    simulationoutput <- simulationoutput %>%
      mutate(scenario = scenarionames_in_plot[scenario],
             scenario = factor(scenario, scenarionames_in_plot))
  }
  
  simulationoutput %>% 
    mutate(tijd = time + as.Date("2020-02-12")) %>%
    filter(tijd >= from & tijd <= to) %>%
    group_by(tijd, scenario, repl) %>%
    summarise(
      totsimulated = sum(inchosp)
    ) %>%
    group_by(tijd, scenario) %>%
    summarise(
      lowerbound = qpois(0.025, quantile(totsimulated, .025)),
      upperbound = qpois(0.975, quantile(totsimulated, .975)),
      midbound = median(totsimulated)
    ) %>%
    ungroup() %>%
    left_join(tibble(tijd = seq(as.Date("2020-02-13"), ANALYSISDATE, 1),
                     datapoints = inc_nice_hosp()), 
              by = c("tijd")) %>%
    mutate(datapunten = "NICE ziekenhuisopnames") %>%
    ggplot(aes(x = tijd, y = midbound)) +
    geom_line(aes(color = scenario), size = 2) +
    geom_ribbon(aes(ymin = lowerbound, ymax = upperbound, fill = scenario, color = NULL), alpha = 0.3) +
    geom_point(mapping = aes(y = datapoints, shape = datapunten), size = 2) +
    theme_light() + 
    labs(x = "Dag", y = "Dagelijks aantal ziekenhuisopnames") +
    labs(title = "Aantal ziekenhuisopnames per dag") +
    scale_color_discrete(name = "scenario") +
    scale_fill_discrete(name = "scenario") +
    scale_x_date(date_breaks = "1 month", date_labels = "1 %b") +
    coord_cartesian(ylim = c(0, y_axis_max)) +
    guides(color = guide_legend(order = 1), 
           fill = guide_legend(order = 1),
           shape = guide_legend(order = 2))
}



plothospprev <- function(simulationoutput, 
                        from = as.Date("2020-03-08"), to,
                        y_axis_max = 4500,
                        scenarios = "all", 
                        scenarionames_in_plot = "as_is") {
  if(scenarios[1] != "all") {
    simulationoutput <- simulationoutput %>%
      filter(scenario %in% scenarios)
  }
  
  if(scenarionames_in_plot[1] != "as_is") {
    simulationoutput <- simulationoutput %>%
      mutate(scenario = scenarionames_in_plot[scenario],
             scenario = factor(scenario, scenarionames_in_plot))
  }
  
  simulationoutput %>% 
    mutate(tijd = time + as.Date("2020-02-12")) %>%
    filter(tijd >= from & tijd <= to) %>%
    group_by(tijd, scenario, repl) %>%
    summarise(
      totsimulated = sum(prevhosp)
    ) %>%
    group_by(tijd, scenario) %>%
    summarise(
      lowerbound = qpois(0.025, quantile(totsimulated, .025)),
      upperbound = qpois(0.975, quantile(totsimulated, .975)),
      midbound = median(totsimulated)
    ) %>%
    ungroup() %>%
    left_join(tibble(tijd = seq(as.Date("2020-02-13"), ANALYSISDATE, 1),
                     datapoints = prev_nice_hosp()), 
              by = c("tijd")) %>%
    mutate(datapunten = "NICE ziekenhuisbezetting") %>%
    ggplot(aes(x = tijd, y = midbound)) +
    geom_line(aes(color = scenario), size = 2) +
    geom_ribbon(aes(ymin = lowerbound, ymax = upperbound, fill = scenario, color = NULL), alpha = 0.3) +
    geom_point(mapping = aes(y = datapoints, shape = datapunten), size = 2) +
    theme_light() + 
    labs(x = "Dag", y = "ziekenhuisbezetting") +
    labs(title = "Aantal bezette ziekenhuisbedden (inclusief IC)") +
    scale_color_discrete(name = "scenario") +
    scale_fill_discrete(name = "scenario") +
    scale_x_date(date_breaks = "1 month", date_labels = "1 %b") +
    coord_cartesian(ylim = c(0, y_axis_max)) +
    guides(color = guide_legend(order = 1), 
           fill = guide_legend(order = 1),
           shape = guide_legend(order = 2))
}


