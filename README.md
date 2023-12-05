|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
|------------------------------------------------------------------------|
| BELANGRIJK: Deze repository bevat synthetische data, met een andere licentie dan de code (CC BY-NC-ND, http://creativecommons.org/licenses/by-nc-nd/4.0/?). Deze synthetische data zijn gemaakt met datasets met individuele persoonsgegevens, waar de echte personen zijn vervangen door 'neppersonen'. Dat is gedaan op zo'n manier dat analyse van de synthetische data ongeveer dezelfde resultaten geeft als analyse van de oorspronkelijke data. Als echte personen zich toch lijken te herkennen in de synthetische data, dan is dat zuiver toeval: wij publiceren geen persoonsgegevens. Statistische eigenschappen van de synthetische data anders dan de resultaten in deze repository zijn niet noodzakelijkerwijs gelijk aan statistische eigenschappen van de oorspronkelijke data. Deze synthetische data mogen niet buiten deze repository gebruikt worden. |
| IMPORTANT: This repository contains synthetic data, with a different license from the code (CC BY-NC-ND, http://creativecommons.org/licenses/by-nc-nd/4.0/?). These synthetic data have been created from datasets with personal information, but with real persons replaced by 'fake persons'. That has been done in such a way that the analysis of the synthetic data gives approximately the same results as analysis of the original data. If real persons seem to recognise themselves in the synthetic data, then that is pure coincidence: we do not publish personal data. Statistical properties of the synthetic data other than the results in this repository do not necessarily reflect the statistical properties of the original data. These synthetic data may not be used outside this repository.                                                          |

# Introduction

This repository contains the code used for data analysis and simulations leading to projections of COVID-19 ICU and hospital admissions, carried out by the Dutch National Institute for Public Health and the Environment (RIVM) on 6 January 2021. It is supplementary material with the publication "Projecting COVID-19 intensive care admissions in the Netherlands for policy advice: February 2020 to January 2021", by Klinkenberg et al (https://doi.org/10.1101/2023.06.30.23291989).

The code was originally used with the newest surveillance data, containing privacy-sensitive information on individual patient level. To let the code in this repository work, synthetic datasets have been created that approximately produce the same parameter estimates. To run the exact same simulations, the original parameter estimates have been provided as well.

## How to use the code?

Save the repository in your local environment and open it as an R Project in RStudio. By opening the file "R/00_masterscript_20210106.R" and running it in order, all analyses are carried out and all simulation functions are defined by sourcing code files elsewhere in the repository, all simulations are run, and results are plotted. The code files themselves contain (brief) comments explaining what is done. The masterscript consists of the following steps:

### Block 1: load libraries

-   The repository was created with the following version of R and packages

```         
platform       x86_64-redhat-linux-gnu     
arch           x86_64                      
os             linux-gnu                   
system         x86_64, linux-gnu           
status                                     
major          4                           
minor          3.0                         
year           2023                        
month          04                          
day            21                          
svn rev        84292                       
language       R                           
version.string R version 4.3.0 (2023-04-21)
nickname       Already Tomorrow 

  deSolve lubridate   forcats   stringr     dplyr     purrr     readr 
   "1.35"   "1.9.2"   "1.0.0"   "1.5.0"   "1.1.2"   "1.0.1"   "2.1.4" 
    tidyr    tibble   ggplot2 tidyverse     stats  graphics grDevices 
  "1.3.0"   "3.2.1"   "3.4.2"   "2.0.0"   "4.3.0"   "4.3.0"   "4.3.0" 
    utils  datasets   methods      base 
  "4.3.0"   "4.3.0"   "4.3.0"   "4.3.0" 
```

### Block 2: analyse data

-   Define ANALYSISDATE and DELAYSTARTDATE. The latter determines which data to include to calculate the reporting delay distribution.
-   Read population data (size and age distribution)
-   Read and analyse the NICE hospital data for all probabilities and lengths-of-stay distributions. Warning messages (NaNs produced) can be ignored, as the final fits look well.
-   Calculate the reporting delay distribution from the NICE data
-   Read and analyse the OSIRIS notification data for the symptom-to-hospital distributions, and define incubation period distribution based on literature\
-   Analyse serological data to estimate hospitalisation probabilities and age-dependent infectivity/susceptibility, and define generation interval distribution
-   Read the contact matrices for all sets of control measures
-   Save the results (the repository contains the saved file)

### Block 3: define additional functions

-   Functions to process simulation options and parameters to correct population and infectivity
-   Functions to process simulation options and parameters for delays and probabilities
-   Functions for the simulations themselves, with different purposes (estimation, simulation). There are two sets of simulation functions: first, 'engine' functions that do the actual simulations and take numerical input for all parameters; second, functions to optimise the likelihood or simulate scenarios, with more readable input parameters (eg names of matrices, named options). The file starts with function definitions to convert some of the readable input of the second set of functions to numerical input required for the first set.
-   Some additional functions required for parameter sampling and plotting simulated output

### Block 4: load newest incidence data

-   The ANALYSISDATE is redefined and the newest data imported. This makes it possible to fit the model and run simulations without re-estimating all datasets in step 1. Not needed here, but used during code updates.
-   Save everything just before fitting to incidence data (the repository contains the saved file)

### Block 5: load original results (only relevant when using synthetic data, not relevant when using original data as in the original code)

-   Replace parameter estimates from synthetic data with original parameter estimates from actual data. Keep synthetic data in memory.

### Block 6: optimise likelihood to estimate stepwise constant transmissibilities and initial state (step 3 of parameter estimation)

-   The function logLik_optimise() estimates the initial state y_0, and the stepwise constant transmissibilities, given the changepoints of these constant transmissibilities, by fitting the simulated daily ICU admissions to reported daily ICU admissions. The changepoints are indicated in the input vector 'periodgroups', containing indicator variables for which stepwise transmissibility is used with which contact matrix. The corresponding contact matrices are given in 'contactcontrol', with corresponding transition times in 'endtimes'. The function logLik_optimise is used to assess the likelihood per set of changepoints and matching transmissibilities. Many options for plausible changepoint sets are evaluated manually and an optimal set is selected based on AIC. When comparing two fits, the fit with fewer parameters is preferred if 2\*negative_log_likelihood + 2\*nr_of_changepoints is no more than two points higher. The negative_log_likelihoods of three fits (sets of changepoints) are given below, with the final model choice. Not all evaluated changepoint sets  are listed here.
-   Save the optim result (the repository contains the saved file)

### Block 7: run the simulations

-   Sample 200 parameter sets of initial values and stepwise transmissibilities given the result of the previous block
-   Run the simulations for different scenarios, given by the contact matrices. Median contact matrices (that were for fitting the simulated to the observed ICU admissions) are used up to 14 days before ANALYSISDATE. The uncertainty about the recent and future contact patterns is larger than the uncertainty about the past, so 200 different contact matrix samples are used after ANALYSISDATE - 14.
-   Save the simulations (the repository contains the saved file of the first 20 runs)

### Block 8: plotting the results

-   Some plots as used to present results to the policy makers.



## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
