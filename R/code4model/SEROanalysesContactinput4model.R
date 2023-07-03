#######################################################
### contacts, infectivities,                        ###
###   susceptibilities, symptomatic fractions,      ###
###   natural history and control                   ###
#######################################################

data_serology_pico1 <- read_csv("data/sero/pico1data_synthetic.csv")

meansampleday1 <- data_serology_pico1 %>%
  filter(region < 6) %>%
  mutate(day = as.Date(algdatuminvul_r1, "%d%b%Y")) %>% 
  pull(day) %>% mean(na.rm = T) %>% round()

data_serology_pico2 <- read_csv("data/sero/pico2data_synthetic.csv")

meansampleday2 <- data_serology_pico2 %>%
  filter(region < 6) %>%
  mutate(day = as.Date(algdatuminvul_r2, "%d%b%Y")) %>% 
  pull(day) %>% mean(na.rm = T) %>% round()

# number of hospitalisations (weighted = scaled to the hospitalisation-to-death rate at t = -Inf)
aantallenhosp <-
  pddata4fit %>%
  filter(ti <= as.numeric(meansampleday2 - 7 - as.Date("2020-02-12")) & ti > 0) %>%
  mutate(weight = mapply(function(x,y) 1/NICEprobabilities$probSe2A[x, y], x = ti, y = age + 1)) %>%
  group_by(age) %>%
  summarise(aantal.weighted = sum(weight),
            aantal.unweighted = n()) 

# summarise serology data as seroprevalences
seroprev1 <- data_serology_pico1 %>%
  filter(region < 6) %>%
  mutate(lftgroep = floor(lftyear_r1/10),
         lftgroep = pmin(lftgroep, 8)) %>%
  group_by(lftgroep) %>%
  summarise(seropos = sum(pico1_pos_combi * pico1_gewicht_alle_var) * n() / sum(pico1_gewicht_alle_var),
            aantal = n()) 
seroprev2 <- data_serology_pico2 %>%
  filter(region < 6) %>%
  mutate(lftgroep = floor(lftyear_r2/10),
         lftgroep = pmin(lftgroep, 8)) %>%
  group_by(lftgroep) %>%
  summarise(seropos = sum(pico2_pos * pico2_gewicht) * n() / sum(pico2_gewicht),
            aantal = n()) 

data_serology <- list(
  serodays = c(pico1 = meansampleday1, pico2 = meansampleday2),
  pico1 = seroprev1,
  pico2 = seroprev2)

# probability from Infection to Severe (Severe = scaled to hospitalisation rate at t = -Inf)
PreAdmissionProbs <- list(
  probI2Se = 
    list(default = aantallenhosp$aantal.weighted/(agedist() * popsize() * seroprev2$seropos / seroprev2$aantal),
         samples = NULL)
)




##################################################
### Relative Susceptibility and Infectiousness ###
##################################################
# age-dependent susceptibility and infectivity

# calculate default vector based on pico2: weigh contact matrices before and during lockdown by 50% each
precontrolmatrix <- readRDS("data/contacts/ContactmatricesD3praktijk_midpoint_24mrt2020.rds")[['99']]
lockdownmatrix <- Reduce(`+`, readRDS("data/contacts/Contactmatrix_D3asEpiPose1_residualincreased_27mei2020.rds"))/200
expectedeigenvector <- (data_serology$pico2$seropos / data_serology$pico2$aantal) * popsize() * agedist()
expectedeigenvector <- expectedeigenvector / sum(expectedeigenvector)
agenormfit <- function(agenormvector) {
  ev1 <- abs(eigen(t(precontrolmatrix * c(exp(agenormvector), 1)) *
                     agedist() * c(exp(agenormvector), 1))$vectors[, 1])
  ev2 <- abs(eigen(t(lockdownmatrix * c(exp(agenormvector), 1)) *
                     agedist() * c(exp(agenormvector), 1))$vectors[, 1])
  ev <- 0.5 * (ev1/sum(ev1) + ev2/sum(ev2))
  sum(log(ev / expectedeigenvector)^2)
}
vectorfit <- c(exp(optim(rep(1,8), agenormfit, method = "CG")$par), 1)
vectorfit <- vectorfit/vectorfit[1]

ContactModelInput <- list(
  LogRelSusInf = list(default = log(vectorfit[2:9]),
                      samples = NULL)
)




#################################################################
### Infectiousness Profile (generation interval distribution) ###
#################################################################
# natural history parameters: discrete generation interval distribution, 
# based on SEEIIR-structure with mean generation interval of 4 days
ressim1 <- ode(y = c(1,0,0,0),
               times = seq(0,50,.01),
               func = function(t, y, parms) {list(c(0, parms[1] * y[1:3]) - parms[1] * y)},
               parms = c(0.875))
cdfgentime <- cumsum(rowSums(ressim1[,4:5]))/sum(ressim1[,4:5])
pdfgentime <- c(0, tail(cdfgentime, -1) - head(cdfgentime, -1))
pgentime <- function(gt) {
  gewichten <- c(seq(0,1,.01), rev(seq(0,.99,.01)))
  gewauc <- sum(gewichten)
  posities <- 100 * seq(gt - 1, gt + 1, .01) + 1
  gewichten <- gewichten[posities > 0]
  posities <- posities[posities > 0]
  100 * sum(gewichten * pdfgentime[posities])/gewauc
}
pgt <- sapply(0:20, pgentime)
ContactModelInput <- c(ContactModelInput,
                       list(
                         InfCurves = matrix(rep(rev(pgt[2:13]/sum(pgt[2:13])), each = 9), nrow = 9)
                       ))


# input for logy0 and infectivities
ContactModelInput <- c(ContactModelInput,
                       list(
                         LogY0 = list(samples = NULL),
                         LogInfectivities = list(samples = NULL),
                         InfectivityDeviations = list(samples = 1)
                       ))



rm(meansampleday1, meansampleday2,
   aantallenhosp, seroprev1, seroprev2,
   data_serology_pico1, data_serology_pico2,
   precontrolmatrix, lockdownmatrix,
   expectedeigenvector, agenormfit, vectorfit,
   ressim1, cdfgentime, pdfgentime,
   pgentime, pgt)

