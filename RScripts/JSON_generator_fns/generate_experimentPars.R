################################################################################
#                                                                              #
################################################################################
if (Sys.info()[1] == "Windows"){
  working_directory <- "D:/Polybox/robot-data/Robot-Shared/TreatmentStrategies"
} else {
working_directory <- here::here()
}
setwd(working_directory)
library(tibble)
library(jsonlite)

################################################################################
# Set Experiment Parameters                                                    #
################################################################################
antibiotics <- data_frame(
  drug_letter = c("A","B","AB"),
  drug_short = c("Nal","Str","NalStr"),
  drug_long = c("nalidixic acid", "streptomycin", "Nal + Strep"),
  MIC_wt = c(5,5,5),
  C_treat = c(20,12.5))

format <- data_frame(
  format = 384,
  n_rep = 4,
  n_rwells = format / n_rep,
  n_bl = 2,
  n_patients = n_rwells - n_bl,
  rwells_bl = list(1:n_bl),
  rwells_p = list((n_bl+1):n_rwells))

treatments <- data_frame(
  plate = 1:6,
  conditions_name = c(
    "no treatment",
    "monotherapy",
    "monotherapy",
    "combination therapy",
    "cycling",
    "mixing"),
  drugs = list(
    "none",     # ctrl, no treatment
    "A",        # treatment with drug A, nalidixic acid @ ~MIC
    "B",        # treatment with drug B, streptomycin @ ~MIC
    c("AB"),    # treatment with drug A+B, both @ ~MIC
    c("A", "B"),    # treatment with drug A for period days, then B for period days
    c("A", "B")     # treat 50% of Pop. with A, 50% with B
    ),
  concentration = list( # not used (yet?)
    1,
    1,
    1,
    1,
    c(1, 1),
    c(1, 1)), 
  period = c(1, 1, 1, 1, 2, 1), # number of days on one drug
  probability = list(1, 1, 1, 1, 1, c(0.5, 0.5)) # probability of receiving A/B
  )

# frequency of the strains in the environment, adds up to n_patients
community <- data_frame(
  strain = c("AB_r","A_r", "B_r", "wt", "UI"),
       n = as.integer(c(5, 5, 1, NA, 20)),
       phenotype = as.integer(c(3, 5, 15, 1, 0))
  )
community[4, ]$n <- format$n_patients - sum(community$n, na.rm = TRUE)

# rates
rates <- data.frame(
  turnover = 0.1,   # average of 10d in hosp.
  infection = 0.3   # -> R0 == 3
  )

################################################################################
# save Data                                                                    #
################################################################################
experiment_pars <- lst(
  format = format,
  antibiotics = antibiotics,
  treatments = treatments,
  community = community,
  rates = rates
  )

write(toJSON(experiment_pars, pretty = TRUE),
  "experiments/templates/default_exp_pars.json")
