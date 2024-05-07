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
antibiotics <- tibble(
  drug_letter = c("AB22", "AB12", "AB21", "AB11", "AB051", "AB105", "AB0505"),
  drug_short = c("NalStr22", "NalStr12", "NalStr21", "NalStr11", "NalStr051", "NalStr105", "NalStr0505"),
  drug_long = c("Nal + Sm", "Nal + Sm", "Nal + Sm", "Nal + Sm", "Nal + Sm", "Nal + Sm", "Nal + Sm"),
  MIC_wt = c(0, 0, 0, 0, 0, 0, 0),
  C_treat = c(0, 0, 0, 0, 0, 0, 0))

format <- data_frame(
  format = 384,
  n_rep = 4,
  n_rwells = format / n_rep,
  n_bl = 2,
  n_patients = n_rwells - n_bl,
  rwells_bl = list(1:n_bl),
  rwells_p = list((n_bl+1):n_rwells))

treatments <- tibble(
  plate = 1:8,
  conditions_name = c(
    "no treatment",
    "combo: 2/2",
    "combo: 1/2",
    "combo: 2/1",
    "combo: 1/1",
    "combo: 0.5/1",
    "combo: 1/0.5",
    "combo: 0.5/0.5"
    ),
  drugs = c(
    "none",
    "AB22",
    "AB12",
    "AB21",
    "AB11",
    "AB051",
    "AB105",
    "AB0505"),
  concentration = c( # not used (yet?)
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1),
  period = c(1, 1, 1, 1, 1, 1, 1, 1), # number of days on one drug
  probability = c(1, 1, 1, 1, 1, 1, 1, 1), # probability of receiving A/B
  description = c(
    "no treatment",
    "Combination therapy\nNal20 + Sm12.5",
    "Combination therapy\nNal10 + Sm12.5",
    "Combination therapy\nNal20 + Sm6.25",
    "Combination therapy\nNal10 + Sm6.25",
    "Combination therapy\nNal5 + Sm6.25",
    "Combination therapy\nNal10 + Sm3.125",
    "Combination therapy\nNal5 + Sm3.125")
  )

# frequency of the strains in the environment, adds up to n_patients
community <- data_frame(
  strain = c("AB_r","A_r", "B_r", "wt", "UI"),
       n = as.integer(c(0, 0, 0, 80, 14)),
       phenotype = as.integer(c(3, 5, 15, 1, 0))
  )


# rates
rates <- data.frame(
  turnover = 0.2,   # average of 5d in hosp.
  infection = 0.3   # -> R0 == 1.5
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
