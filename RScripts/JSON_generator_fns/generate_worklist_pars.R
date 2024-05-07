################################################################################
#                                                                              #
################################################################################
if (Sys.info()[1] == "Windows") {
  workDir <- "D:/Polybox/robot-data/Robot-Shared/TreatmentStrategies"
} else {
  workDir <-
    "/Users/daniel/polybox/Robot-Shared/TreatmentStrategies/"
}
setwd(workDir)
library(tibble)
library(jsonlite)

pars <- fromJSON("experiments/templates/default_exp_pars.json")

vols <- tibble(
  strain = as.integer(5),
  abx = as.integer(5),
  media = as.integer(40),
  total = strain + abx + media)

lc <- tibble(
  strain_dis = "Minimal CD",
  strain_return = "Minimal CD",
  strain_asp = "Minimal CD",
  strain_asp_res = "Minimal CD trough",
  bl_dis = "Minimal FD trough",
  abx_dis = "Minimal CD",
  abx_asp = "Minimal CD")

positions <- tibble(
  strains_init_1 = "MP2pos_bck",
  strains_init_2 = "MP2pos_mid",
  strains_transfer = "MP3pos_bck",
  bl_transfer = "trough100_mid",
  abx =  "MP3pos_bck",
  source = "MP2pos_mid",
  destination = "MP3pos_mid")

strainsTransfer  <- tibble(
  rep = rep(1:8, each = length(pars$community$strain) + 1),
  strain = rep(c(pars$community$strain, "bl"), 8),
  rack = positions$strains_transfer,
  col = rep(1:12, 4),
  row = rep(c("A", "C", "E", "G"), each = 12),
  well = rrobot::rc_to_well(row, col, 8, 12),
  volume = 0)

data.frame(strainsTransfer)

strainsSetup <- tibble(
  rep = rep(1:pars$format$n_rep, each = length(pars$community$strain) + 1),
  strain = rep(c(pars$community$strain, "bl"), pars$format$n_rep),
  rack = rep(c("MP2pos_bck", "MP2pos_mid"), each = 12),
  col = rep(1:12, 2),
  volume = 0)

data.frame(strainsSetup)

abxSource <- tibble(
  drug = c("none", pars$antibiotics$drug_letter),
  rack = rep(positions$abx, length(drug)),
  col = 1:(length(drug)),
  volume = integer(1))

gwlParsJson <- lst(
  vols,
  lc,
  positions,
  strainsTransfer,
  strainsSetup,
  abxSource
  )

write(toJSON(gwl_pars_json, pretty = TRUE),
  "experiments/templates/default_gwl_pars.json")
