################################################################################
# Sript to generate transfers                                                  #
# usage:                                                                       #
#   Rscript.exe [params] experiment_id                                         #
################################################################################
suppressMessages(library(treatR))

if (Sys.info()[1] == "Windows") {
  working_directory <- "C:/Users/COMPUTER/polybox/Robot-Shared/TSPlasmids"
} else {
  working_directory <- here::here()
}
setwd(working_directory)

# get arguments
args <- parse_args(c("experiment_id"), "i")

# load experiment parameter
experiment <- read_csv("experiments/experiments.csv",
  col_types = cols(
    id = col_integer(),
    created = col_datetime(format = ""),
    output_directory = col_character(),
    exp_name = col_character(),
    n_plates = col_integer()
  )) %>%
  filter(id == args$experiment_id)

# set up logging
scriptName <- "new_transfer"

if (!interactive()) {
  logfile <- file(str_c("experiments/", experiment$output_directory, "/logs/", scriptName, ".log"), open = "ab")
  sink(logfile, append = TRUE)
  sink(logfile, type = "message", append = TRUE)
}

log_string <- paste0(format(Sys.time()), " | ", scriptName, ".R ")
cat(stringr::str_pad(log_string, 79, side = "right", "-"), "\n")

################################################################################
# Load                                                                         #
################################################################################
cat("Experiment: \n")
print(experiment)

pars <- fromJSON(str_c("experiments/", experiment$output_directory,
  "/exp_pars.json"))
transfer_file <- str_c("experiments/", experiment$output_directory,
  "/transfers.csv")

# load or create transfer file
if (file.exists(transfer_file)){
  old_transfers <- read_csv(transfer_file,
    col_types = cols(
      n = col_integer(),
      run = col_integer(),
      created = col_datetime(format = ""),
      type = col_character()))
  transfer_n  <- old_transfers$n[1] + 1
  run <- old_transfers$run[1]
  transfer_type  <- "transfer"
} else {
  old_transfers <- tibble(
      n = integer(0),
      run = integer(0),
      created = character(0),
      type = character(0))
  transfer_n  <- 0
  run <- 1
  transfer_type  <- "setup"
}
# add new transfer
new_transfer <- tibble(
  n = as.integer(transfer_n),
  run = as.integer(run),
  created = Sys.time(),
  type = transfer_type)

transfers <- rbind(new_transfer, old_transfers)

################################################################################
# setup plates                                                                 #
################################################################################
plate <- read_csv("experiments/templates/plate_template.csv",
  col_types = cols(
    well = col_integer(),
    row = col_character(),
    col = col_integer(),
    tip = col_integer(),
    rep = col_integer(),
    rwell = col_integer()
))

plates <- rep(lst(plate), length(pars$treatments$plate))
names(plates) <- pars$treatments$plate
plates <- bind_rows(plates, .id = "plate")
plates$plate <- as.integer(plates$plate)

################################################################################
# assign treatment, turnover, transfer and infection                           #
################################################################################
# turnover
forceTurnoverFile <- file.path("experiments", experiment$output_directory, str_c("t", transfer_n, "forceTurnover.txt"))
if (file.exists(forceTurnoverFile)) {
  forceTurnoverWells <- scan(forceTurnoverFile, quiet = TRUE)
} else {
  forceTurnoverWells <- NULL
}

new_turnover_strain <- gen_turnover(pars, transfer = transfer_n, forceTurnover = forceTurnoverWells)
cat("turnover:")
table(new_turnover_strain$turnover_strain)
if (!is.null(forceTurnoverWells)) {
  cat("!! Forced turnover for wells: \n")
  cat("  ", forceTurnoverWells)
}
cat("\n")
# transfer every well that is not turned over
new_transfer_to_well <- tibble(
      rwell = 1:pars$format$n_rwells,
      transfer_to_well = rwell)
new_transfer_to_well$transfer_to_well[
  !is.na(new_turnover_strain$turnover_strain)] <- NA
# infection
new_infection_to_well <- gen_infection(pars, transfer = transfer_n, ignoreRWell = forceTurnoverWells)

cat("Infection: (", sum(!is.na(new_infection_to_well$infection_to_well)), ")")
table(new_infection_to_well$infection_to_well)
cat("\n")
# treatment
 plate_treatments <- pars$treatments$plate
new_treatment_with <- map(
  plate_treatments,
  ~gen_treatment(pars = pars, condition = .x, transfer = transfer_n)) %>%
  bind_rows(.id = "plate") %>%
  mutate(plate = as.integer(plate))
cat("Treatment: \n")
for (i in pars$treatments$plate){
  cat(" - Plate ", i)
  print(
    table(new_treatment_with$treatment_with[new_treatment_with$plate == i]))
}
cat("\n")
#add to plate file
plates <- plates %>%
  # transfer number
  mutate(transfer_n = as.integer(transfer_n)) %>%
  # turnover
  left_join(new_turnover_strain, by = "rwell") %>%
  # transfer
  left_join(new_transfer_to_well, by = "rwell") %>%
  # infection
  left_join(new_infection_to_well, by = "rwell") %>%
  # treatment
  left_join(new_treatment_with, by = c("plate", "rwell"))

################################################################################
# save                                                                         #
################################################################################
#save transfers
write_csv(transfers, file = transfer_file)
#save plate file
plates_path <- str_c("experiments/", experiment$output_directory,
  "/platefiles/transfer_", transfer_n, ".csv")
write_csv(plates, file = plates_path)
repeat_state <- tibble(repeatState = FALSE)
write_csv(repeat_state, file = str_c("experiments/", experiment$output_directory, "/repeat_state.csv"))

### LOG OUTPUT
cat("new transfer: \n")
print(new_transfer)

cat("new plate file saved to: ", plates_path, "\n")

### LOG OUTPUT END
log_string <- paste0(format(Sys.time()), " | ", scriptName, ".R done ")
cat(stringr::str_pad(log_string, 79, side = "right", "-"), "\n\n")
sink()
