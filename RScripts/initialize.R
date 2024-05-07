################################################################################
# Sript to initialize a new experiment                                         #
# usage:                                                                       #
#   Rscript.exe exp_name                                                       #
#     exp_name: name of experiment, YYYYMMDD to use system data                #
#     exp_pars: optional, path to json file with experiment parameters         #
#     gwl_pars: optional, path to json file with worklist parameters           #
################################################################################
suppressMessages(library(treatR))

# set up logging
scriptName <- "initalize"

if (Sys.info()[1] == "Windows") {
  working_directory <- "C:/Users/COMPUTER/polybox/Robot-Shared/TSPlasmids"
} else {
  working_directory <- here::here()
}
setwd(working_directory)

if (!interactive()) {
  logfile <- file(paste0("experiments/", scriptName, ".log"), open = "ab")
  sink(logfile, append = TRUE)
  sink(logfile, type = "message", append = TRUE)
}

log_string <- paste0(format(Sys.time()), " | ", scriptName, ".R ")
cat(stringr::str_pad(log_string, 79, side = "right", "-"), "\n")

################################################################################
# initialize experiment                                                        #
################################################################################
# get arguments
args <- parse_args(c("exp_name", "exp_pars", "gwl_pars"), "ccc")

if (args$exp_name == "YYYYMMDD") {
  #if the exp_name is the date, only use exp_date for name
  args$exp_name <- "output"
}

# experiment name
exp_date <- as.integer(format(Sys.Date(), "%Y%m%d"))
output_directory <- str_c(exp_date, "_", args$exp_name)

# check and load previous experiments
if (file.exists("experiments/experiments.csv")){
  old_experiments <- read_csv("experiments/experiments.csv",
    col_types = cols(
      id = col_integer(),
      created = col_datetime(format = ""),
      output_directory = col_character(),
      exp_name = col_character(),
      n_plates = col_integer()
    ))
  new_id <- as.integer(max(old_experiments$id) + 1)
} else {
  old_experiments <- tibble(
    id = integer(0),
    created = date(),
    output_directory = character(0),
    exp_name = character(0),
    n_plates = integer(0))
  new_id <- 1
}

# initialize experiment if it doesn't exist yet
if (output_directory %in% old_experiments$output_directory) {
  stop("Experiment already created! Aborting...\n\n",
    stringr::str_pad(str_c(log_string, "end "), 79,
    side = "right", "-"), "\n\n")
} else {
  #create experiment data structure
  dir.create(str_c("experiments/", output_directory, "/logs"),
    recursive = TRUE)
  dir.create(str_c("experiments/", output_directory, "/analysis"))
  dir.create(str_c("experiments/", output_directory, "/platefiles"))
  dir.create(str_c("experiments/", output_directory, "/worklists"))
  dir.create(str_c("experiments/", output_directory, "/barcode_files"))
  dir.create(str_c("experiments/", output_directory, "/OD/xml"),
    recursive = TRUE)
  dir.create(str_c("experiments/", output_directory, "/pickolo/csv"),
    recursive = TRUE)
  dir.create(str_c("experiments/", output_directory, "/pickolo/img"))
  dir.create(str_c("experiments/", output_directory, "/pickolo/json"))
  dir.create(str_c("experiments/", output_directory, "/pickolo/strainCheckPlates"))

  # setup parameters
  if (is.na(args$exp_pars)){
    # default pars
    pars <- fromJSON("experiments/templates/default_exp_pars.json")
  } else {
    # custom pars
    pars <- fromJSON(str_c("experiments/templates/", args$exp_pars))
  }
  # choose blank wells
  pars$format$rwells_bl <- list(
      sample(pars$format$n_rwells, pars$format$n_bl))
  pars$format$rwells_p <- list(
      (1:pars$format$n_rwells)[! 1:pars$format$n_rwells %in% pars$format$rwells_bl[[1]]])
  # save pars
  write(toJSON(pars, pretty = TRUE),
    str_c("experiments/", output_directory, "/exp_pars.json"))
  # experiment data
  new_experiment <- tibble(
    id = as.integer(new_id),
    created = Sys.time(),
    output_directory = output_directory,
    exp_name = args$exp_name,
    n_plates = max(pars$treatment$plate))

  experiments <- rbind(new_experiment, old_experiments)
  write_csv(experiments, file = "experiments/experiments.csv")

  #gwl pars
  if (is.na(args$gwl_pars)){
    # default pars
    file.copy(
      from = "experiments/templates/default_gwl_pars.json",
      to = str_c("experiments/", output_directory, "/gwl_pars.json"))
  } else {
    # custom pars
    file.copy(
      from = args$gwl_pars,
      to = str_c("experiments/", output_directory, "/gwl_pars.json"))
  }

  # start timekeeping
  times <- tibble(
    time = Sys.time(),
    type = "script",
    event = "initialize.R",
    run = as.integer(1),
    plate = as.integer(0),
    barcode = "none",
    replicate = as.integer(0),
    transfer = as.integer(0),
    comment = "none"
    )
  write_csv(times, file = str_c("experiments/", output_directory, "/times.csv"))

  ### LOG OUTPUT
  cat("new experiment: \n")
  print(new_experiment, n = Inf)
}

### LOG OUPTUT END
log_string <- paste0(format(Sys.time()), " | ", scriptName, ".R done ")
cat(stringr::str_pad(log_string, 79, side = "right", "-"), "\n\n")
sink()
