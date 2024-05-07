################################################################################
# Sript to record time & Barcode                                               #
# usage:                                                                       #
#   Rscript.exe experiment_id type event run plate barcode replicate           #
#     transfer comment                                                         #
################################################################################
suppressMessages(library(treatR))

if (Sys.info()[1] == "Windows") {
  working_directory <- "C:/Users/COMPUTER/polybox/Robot-Shared/TSPlasmids"
} else {
  working_directory <- here::here()
}
setwd(working_directory)

# get arguments
args <- parse_args(
  c("experiment_id", "type", "event", "run", "plate", "barcode", "replicate",
  "transfer", "comment"), "icciiciic")

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
scriptName <- "record_time"

if (!interactive()) {
  logfile <- file(str_c("experiments/", experiment$output_directory, "/logs/", scriptName, "_t", args$transfer, ".log"), open = "ab")
  sink(logfile, append = TRUE)
  sink(logfile, type = "message", append = TRUE)
}

log_string <- paste0(format(Sys.time()), " | ", scriptName, ".R ")
cat(stringr::str_pad(log_string, 79, side = "right", "-"), "\n")

################################################################################
# Load                                                                         #
################################################################################

# new timestamp
new_time <- mutate(args, time = Sys.time()) %>%
  select(time, type:comment)

# old timestamps
times_path <- str_c("experiments/", experiment$output_directory, "/times.csv")

if (file.exists(times_path)){
  times <- read_csv(times_path,
    col_types = cols(
      time = col_datetime(format = ""),
      type = col_character(),
      event = col_character(),
      run = col_integer(),
      plate = col_integer(),
      barcode = col_character(),
      replicate = col_integer(),
      transfer = col_integer(),
      comment = col_character()))
} else {
  times <- data_frame(
    time = Sys.time(),
    type = "newTimeFile",
    event = "recordTime.R",
    run = as.integer(1),
    plate = as.integer(0),
    barcode = "none",
    replicate = as.integer(0),
    transfer = as.integer(0),
    comment = "none")
  cat("New times.csv file created.\n")
}


# combine
times_out <- rbind(new_time, times)

# save
write_csv(times_out, file = times_path)

### LOG OUTPUT
cat("Timestamp:", format(new_time$time), "added.\n")

### LOG OUTPUT END
log_string <- paste0(format(Sys.time()), " | ", scriptName, ".R done.")
cat(stringr::str_pad(log_string, 79, side = "right", "-"), "\n\n")
sink()
