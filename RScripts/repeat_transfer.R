################################################################################
# Sript to repeat transfers                                                    #
# usage:                                                                       #
#   Rscript.exe [params] experiment_id transfer_n                              #
################################################################################
suppressMessages(library(treatR))

if (Sys.info()[1] == "Windows") {
  working_directory <- "C:/Users/COMPUTER/polybox/Robot-Shared/TSPlasmids"
} else {
  working_directory <- here::here()
}
setwd(working_directory)

# get arguments
args <- parse_args(c("experiment_id", "transfer_n"), "ii")

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
scriptName <- "repeat_transfer"

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
print.data.frame(experiment)

pars <- fromJSON(str_c("experiments/", experiment$output_directory,
  "/exp_pars.json"))
transfer_file <- str_c("experiments/", experiment$output_directory,
  "/transfers.csv")

# load transfer file
transfers <- read_csv(transfer_file,
  col_types = cols(
    n = col_integer(),
    run = col_integer(),
    created = col_datetime(format = ""),
    type = col_character()
    ))

next_transfer <- filter(transfers, n == args$transfer_n) %>%
  filter(run == max(run)) %>%
  mutate(
    created = Sys.time(),
    run = as.integer(run + 1))

transfers <- rbind(next_transfer, transfers)
new_run <- max(transfers$run)
################################################################################
# move files from later transfers                                              #
################################################################################
later_transfers <- filter(transfers, run < max(run), n >= args$transfer_n)
# OD
old_path <- str_c("experiments/", experiment$output_directory, "/OD/xml")
new_path <- str_c("experiments/", experiment$output_directory,
  "/OD/xml/run", new_run - 1)
dir.create(new_path)

od_files <- list.files(path = old_path,
  pattern = str_c("t[", str_c(later_transfers$n, collapse = ","),"]"))
od_files <- c(od_files, # OD after infection for plates that are reused
  list.files(path = old_path,
    pattern = str_c("t",args$transfer_n-1,".*AI")))

file.rename(
  from = str_c(old_path, "/", od_files),
  to = str_c(new_path, "/", od_files))

#pickolo
old_path <- str_c("experiments/", experiment$output_directory, "/pickolo")
new_path <- str_c("experiments/", experiment$output_directory,
  "/pickolo/run", new_run - 1)
dir.create(new_path)
for(i in c("csv","img","json")) {
  dir.create(str_c(new_path, "/", i))
}

pickolo_files <- list.files(
  path = old_path,
  pattern = str_c("t[",
    str_c(c(args$transfer_n-1, later_transfers$n), collapse = ","),"]"),
  recursive = TRUE)

file.rename(
  from = str_c(old_path, "/", pickolo_files),
  to = str_c(new_path, "/", pickolo_files))

# platefiles
old_path <- str_c("experiments/", experiment$output_directory, "/platefiles")
new_path <- str_c("experiments/", experiment$output_directory,
  "/platefiles/run", new_run - 1)
dir.create(new_path)

plate_files <- list.files(
  path = old_path,
  pattern = str_c("_[", str_c(later_transfers$n+1, collapse = ","),"]"))

file.rename(
  from = str_c(old_path, "/", plate_files),
  to = str_c(new_path, "/", plate_files))

################################################################################
# save                                                                         #
################################################################################
#save transfers
write_csv(transfers, file = transfer_file)
repeat_state <- data_frame(repeatState = TRUE)
write_csv(repeat_state, file = str_c("experiments/", experiment$output_directory,
  "/repeat_state.csv"))
### LOG OUTPUT
cat("repeating transfer:\n")
print.data.frame(next_transfer)

### LOG OUTPUT END
log_string <- paste0(format(Sys.time()), " | ", scriptName, ".R done.")
cat(stringr::str_pad(log_string, 79, side = "right", "-"), "\n\n")
sink()
