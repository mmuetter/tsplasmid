################################################################################
# Script to make barcode files for transfer for reading in EvoWare             #
#                                                                              #
# usage:                                                                       #
#   Rscript.exe [params] experiment_id transfer_n run_n                        #
################################################################################
suppressMessages(library(treatR))

if (Sys.info()[1] == "Windows") {
  working_directory <- "C:/Users/COMPUTER/polybox/Robot-Shared/TSPlasmids"
} else {
  working_directory <- here::here()
}
setwd(working_directory)

# get arguments
args <- parse_args(c("experiment_id", "transfer_n", "run_n"), "iii")

# load experiment parameter
experiment <- read_csv("experiments/experiments.csv",
  col_types = cols(
    id = col_integer(),
    created = col_datetime(format = ""),
    output_directory = col_character(),
    exp_name = col_character(),
    n_plates = col_integer())) %>%
  filter(id == args$experiment_id)

# set up logging
scriptName <- "generate_bc_files_endread"

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

# load barcodes
barcodes <- read_csv(str_c("experiments/", experiment$output_directory,
  "/barcodes.csv"),
  col_types = cols(
    transfer = col_integer(),
    run = col_integer(),
    desc = col_character(),
    plate_n = col_integer(),
    barcode = col_character(),
    ab = col_character()))

barcode_file_path <- str_c("experiments/", experiment$output_directory,
  "/barcode_files")

# generate barcode files for agar plates to be read
if (args$transfer_n > 1){
  # source plates
  gen_barcodefile(barcodes, barcode_file_path,
    "bc_last_agar_source", args$transfer_n, "max", 1, "assay")
  # agar plates
  for (iplate in 1:experiment$n_plates){
    gen_barcodefile(barcodes, barcode_file_path,
      str_c("bc_last_agar_", iplate), args$transfer_n, "max", iplate, "agar")
  }
}

### LOG OUTPUT END
log_string <- paste0(format(Sys.time()), " | ", scriptName, ".R done.")
cat(stringr::str_pad(log_string, 79, side = "right", "-"), "\n\n")
sink()
