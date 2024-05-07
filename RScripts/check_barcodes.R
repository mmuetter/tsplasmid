################################################################################
# Sript to check for obvious errors in barcodes in StoreX                      #
# usage:                                                                       #
#   Rscript.exe [params] experiment_id transfer run_n                          #
################################################################################
suppressMessages(library(treatR))


if (Sys.info()[1] == "Windows") {
  working_directory <- "C:/Users/COMPUTER/polybox/Robot-Shared/TSPlasmids"
  storex_file <- "C:/ProgramData/Tecan/EVOware/output/StoreX_PosList.txt"
} else {
  working_directory <- here::here()
  storex_file <- here::here("RScripts/StoreX_PosListExample.txt")
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
    exp_name = col_character())) %>%
  filter(id == args$experiment_id)

# set up logging
scriptName <- "check_barcodes"

if (!interactive()) {
  logfile <- file(str_c("experiments/", experiment$output_directory, "/logs/", scriptName, "_t", args$transfer_n ,".log"), open = "ab")
  sink(logfile, append = TRUE)
  sink(logfile, type = "message", append = TRUE)
}

log_string <- paste0(format(Sys.time()), " | ", scriptName, ".R ")
cat(stringr::str_pad(log_string, 79, side = "right", "-"), "\n")

################################################################################
# Load                                                                         #
################################################################################

cat("Experiment:\n")
print.data.frame(experiment)
cat("\n")

read_storex <- function(storex_file) {
    storex_raw <- read.table(file = storex_file, stringsAsFactors = FALSE)
    storex <- as_tibble(str_split(storex_raw[, 1], ",|=", simplify = TRUE)) %>%
        mutate(
          cartridge = suppressWarnings(as.integer(V1)),
          position = as.integer(V2),
          barcode = as.character(V3)) %>%
        filter(!is.na(cartridge)) %>%
        select(cartridge, position, barcode)
    return(storex)
}

storex_bc <- read_storex(storex_file)
plates_bc <- read_csv(
  str_c("experiments/", experiment$output_directory, "/barcodes.csv"),
  col_types = cols(
    transfer = col_integer(),
    run = col_integer(),
    desc = col_character(),
    plate_n = col_integer(),
    barcode = col_character(),
    ab = col_character()))

# output
output <- list()
# have all barcodes been read properly?
output$readOK <- all(!str_detect(storex_bc$barcode, "Gen"))
# are all barcodes that are needed present?
assayplatesNew <- filter(plates_bc,
  transfer == args$transfer_n, run == args$run_n, desc == "assay")$barcode

if (args$transfer_n > 0){
  on <- filter(plates_bc,
    transfer == args$transfer_n, run == args$run_n, desc == "on")$barcode
  newAgar <- filter(plates_bc,
    transfer == args$transfer_n - 1, run == max(run), desc == "agar")$barcode
  assayplatesOld <- filter(plates_bc,
    transfer == args$transfer_n - 1, run == args$run_n, desc == "assay")$barcode
  checkPlates <- filter(plates_bc,
    transfer == args$transfer_n, run == args$run_n, desc == "test")$barcode
} else {
  on <- NULL
  newAgar <- NULL
  assayplatesOld <- NULL
  checkPlates <- NULL
}

if (args$transfer_n > 1){
  oldAgar <- filter(plates_bc,
    transfer == args$transfer_n - 2, run == max(run), desc == "agar")$barcode 
  oldCheckPlates <- filter(plates_bc,
    transfer == args$transfer_n - 1, run == args$run_n, desc == "test")$barcode
} else {
  oldAgar <- NULL
  oldCheckPlates <- NULL
}

neededBarcodes <- c(assayplatesNew, assayplatesOld, on, newAgar, oldAgar, checkPlates, oldCheckPlates)
output$missingBC <- filter(plates_bc, barcode %in%
    neededBarcodes[!as.character(neededBarcodes) %in% storex_bc$barcode])

output$duplicateBarcodes <- filter(storex_bc, barcode != "") %>%
  filter(duplicated(barcode) | duplicated(barcode, fromLast = TRUE))

cat("All barcodes read ok? ", output$readOK, "\n")

if (dim(output$duplicateBarcodes)[1] > 0){
  cat("DUPLICATE BARCODES!!\n")
  print.data.frame(output$duplicateBarcodes)
  cat("\n")
}

if (dim(output$missingBC)[1] > 0){
  cat("Some missing barcodes:\n")
  print.data.frame(output$missingBC)
  cat("\n")
} else {
  cat("no missing barcodes, as far as I know!\n")
}

### LOG OUTPUT END
log_string <- paste0(format(Sys.time()), " | ", scriptName, ".R done.")
cat(stringr::str_pad(log_string, 79, side = "right", "-"), "\n\n")
sink()
