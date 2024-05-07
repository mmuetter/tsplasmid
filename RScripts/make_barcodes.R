################################################################################
# Sript make new barcodes                                                      #
# usage:                                                                       #
#   Rscript.exe experiment_id transfer_n run_n (basedate)                      #
################################################################################
suppressMessages(library(treatR))

if (Sys.info()[1] == "Windows") {
  working_directory <- "C:/Users/COMPUTER/polybox/Robot-Shared/TSPlasmids"
} else {
  working_directory <- here::here()
}
setwd(working_directory)

# get arguments
args <- parse_args(c("experiment_id", "transfer_n", "run_n", "basedate"),
  "iiii")
# Manual mode
  # args$experiment_id <- as.integer(14)
  # args$transfer_n <- as.integer(0)
  # args$run_n <- as.integer(1)
  # args$basedate <- 20220127

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
scriptName <- "make_barcodes"

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

if (is.na(args$basedate)) {
  args$basedate <- as.integer(format(Sys.Date(), "%Y%m%d"))
}

pars <- fromJSON(str_c("experiments/", experiment$output_directory,
  "/exp_pars.json"))
barcode_file <- str_c("experiments/", experiment$output_directory,
  "/barcodes.csv")

# load/make barcode file
if (file.exists(barcode_file)) {
  barcodes <- read_csv(file = barcode_file,
    col_types = cols(
      transfer = col_integer(),
      run = col_integer(),
      desc = col_character(),
      plate_n = col_integer(),
      barcode = col_integer(),
      ab = col_character())
    )
} else {
  barcodes <- tibble(
    transfer = integer(0),
    run = integer(0),
    desc = character(0),
    plate_n = integer(0),
    abx = character(0),
    barcode = integer(0))
}

barcodes_present <- dim(filter(barcodes,
  transfer == args$transfer_n,
  run == args$run_n,
  desc == "assay"))[1] > 0

if (barcodes_present){
  warning(str_c("Barcodes for transfer ",
    args$transfer_n, ", run ", args$run_n, " allready present!"))
} else {
  n_plates <- max(pars$treatments$plate)

  new_barcodes <- tibble(
    transfer = as.integer(
        c(
          rep(c(args$transfer_n, args$transfer_n + 1), c(n_plates * 5, 1)),
          rep(args$transfer_n, 4)
        )
      ),
    run = args$run_n,
    desc = rep(
      c("assay", "agar", "on", "test"),
      c(n_plates, n_plates * 4, 1, 4)),
    plate_n = as.integer(c(rep(1:n_plates, 5), 1, rep(1, 4))),
    barcode = as.integer(args$basedate * 100 + 1:(n_plates * 5 + 1 + 4)),
    ab = c(
        rep(
          c("N", "N", "A", "B", "AB", "N"),
          c(rep(n_plates, 5), 1)
        ),
        c("N", "A", "B", "AB")
      )
  )

  if (args$transfer_n == 0) {
    new_barcodes <- new_barcodes %>%
      filter(plate_n == 1)
  }

  cat("generated new barcodes",
    first(new_barcodes$barcode, 1), "to", last(new_barcodes$barcode),
    "\n")

  if(args$run_n > 1){
    #check if barcodes for agar plates to spot are present
    spot_barcodes_present <- dim(filter(barcodes,
      transfer == args$transfer_n-1,
      run == args$run_n,
      desc == "agar"))[1] > 0
    if(!spot_barcodes_present){
      new_barcodes2 <- tibble(
        transfer = as.integer(c(
          rep(args$transfer_n-1, n_plates * 4),args$transfer_n)),
        run = args$run_n,
        desc = rep(
          c("agar", "on"),
          c(n_plates * 4, 1)),
        plate_n = as.integer(c(rep(1:n_plates, 4), 1)),
        barcode = as.integer(args$basedate * 100 + 32:(32 + n_plates * 4)),
        ab = rep(
          c("N", "A", "B", "AB", "N"),
          c(rep(n_plates, 4), 1)))
      cat("generated new barcodes",
        first(new_barcodes2$barcode, 1), "to", last(new_barcodes2$barcode),
        "\n")
      new_barcodes <- rbind(new_barcodes, new_barcodes2)
    }
  }
  barcodes <- rbind(barcodes, new_barcodes)
  write_csv(barcodes, file = barcode_file)
}

### LOG OUTPUT END
log_string <- paste0(format(Sys.time()), " | ", scriptName, ".R done.")
cat(stringr::str_pad(log_string, 79, side = "right", "-"), "\n\n")
sink()
