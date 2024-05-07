################################################################################
# Sript to generate setup worklists                                            #
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
    exp_name = col_character())) %>%
  filter(id == args$experiment_id)

# set up logging
scriptName <- "mkgwl_setup"

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
gwl_pars <- fromJSON(str_c("experiments/", experiment$output_directory,
  "/gwl_pars.json"))

plates <- read_csv(
  str_c("experiments/", experiment$output_directory,
    "/platefiles/transfer_0.csv"),
  col_types = cols(
    plate = col_integer(),
    well = col_integer(),
    row = col_character(),
    col = col_integer(),
    tip = col_integer(),
    rep = col_integer(),
    rwell = col_integer(),
    transfer_n = col_integer(),
    turnover_strain = col_character(),
    transfer_to_well = col_integer(),
    infection_to_well = col_integer(),
    treatment_with = col_character())
)
################################################################################
# generate setup worklists                                                     #
################################################################################

nMultiPipette <- 6
nAddVolumes <- 4
pt_speed_move <- 30
pt_speed_in <- 5
pt_speed_out <- 10


cat("creating plate_setup worklists: ")

# push up all except 4 pins
  init(filename = str_c("experiments/", experiment$output_directory, "/worklists/prep4PinTool"),
    LiquidClass = gwl_pars$LC$strain_dis)
  # prepare pintool
  iPTdeactivate(filter(plates, plate == 1, rep == 1, rwell != 1)$well)
  write.gwl(gwl, quietly = TRUE)

# identical plates (difference only in treatment)
plate <- filter(plates, plate == 1) %>%
  filter(
    !is.na(turnover_strain) &
    !(turnover_strain %in% c("bl", "UI")),
    rep == 1
  ) %>%
  arrange(turnover_strain, well)

init(filename = str_c("experiments/", experiment$output_directory, "/worklists/plate_setup"),
  LiquidClass = gwl_pars$LC$strain_dis)

turnoverStrains <- unique(plate$turnover_strain)

for (istrain in turnoverStrains) {
  adv_gwl_comment(str_c("Setup. Strain ", istrain))
  istrain_source <- gwl_pars$strains_transfer384 %>%
    filter(strain == istrain, rep == 1)
  plateStrain <- plate %>%
    filter(turnover_strain == istrain)
  for (iwell in plateStrain$well){
    iPTdip(
      rack = istrain_source$rack,
      well = istrain_source$well,
      speed_in = pt_speed_in,
      speed_out = pt_speed_out)
    iPTdip(
      rack = gwl_pars$positions$destination,
      well = iwell,
      speed_in = pt_speed_in,
      speed_out = pt_speed_out)
  }
  if (istrain != last(turnoverStrains)) {
      washPinTool()
  }
}
write.gwl(gwl, quietly = TRUE)
cat("done.\n")

### LOG OUTPUT END
log_string <- paste0(format(Sys.time()), " | ", scriptName, ".R done ")
cat(stringr::str_pad(log_string, 79, side = "right", "-"), "\n\n")
sink()
