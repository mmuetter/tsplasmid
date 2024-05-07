################################################################################
# Sript to generate transfer worklists                                         #
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
# args$experiment_id <- 14

# load experiment parameter
experiment <- read_csv("experiments/experiments.csv",
  col_types = cols(
    id = col_integer(),
    created = col_datetime(format = ""),
    output_directory = col_character(),
    exp_name = col_character())) %>%
  filter(id == args$experiment_id)

# set up logging
scriptName <- "mkgwl_transfer"

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

cat("creating worklists for transfer", args$transfer_n, "\n")

cat("Experiment: \n")
print.data.frame(experiment)

pars <- fromJSON(str_c("experiments/", experiment$output_directory, "/exp_pars.json"))
gwl_pars <- fromJSON(str_c("experiments/", experiment$output_directory, "/gwl_pars.json"))
plate_file <- str_c("experiments/", experiment$output_directory, "/platefiles/transfer_", args$transfer_n, ".csv")

plates <- read_csv(plate_file,
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
    treatment_with = col_character()
  )
)
################################################################################
# generate transfer worklists                                                  #
################################################################################
# 1. TURNOVER: all wells that are replaced are inoculated from source (5µl)
# 2. TRANSFER: all other wells are transfered from old plate by iPT384
# 3. INFECTION selected wells who are not turned over are transfered by
#    iPT384
# 4. TREATMENT

nMultiPipette <- 6
nAddVolumes <- 3
pt_speed_move <- 30
pt_speed_in <- 5
pt_speed_out <- 10

# 1. TURNOVER
  cat("creating plate_turnover worklists...")
  # push up all except 4 pins
  init(filename = str_c("experiments/", experiment$output_directory, "/worklists/prep4PinTool"),
    LiquidClass = gwl_pars$LC$strain_dis)
  # prepare pintool
  iPTdeactivate(filter(plates, plate == 1, rep == 1, rwell != 1)$well)
  write.gwl(gwl, quietly = TRUE)
  
  # identical plates but different source positions per plate
  # same column on source plate, for more volume
  plate <- filter(plates, plate == 1) %>%
    filter(
      !is.na(turnover_strain) &
      !(turnover_strain %in% c("bl", "UI")),
      rep == 1
    ) %>%
    arrange(turnover_strain, well)

  for (iplate in pars$treatments$plate) {
    init(filename = str_c("experiments/", experiment$output_directory, "/worklists/plate_", iplate, "_turnover"),
      LiquidClass = gwl_pars$LC$strain_dis)

    turnoverStrains <- unique(plate$turnover_strain)

    for (istrain in turnoverStrains) {
      adv_gwl_comment(str_c("Turnover. Strain ", istrain))
      istrain_source <- gwl_pars$strains_transfer384 %>%
        filter(strain == istrain, rep == 1)
      plateStrain <- plate %>%
        filter(turnover_strain == istrain)
      for (iwell in plateStrain$well) {
        iPTdip(
          rack = istrain_source$rack,
          well = istrain_source$well + 48 * (iplate - 1),
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
  }
  cat("done.\n")

# 2. TRANSFER by iPT transfer
  # identical plates (difference only in treatment)
  init(filename = str_c("experiments/", experiment$output_directory, "/worklists/plate_transfer_A"),
    LiquidClass = gwl_pars$LC$strain_dis)
  adv_gwl_comment(str_c("TRANSFER #", args$transfer_n, " prep iPT"))
  
  plate_nottransfered <- filter(plates, plate == 1, is.na(transfer_to_well))
  pins_to_deactivate <- filter(plate_nottransfered, rep == 1)$well
  
  #ready pin tool
  iPTdeactivate(pins_to_deactivate)
  write.gwl(gwl, quietly = TRUE)

  #move pin tool
  init(filename = str_c("experiments/", experiment$output_directory, "/worklists/plate_transfer_B"),
    LiquidClass = gwl_pars$LC$strain_dis)
  adv_gwl_comment(str_c("TRANSFER #", args$transfer_n, " move iPT"))
  adv_gwl_comment("iPT transfer: pick up")
  iPTdip(
    rack = "MP2pos_mid",
    well = 1,
    speed_move = pt_speed_move,
    speed_in = pt_speed_in,
    speed_out = pt_speed_out)
  adv_gwl_comment("iPT transfer: deliver")
  iPTdip(
    rack = "MP3pos_mid",
    well = 1,
    speed_move = pt_speed_move,
    speed_in = pt_speed_in,
    speed_out = pt_speed_out)
  write.gwl(gwl, quietly = TRUE)
  cat("plate transfer gwls done.\n")

# 3. INFECTION by 4pin pintool, using tips individually
  # identical plates (difference only in treatment)
  plate <- filter(plates, plate == 1)
  
  # move Pintool
  init(filename = str_c("experiments/", experiment$output_directory, "/worklists/plate_infection"),
    LiquidClass = gwl_pars$LC$strain_dis)
  #wells from - to
  rwell_to_384 <- select(plate, rep, rwell, well) %>%
    rename(infection_to_well_384 = well)
  plate_infection <- plate %>%
    filter(!is.na(infection_to_well), rep == 1) %>%
    arrange(col, row) %>%
    left_join(rwell_to_384,
      by = c("rep" = "rep", "infection_to_well" = "rwell"))
  infection_wells <- seq_along(plate_infection$well)
  for (i in infection_wells) {
    adv_gwl_comment(
      str_c("INFECTION. From ", plate_infection$well[i],
      " to ", plate_infection$infection_to_well_384[i],
      " (", i, "/", length(plate_infection$well), ")"))
    iPTdip(
      rack = gwl_pars$positions$source,
      well = plate_infection$well[i],
      speed_in = pt_speed_in,
      speed_out = pt_speed_out)
    iPTdip(
      rack = gwl_pars$positions$destination,
      well = plate_infection$infection_to_well_384[i],
      speed_in = pt_speed_in,
      speed_out = pt_speed_out)
    if (i != last(infection_wells)) {
      washPinTool()
    }
  }
  write.gwl(gwl, quietly = TRUE)
  cat("plate infection gwl done.\n")
  

# 4. TREATMENT
  # each plate is different
  cat(str_c("creating treatment worklists: plate "))
  for (iplate in pars$treatments$plate) {
    cat(str_c(iplate, " "))
    plate <- filter(plates, plate == iplate)
    init(filename = str_c("experiments/", experiment$output_directory, "/worklists/plate_", iplate, "_treatment"),
       WTtemplate = "basic", LiquidClass = gwl_pars$LC$abx_dis)

    which_treatments <- gwl_pars$abx_source$drug[gwl_pars$abx_source$drug %in% plate$treatment_with]

    for (itreatment in which_treatments) {
      adv_gwl_comment(str_c("T", args$transfer_n, " | plate ", iplate, ", Treatment: ", itreatment))
      i_abx_source <- filter(gwl_pars$abx_source, drug == itreatment)
      plate_treatment <- filter(plate, treatment_with == itreatment) %>%
        arrange(tip, col, row) %>%
        group_by(tip) %>%
        mutate(pipGroup = rep(1:ceiling(n() / nMultiPipette), each = nMultiPipette, length.out = n()))

      plate_treatment_list <- plate_treatment %>%
        group_by(pipGroup) %>%
        group_split()

      pipGroups <- seq_along(plate_treatment_list)
      for (i in pipGroups) {
        tip_volume <- plate_treatment_list[[i]] %>%
          group_by(tip) %>%
          summarise(
            n_wells = length(tip),
            volume = n_wells * gwl_pars$vols$abx + nAddVolumes * gwl_pars$vols$abx,
            .groups = "drop_last"
            )
        masks <- generateMasks(tip_volume$tip, tip_volume$volume)
        # aspirate abx
        LC_abx_asp <- gwl_pars$LC$abx_asp
        wells_aspirate <- rc_to_well(1:8, i_abx_source$col, 8, 12) * masks$tip_mask

        adv_aspirate(
          tipMask       = masks$tip_mask,
          volumes       = masks$vol_mask,
          RackLabel     = i_abx_source$rack,
          wellSelection = wells_aspirate,
          liquidClass   = LC_abx_asp,
          ncol = 12, nrow = 8)
        # return 2 volumes
        adv_dispense(
          tipMask       = masks$tip_mask,
          volumes       = 2 * gwl_pars$vols$abx,
          RackLabel     = i_abx_source$rack,
          wellSelection = wells_aspirate,
          liquidClass   = LC_abx_asp,
          ncol = 12, nrow = 8)
        # dispense in columns
        column_vector <- sort(unique(plate_treatment_list[[i]]$col))

        for (icol in column_vector) {
          for (rows_set in list(c(T, F), c(F, T))) {
            tip_volume_col <- plate_treatment_list[[i]] %>%
              filter(col == icol, row %in% LETTERS[1:16][rows_set]) %>%
            group_by(tip) %>%
            summarise(
              n_wells = length(tip),
              volume = n_wells * gwl_pars$vols$strain,
              well = well,
              .groups = "drop_last")
            masks_col <- generateMasks(tip_volume_col$tip, tip_volume_col$volume)
            adv_dispense(
              tipMask       = masks_col$tip_mask,
              volumes       = masks_col$vol_mask,
              RackLabel     = gwl_pars$positions$destination,
              wellSelection = tip_volume_col$well,
              ncol = 24, nrow = 16)
          }
        }
        if (i != last(pipGroups)) {
          #only wash in between not at end (end -> washing in evoware)
          sterile_wash_c()
        }
      }
      # clean tips between treatments, but not at end (to allow for
      # simultanous plate movement while cleaning at the end)
      if (itreatment != last(which_treatments)) {
        sterile_wash_c()
      }
    }
    write.gwl(gwl, quietly = TRUE)
  }
cat("\n")
cat("treatment gwls done.\n")

### LOG OUTPUT END
log_string <- paste0(format(Sys.time()), " | ", scriptName, ".R done.")
cat(stringr::str_pad(log_string, 79, side = "right", "-"), "\n\n")
sink()
