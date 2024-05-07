if (Sys.info()[1] == "Windows"){
  working_directory <- "D:/Polybox/robot-data/Robot-Shared/TreatmentStrategies"
} else {
  working_directory <-
    "/Users/Daniel/polybox/Shared/Robot-Shared/TreatmentStrategies"
}
setwd(working_directory)
library(tidyverse)

plate <- data_frame(
    well = 1:384,
    row = rep(LETTERS[1:16], each = 24),
    col = rep(1:24, 16),
    tip = rep(1:8, each = 24 * 2),
    rep = as.integer(rep(c(rep(c(1,2),12),rep(c(3,4),12)),8))) %>%
  group_by(rep) %>%
  mutate(rwell = 1:96)

write_csv(plate, file = "experiments/templates/plate_template.csv")
