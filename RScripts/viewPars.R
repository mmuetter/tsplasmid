suppressMessages(library(treatR))

working_directory <- here::here()
setwd(working_directory)

experiment_id <- 2

experiment <- read_csv("experiments/experiments.csv",
  col_types = cols(
    id = col_integer(),
    created = col_datetime(format = ""),
    output_directory = col_character(),
    exp_name = col_character())) %>%
  filter(id == experiment_id)
cat("Experiment: \n")
print(experiment)

pars <- fromJSON(str_c("experiments/", experiment$output_directory,
  "/exp_pars.json"))
print(pars)
gwl_pars <- fromJSON(str_c("experiments/", experiment$output_directory,
  "/gwl_pars.json"))
print(gwl_pars)

plate <- data_frame(
    well = 1:96,
    row = rep(LETTERS[1:8], each = 12),
    col = rep(1:12,8)) %>%
  left_join(data_frame(row = LETTERS[1:8], row_int = 8:1), by = "row") %>%
  left_join(
    select(gwl_pars$strains_transfer, strain, well), by = "well")

ggplot(plate,
  aes(col, row_int, color = strain, label = strain)) +
  geom_point(size = 4) +
  geom_label() +
  scale_y_continuous("row", breaks = 8:1, labels = LETTERS[1:8]) +
  scale_x_continuous("column", breaks = 1:12)
