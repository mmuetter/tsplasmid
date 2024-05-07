library(treatR)

setwd("/Users/daniel/polybox/Shared/Robot-Shared/treatmentstrategies/RScripts/plateViewer")

load("all_data.Rdata")

phenotype_table <- fromJSON("../../experiments/templates/phenotypes.json")


labels_phenotypes <- phenotype_table$phenotypes$let
colors_phenotypes <- phenotype_table$phenotypes$color
names(colors_phenotypes) <- phenotype_table$phenotypes$int
names(labels_phenotypes) <- phenotype_table$phenotypes$int

as.integer(cut(34, c(0:(4 + 1) * 8 + 1.5)))

all_data
stretch <- 8
plate_data <- filter(all_data, plate == 2)
onewell <- filter(plate_data, rwell == 63) %>%
  mutate(
    x = stretch * transfer_n + ceiling(rep %% 2.5),
    y = (rep < 3) * 4 + (rep > 2) * 5) %>%
  arrange(transfer_n, re)


treatment_data <- data_frame(
    x = onewell$transfer_n[c(T, F, F, F)] * stretch + 1.5,
    y = 6,
    labels = onewell$treatment_with[c(T, F, F, F)])

infecting_wells <- filter(onewell, !is.na(infection_from_well)) %>%
  select(transfer_n, infection_from_well) %>%
  mutate(transfer_n = transfer_n - 1) %>%
  distinct() %>%
  left_join(plate_data,
    by = c("transfer_n", "infection_from_well" = "rwell")) %>%
  mutate(
    x = stretch * transfer_n + ceiling(rep %% 2.5) + 2.5,
    y = (rep < 3) * 7 + (rep > 2) * 8,
    y_label = 8.8)

infected_wells <- filter(onewell, !is.na(infection_to_well)) %>%
  select(transfer_n, infection_to_well) %>%
  mutate(transfer_n = transfer_n + 1) %>%
  distinct() %>%
  left_join(plate_data,
    by = c("transfer_n", "infection_to_well" = "rwell")) %>%
  mutate(
    x = stretch * transfer_n + ceiling(rep %% 2.5),
    y = (rep < 3) * 1 + (rep > 2) * 2,
    y_label = 1.5)

turnover_data <- filter(onewell, !is.na(turnover_strain)) %>%
  select(transfer_n, turnover_strain) %>%
  distinct() %>%
  mutate(x = stretch * transfer_n + 0.5,
    y = 7.5) %>%
  left_join(phenotype_table$strains, by = c("turnover_strain" = "strain"))

transfer_data <- filter(onewell, !is.na(transfer_to_well)) %>%
  select(transfer_n) %>%
  distinct() %>%
  mutate(x = stretch * (transfer_n - 1) + 1.5,
    y = 4.5) 

selected_transfer  <-  3


chain_plot <- ggplot(onewell,
      aes(x, y, fill = as.character(phenotype), label = labels)) +
      #highlight selected transfer
      annotate("rect",
        xmin = stretch_y / 2 + selected_transfer * stretch_x - stretch_y/2,
        xmax = stretch_y / 2 + selected_transfer * stretch_x + stretch_y/2,
        ymin = -Inf, ymax = Inf,
        alpha = .2) +
      #turnover
      geom_segment(data = turnover_data,
        aes(x = x, y = y, xend = x + 1, yend = y - stretch_y),
        inherit.aes = FALSE) +
      geom_point(data = turnover_data,
        aes(x, y, color = as.character(phenotype)),
        size = 10,
        inherit.aes = FALSE) +
      geom_text(data = turnover_data, aes(x, y, label = turnover_strain),
        inherit.aes = FALSE) +
      #transfer
      geom_segment(data = transfer_data,
        aes(x = x, y = y, xend = x + stretch_x, yend = y),
        inherit.aes = FALSE) +
      #infecting wells
      geom_label(data = distinct(select(
        infecting_wells, transfer_n, y_label, infection_from_well)),
        aes(transfer_n * stretch_x + 4, y_label, label = infection_from_well),
        inherit.aes = FALSE) +
      geom_segment(data =
          group_by(infecting_wells, infection_from_well, transfer_n) %>%
          summarize(
            x_line = mean(x),
            x_line_end = x_line + stretch_x - 0.75 * stretch_y,
            y_line = mean(y),
            y_line_end = y_line - stretch_y),
        aes(x = x_line, y = y_line,
          xend = as.double(x_line_end), yend = as.double(y_line_end)),
        inherit.aes = FALSE) +
      geom_tile(data = infecting_wells,
        size = 0.25, colour = "black") +
      geom_text(data = infecting_wells) +
      #infected wells
      geom_label(data = distinct(select(
        infected_wells, transfer_n, y_label, infection_to_well)),
        aes(transfer_n * stretch_x, y_label, label = infection_to_well),
        inherit.aes = FALSE) +
      geom_segment(data =
          group_by(infected_wells, infection_to_well, transfer_n) %>%
          summarize(
            x_line = mean(x) - stretch_x,
            x_line_end = x_line + stretch_x,
            y_line_end = mean(y),
            y_line = y_line_end + stretch_y),
        aes(x = x_line, y = y_line,
          xend = as.double(x_line_end), yend = as.double(y_line_end)),
        inherit.aes = FALSE) +
      geom_tile(data = infected_wells,
        size = 0.25, colour = "black") +
      geom_text(data = infected_wells) +
      #treatment
      geom_label(data = treatment_data,
        aes(x, y, label = labels),
        inherit.aes = FALSE) +
      #main well
      geom_tile(size = 0.25, colour = "black") +
      geom_text() +
      scale_x_continuous("",
        limits =
          c(0, 2 * stretch_y / 2 + stretch_x +
            max(onewell$transfer_n) * stretch_x),
        breaks = 0:(max(onewell$transfer_n) + 1) * stretch_x + stretch_y / 2,
        labels = 0:(max(onewell$transfer_n) + 1)) +
      scale_y_continuous("",
        limits = c(0, 3 * stretch_y),
        breaks = 3 * stretch_y / 2,
        labels = unique(onewell$rwell)) +
      scale_fill_manual(values = colors_phenotypes,
        labels = labels_phenotypes,
        guide = guide_legend(title = NULL)) +
      scale_color_manual(values = colors_phenotypes,
        labels = labels_phenotypes,
        guide = guide_legend(title = NULL)) +
      theme_bw_custom_hm
    chain_plot
