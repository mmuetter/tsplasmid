library(treatR)

experiments <- read_csv("../../experiments/experiments.csv",
  col_types = cols(
    id = col_integer(),
    created = col_datetime(format = ""),
    output_directory = col_character(),
    exp_name = col_character()))

phenotype_table <- fromJSON("../../experiments/templates/phenotypes.json")

exp_choices <- experiments$id
names(exp_choices) <- experiments$output_directory

theme_bw_custom_hm <- theme_bw() +
  theme(
    # overall text size
    text = element_text(size = 12),
    # legend
    legend.position = "right",
    legend.key = element_rect(color = NA),
    # facet labels
    strip.background = element_rect(fill = NA, color = NA),
    strip.text.x = element_text(size = 12, hjust = 0),
    # panel
    panel.border = element_rect(color = NA),
    panel.grid = element_blank(),
    # axis
    axis.line.x  = element_blank(),
    axis.line.y  = element_blank(),
    axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
    axis.title.y = element_text(margin = margin(0, 10, 0, 0))
  )

#color
library(RColorBrewer)
hmPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space = "Lab")

labels_phenotypes <- phenotype_table$phenotypes$let
colors_phenotypes <- phenotype_table$phenotypes$color
names(colors_phenotypes) <- phenotype_table$phenotypes$int
names(labels_phenotypes) <- phenotype_table$phenotypes$int

labels_changes <- phenotype_table$changes_legend$class_long
colors_changes <- phenotype_table$changes_legend$color
names(colors_changes) <- phenotype_table$changes_legend$class_short
names(labels_changes) <- phenotype_table$changes_legend$class_short

#chain_plot
stretch_x <- 8
stretch_y <- 3

base_coords <- data_frame(
      rep = 1:4,
      x_base = c(1, 2, 2, 1),
      y_base = c(1, 1, 0, 0))

y_heights <- (1/3 + c(0, 1, 2)) * stretch_y
label_height <- 1.8
box_size <- stretch_y
