library(treatR)
library(grid)

experiments <- read_csv("../../experiments/experiments.csv",
  col_types = cols(
    id = col_integer(),
    created = col_datetime(format = ""),
    output_directory = col_character(),
    exp_name = col_character()))

exp_choices <- experiments$id
names(exp_choices) <- experiments$output_directory

img_opts <- list(line_size = 0.1, title_size = 5,
  label_size = 2.5, axis_text_size = 15)

colony_correct <- function(x, y, pickoloData, colfix) {
  row_c <- as.character(cut(y,
    breaks = pickoloData$pars$row_breaks,
    labels = pickoloData$pars$row_lab))
  col_c <- as.integer(as.character(cut(x,
    breaks = pickoloData$pars$col_breaks,
    labels = pickoloData$pars$col_lab)))
  np <- subset(pickoloData$growth, row == row_c & col == col_c)

  if (dim(np)[1] < 1) {
    return(colfix)
  }

  cf <- tibble(x = x, y = y, row = row_c, col = col_c,
    well = rc_to_well(row_c, col_c, 16, 24),
    auto = np$auto, manual = np$manual, change = "no",
    csv_path = pickoloData$pars$csv_path)
  previouslySelected <- which(colfix$row == row_c & colfix$col == col_c)

  if (length(previouslySelected) != 0) {
    colfix[previouslySelected, ]$manual <-
     !colfix[previouslySelected, ]$manual
    fix_manual <- FALSE
  } else {
    if (!is.na(cf$manual)) {
      cf$manual <- !cf$manual
      fix_manual <- TRUE
    } else {
      cf$manual <- !cf$auto
      fix_manual <- FALSE
    }
    colfix <- rbind(colfix, cf)
  }
  nl <- dim(colfix)[1]

  if (colfix$auto[nl] == TRUE & colfix$manual[nl] == FALSE) {
    colfix$change[nl] <- "-"
  } else if (colfix$auto[nl] == FALSE & colfix$manual[nl] == TRUE) {
    colfix$change[nl] <- "+"
  }

  if (!fix_manual) {
    colfix <- colfix[colfix$auto != colfix$manual, ]
  }

  return(colfix)
}

colfix_blank <- function() {
  df <- tibble(
    x = numeric(0),
    y = numeric(0),
    row = character(0),
    col = integer(0),
    well = integer(0),
    auto = logical(0),
    manual = logical(0),
    change = character(0),
    csv_path = character(0))
  return(df)
}

plotColonies2 <- function(
  pickoloData,
  img_path,
  img = readJPEG(img_path),
  plotLabels = FALSE,
  colChange = NULL,
  line_size = 0.1,
  title_size = 4, label_size = 0.8,
  axis_text_size = 5,
  pars = pickoloData$pars) {

  growth <- pickoloData$growth

  filename <- tail(str_split(img_path, "/")[[1]], 1)

  growth$row_int <-
    match(as.character(growth$row), LETTERS)

  select <- growth$auto
  select[!is.na(growth$manual)] <- na.omit(growth$manual)

  platePlot <- ggplot(pickoloData$data,
    aes(x, y, color = as.character(selected))) +
    annotation_custom(
      rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc")),
      0, dim(img)[2], 0, -dim(img)[1]) +
    scale_y_reverse(name = element_blank(),
      limits = c(dim(img)[1], 0),
      breaks = pars$row_breaks,
      labels = c(pars$row_lab, "")) +
    scale_x_continuous(name = element_blank(),
      limits = c(0, dim(img)[2]),
      breaks = pars$col_breaks,
      labels = c(pars$col_lab, "")) +
    geom_hline(yintercept = pars$row_wells,
      alpha = 0.5, size = line_size) +
    geom_vline(xintercept = pars$col_wells,
      alpha = 0.5, size = line_size) +
    theme_bw() +
    theme(plot.background = element_blank(),
      legend.position = "none",
      axis.text.y = element_text(vjust = 1, size = axis_text_size),
      axis.text.x = element_text(hjust = 0, size = axis_text_size),
      panel.grid.minor = element_blank()) +
    # geom_text(data = NULL, x = 800, y = 10, label = filename,
    #   color = "black", size = title_size) +
    scale_colour_manual(values = c(
      "FALSE" = rgb(1, 0, 0, 1),
      "TRUE" = rgb(0, 0, 0, 1),
      "none" = rgb(0, 0, 0, 0),
      "-" = rgb(1, 0, 0, 1),
      "+" = rgb(0, 1, 0, 1))) +
    labs(title = str_remove(filename, ".jpg"))

  if (sum(select) != 0) {
    df <- data.frame(
      x = pars$col_wells[as.numeric(growth[select, ]$col)],
      y = pars$row_wells[as.numeric(growth[select, ]$row_int)],
      xWell = pars$col_wells[2] - pars$col_wells[1],
      yWell = pars$row_wells[2] - pars$row_wells[1])

    platePlot <- platePlot +
      geom_rect(data = df,
        mapping = aes(
          xmin = x + xWell / 10,
          xmax = x + xWell - xWell / 10,
          ymin = y + yWell / 10,
          ymax = y + yWell - xWell / 10,
          color = selected),
        size = line_size * 3, color = "black", alpha = 0)
  }

  if (plotLabels){
    platePlot <- platePlot +
      geom_text(
        label = paste0(pickoloData$data$row, pickoloData$data$col),
        size = label_size)
  }

  if (!is.null(colChange)){
    if (dim(colChange)[1] > 0){
    colChange$row_int <- match(as.character(colChange$row), LETTERS)
    colChange_df <- data.frame(
      x = pars$col_wells[as.numeric(colChange$col)],
      y = pars$row_wells[as.numeric(colChange$row_int)],
      xWell = pars$col_wells[2] - pars$col_wells[1],
      yWell = pars$row_wells[2] - pars$row_wells[1],
      change = colChange$change)

    platePlot <- platePlot +
    geom_rect(data = colChange_df,
      mapping = aes(
        xmin = x + xWell / 10,
        xmax = x + xWell - xWell / 10,
        ymin = y + yWell / 10, ymax = y + yWell - xWell / 10, color = change),
      size = line_size * 3, alpha = 0)
    }
  }
  return(platePlot)
}
