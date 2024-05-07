shinyServer(function(input, output, session){
  ##############################################################################
  # Observers                                                                  #
  ##############################################################################
  observe({
    updateNumericInput(session, "plate", max = exp()$n_plates)
    updateNumericInput(session, "transfer", max = max(plates_data()$transfer_n))
  })

  ##############################################################################
  # paths                                                                      #
  ##############################################################################
  exp <- reactive({
    exp <- get_experiment(exp_id = input$exp, dir = "../../experiments/")
    print(exp)
    return(exp)
  })

  od_path <- reactive({
    od_path <- str_c("../../experiments/", exp()$output_directory,
      "/od/xml")
    return(od_path)
    })

  pickolo_path <- reactive({
    pickolo_path <- str_c("../../experiments/", exp()$output_directory,
      "/pickolo/json")
    return(pickolo_path)
  })

  plates_path <- reactive({
    plates_path <- str_c("../../experiments/", exp()$output_directory,
      "/platefiles")
    return(plates_path)
  })


  ##############################################################################
  # reactive Data & Observers                                                  #
  ##############################################################################
  values <- reactiveValues(well = data_frame(
      col = 0.5,
      row = 0.5,
      rwell = 1))

  observe({
    if (!is.null(input$click)) {
      well <- nearPoints(all_data(), input$click,
        threshold = 50, maxpoints = 1) %>%
        select(col, row_int, rwell)
      if (well$col %% 2 == 0){
        well$col <- well$col - 1
      }
      if (well$row_int %% 2 == 0){
        well$row_int <- well$row_int - 1
      }
      values$well <- data_frame(
        col = well$col - 0.5,
        row = well$row_int - 0.5,
        rwell = well$rwell)
    }
  })


  ##############################################################################
  # data                                                                       #
  ##############################################################################
  exp_data <- reactive({
    exp_data <- get_experiment_data(exp())
  })


  od_data <- reactive({
    return(exp_data()$od_data)
  })

  pickolo_data <- reactive({
    return(exp_data()$pickolo_data)
  })

  plates_data <- reactive({
    return(exp_data()$plates_data)
  })

  all_data <- reactive({
    all_data <- exp_data()$all_data
    return(all_data)
  })

  plate_data <- reactive({
    plate_data <- filter(all_data(), plate == input$plate) %>%
      left_join(base_coords, by = "rep")
    return(plate_data)
  })

  onewell <- reactive({
    onewell <- filter(plate_data(), rwell == as.integer(values$well[1, 3])) %>%
      mutate(
        x = x_base + stretch_x * transfer_n,
        y = y_base + y_heights[2])
    return(onewell)
  })

  treatment_data <- reactive({
    treatment_data <- data_frame(
      x = onewell()$transfer_n[c(T, F, F, F)] * stretch_x + box_size / 2,
      y = y_heights[2] + label_height,
      labels = onewell()$treatment_with[c(T, F, F, F)])
    return(treatment_data)
  })

  infecting_wells <- reactive({
    infecting_wells <- filter(onewell(), !is.na(infection_from_well)) %>%
      select(transfer_n, infection_from_well) %>%
      mutate(transfer_n = transfer_n - 1) %>%
      distinct() %>%
      left_join(plate_data(),
        by = c("transfer_n", "infection_from_well" = "rwell")) %>%
      mutate(
        x = x_base + stretch_x * transfer_n + 2 / 3 * box_size,
        y = y_base + y_heights[3])
    return(infecting_wells)
  })

  infected_wells <- reactive({
    infected_wells <- filter(onewell(), !is.na(infection_to_well)) %>%
      select(transfer_n, infection_to_well) %>%
      mutate(transfer_n = transfer_n + 1) %>%
      distinct() %>%
      left_join(plate_data(),
        by = c("transfer_n", "infection_to_well" = "rwell")) %>%
      mutate(
        x = x_base + stretch_x * transfer_n - 2 / 3 * box_size,
        y = y_base + y_heights[1])
    return(infected_wells)
  })

  turnover_data <- reactive({
    turnover_data <- filter(onewell(), !is.na(turnover_strain)) %>%
      select(transfer_n, turnover_strain) %>%
      distinct() %>%
      mutate(
        x = stretch_x * transfer_n + 1 / 3 * box_size,
        y = 2.5 * stretch_y) %>%
      left_join(phenotype_table$strains, by = c("turnover_strain" = "strain"))
    return(turnover_data)
  })

  transfer_data <- reactive({
    transfer_data <- filter(onewell(), !is.na(transfer_to_well)) %>%
      select(transfer_n) %>%
      distinct() %>%
      mutate(
        x = stretch_x * (transfer_n - 1) + box_size / 2,
        y = 3/2 * stretch_y)
    return(transfer_data)
  })

  ##############################################################################
  # output - Plots                                                             #
  ##############################################################################

  output$plot_expectations <- renderPlot({
    expectations_plot <- filter(plate_data(), transfer_n == input$transfer) %>%
      ggplot(
        aes(col, row_int,
          fill = class_short,
          label = class_short)) +
        geom_tile() +
        geom_text() +
        scale_y_reverse("",
          expand = c(0, 0), breaks = 1:16, labels = LETTERS[1:16]) +
        scale_x_continuous("",
          expand = c(0, 0), breaks = 1:24) +
        geom_hline(yintercept = 1:7 * 2 + 0.5) +
        geom_vline(xintercept = 1:11 * 2 + 0.5) +
        scale_fill_manual(values = colors_changes,
          labels = labels_changes,
          guide = guide_legend(title = NULL)) +
        labs(title = "Classification") +
        annotate("rect",
          xmin = values$well$col, xmax = values$well$col + 2,
          ymin = values$well$row, ymax = values$well$row + 2,
          color = "magenta", fill = NA, size = 3) +
        theme_bw_custom_hm
    expectations_plot
  })

  output$plot_plate <- renderPlot({
    plate_plot <- filter(plate_data(), transfer_n == input$transfer) %>%
      ggplot(
        aes(col, row_int,
          fill = as.character(mix_phenotype),
          label = labels)) +
        geom_tile() +
        geom_text() +
        scale_y_reverse("",
          expand = c(0, 0), breaks = 1:16, labels = LETTERS[1:16]) +
        scale_x_continuous("",
          expand = c(0, 0), breaks = 1:24) +
        geom_hline(yintercept = 1:7 * 2 + 0.5) +
        geom_vline(xintercept = 1:11 * 2 + 0.5) +
        scale_fill_manual(values = colors_phenotypes,
          labels = labels_phenotypes,
          guide = guide_legend(title = NULL)) +
        labs(title = "expected phenotypes") +
        annotate("rect",
          xmin = values$well$col, xmax = values$well$col + 2,
          ymin = values$well$row, ymax = values$well$row + 2,
          color = "magenta", fill = NA, size = 3) +
        theme_bw_custom_hm
    plate_plot
  })

  output$plot_od <- renderPlot({
    od_plot <- filter(plate_data(), transfer_n == input$transfer) %>%
      distinct(plate, well, transfer_n, .keep_all = TRUE) %>% 
      ggplot(
        aes(col, row_int, fill = OD, label = signif(OD, 2))) +
        geom_tile() +
        geom_text() +
        scale_y_reverse("",
          expand = c(0, 0), breaks = 1:16, labels = LETTERS[1:16]) +
        scale_x_continuous("",
          expand = c(0, 0), breaks = 1:24) +
        scale_fill_gradientn(colours = hmPalette(100)) +
        geom_hline(yintercept = 1:7 * 2 + 0.5) +
        geom_vline(xintercept = 1:11 * 2 + 0.5) +
        labs(title = "Optical density 595nm") +
        annotate("rect",
          xmin = values$well$col, xmax = values$well$col + 2,
          ymin = values$well$row, ymax = values$well$row + 2,
          color = "magenta", fill = NA, size = 3) +
        theme_bw_custom_hm
    od_plot
  })

  output$plot_pickolo <- renderPlot({
    pickolo_plot <- filter(plate_data(), transfer_n == input$transfer) %>%
      ggplot(
        aes(col, row_int,
          fill = as.character(phenotype),
          label = intToBin(phenotype))) +
        geom_tile() +
        geom_text() +
        scale_y_reverse("",
          expand = c(0, 0), breaks = 1:16, labels = LETTERS[1:16]) +
        scale_x_continuous("",
          expand = c(0, 0), breaks = 1:24) +
        geom_hline(yintercept = 1:7 * 2 + 0.5) +
        geom_vline(xintercept = 1:11 * 2 + 0.5) +
        scale_fill_manual(values = colors_phenotypes,
          labels = labels_phenotypes,
          guide = guide_legend(title = NULL)) +
        annotate("rect",
          xmin = values$well$col, xmax = values$well$col + 2,
          ymin = values$well$row, ymax = values$well$row + 2,
          color = "magenta", fill = NA, size = 3) +
        labs(title = "Colonies") +
        theme_bw_custom_hm
    pickolo_plot
  })

  output$chain_plot <- renderPlot({
    chain_plot <- ggplot(onewell(),
      aes(x, y, fill = as.character(phenotype), label = labels)) +
      #highlight selected transfer
      annotate("rect",
        xmin = input$transfer * stretch_x,
        xmax = input$transfer * stretch_x + box_size,
        ymin = -Inf, ymax = Inf,
        alpha = .2) +
      #turnover
      geom_segment(data = turnover_data(),
        aes(x = x, y = y,
          xend = x + 1 / 6 * box_size, yend = y - stretch_y),
        inherit.aes = FALSE) +
      geom_point(data = turnover_data(),
        aes(x, y, color = as.character(phenotype)),
        size = 10,
        inherit.aes = FALSE) +
      geom_text(data = turnover_data(), aes(x, y, label = turnover_strain),
        inherit.aes = FALSE) +
      #transfer
      geom_segment(data = transfer_data(),
        aes(x = x, y = y, xend = x + stretch_x, yend = y),
        inherit.aes = FALSE) +
      # infecting wells
      geom_label(data = infecting_wells() %>% group_by(transfer_n) %>%
          summarize(
            x = first(transfer_n) * stretch_x + box_size / 2 + 2 / 3 * box_size,
            y_label = y_heights[3] + label_height,
            infection_from_well = first(infection_from_well),
            .groups = "drop_last"),
         aes(x, y_label, label = infection_from_well),
         inherit.aes = FALSE) +
      geom_segment(data =
          group_by(infecting_wells(), infection_from_well, transfer_n) %>%
          summarize(
            x_line = mean(x),
            x_line_end = x_line + stretch_x - 2 / 3 * box_size,
            y_line = mean(y),
            y_line_end = y_line - stretch_y,
            .groups = "drop_last"),
        aes(x = as.double(x_line), y = as.double(y_line),
          xend = as.double(x_line_end), yend = as.double(y_line_end)),
        inherit.aes = FALSE) +
      geom_tile(data = infecting_wells(),
        size = 0.25, colour = "black") +
      geom_text(data = infecting_wells()) +
      #infected wells
      geom_label(data = infected_wells() %>% group_by(transfer_n) %>%
          summarize(
            x = first(transfer_n) * stretch_x -
              2 / 3 * box_size + 1 / 2 * box_size,
            y_label = y_heights[1] + label_height,
            infection_to_well = first(infection_to_well),
            .groups = "drop_last"),
        aes(x, y_label, label = infection_to_well),
        inherit.aes = FALSE) +
      geom_segment(data =
          group_by(infected_wells(), infection_to_well, transfer_n) %>%
          summarize(
            x_line = mean(x) - stretch_x + 2 / 3 * box_size,
            x_line_end = x_line + stretch_x - 2 / 3 * box_size,
            y_line_end = mean(y),
            y_line = y_line_end + stretch_y,
            .groups = "drop_last"),
        aes(x = as.double(x_line), y = as.double(y_line),
          xend = as.double(x_line_end), yend = as.double(y_line_end)),
        inherit.aes = FALSE) +
      geom_tile(data = infected_wells(),
        size = 0.25, colour = "black") +
      geom_text(data = infected_wells()) +
      # treatment
      geom_label(data = treatment_data(),
        aes(x, y, label = labels),
        inherit.aes = FALSE) +
      # main well
      geom_tile(size = 0.25, colour = "black") +
      geom_text() +
      scale_x_continuous("",
        limits =
          c(0, 2 * stretch_y / 2 + stretch_x + 
            max(onewell()$transfer_n) * stretch_x),
        breaks = 0:(max(onewell()$transfer_n) + 1) * stretch_x + stretch_y / 2,
        labels = 0:(max(onewell()$transfer_n) + 1)) +
      scale_y_continuous("",
        limits = c(0, 3 * stretch_y),
        breaks = 3 * stretch_y / 2,
        labels = unique(onewell()$rwell)) +
      scale_fill_manual(values = colors_phenotypes,
        labels = labels_phenotypes,
        guide = guide_legend(title = NULL),
        na.value = "white") +
      scale_color_manual(values = colors_phenotypes,
        labels = labels_phenotypes,
        guide = guide_legend(title = NULL)) +
      theme_bw_custom_hm
    chain_plot
  })
  ##############################################################################
  # output                                                                     #
  ##############################################################################
  output$selectedWell <- renderDataTable({
        group_by(infected_wells(), infection_to_well, transfer_n) %>%
          summarize(
            x_line = mean(x) - stretch_x + 2 / 3 * box_size,
            x_line_end = x_line + stretch_x - 2 / 3 * box_size,
            y_line_end = mean(y),
            y_line = y_line_end + stretch_y,
            .groups = "drop_last")
  })

  output$test <- renderPrint({
    input$click_chain
    #nearPoints(onewell(), input$click_chain, threshold = 50, maxpoints = 1)
    #sel_w()
  })
})
