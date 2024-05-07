shinyServer(function(input, output, session) {
  ##############################################################################
  # dynamic UI                                                                 #
  ##############################################################################

  output$saveRedrawButton <- renderUI({
    if (rvals$saveOrRedraw){
      buttonLabel <- "save (s)"
    } else {
      buttonLabel <- "redraw (s)"
    }
    actionButton("saveRedrawButton", label = buttonLabel)
  })

  ##############################################################################
  # Observers                                                                  #
  ##############################################################################
  observe({
    rvals$reloadData

    pickoloFiles <- tools::file_path_sans_ext(list.files(
      str_c(pickolo_path(), "/csv")))

    pickoloDf <- str_split(pickoloFiles, "_|t|P", simplify = TRUE) %>%
      as_tibble(.name_repair = ~ make.names(str_c("V", 1:length(.)))) %>%
      select(-V1, -V3) %>%
      transmute(name = pickoloFiles, transfer = as.integer(V2), plate = as.integer(V4), drug = V5, bc = V6) %>%
      arrange(transfer, plate)

    if (!input$show_done) {
      jsonFiles <- tools::file_path_sans_ext(list.files(
        str_c(pickolo_path(), "/json")))
      if (length(jsonFiles) > 0) {
        doneFiles <- str_split(jsonFiles, "_|t|P", simplify = TRUE)[, 6]
        pickoloDf <- pickoloDf %>% filter(!(bc %in% doneFiles))
      }
      inputTransfer <- pickoloDf$transfer[1]
      inputPlate <- pickoloDf$name[1]
    } else {
      if (is.na(input$transfer) | input$transfer == "") {
        inputTransfer <- pickoloDf$transfer[1]
      } else {
        inputTransfer <- input$transfer
      }
      if (is.na(input$plate) | input$plate == "") {
        inputPlate <- pickoloDf$name[1]
      } else {
        inputPlate <- input$plate
      }
    }

    selectablePlates <- pickoloDf$name[pickoloDf$transfer == inputTransfer]

    updateSelectizeInput(session, "transfer",
      choices = unique(pickoloDf$transfer),
      selected = inputTransfer)
    updateSelectizeInput(session, "plate",
      choices = selectablePlates,
      selected = inputPlate)
  })

  observeEvent(input$click_agar, {
    rvals$colfix <- colony_correct(
      input$click_agar$x,
      input$click_agar$y,
      pickolo_data(),
      colfix = rvals$colfix)
    if (rvals$saveOrRedraw) {
      rvals$saveOrRedraw <- FALSE
    }
  })

  observeEvent(input$key_pressed[2], {
    if (input$key_pressed[1] == 115) {
      #S
      if (rvals$saveOrRedraw) {
        rvals$saveButton <- rvals$saveButton + 1
      } else {
        rvals$redrawButton <- rvals$redrawButton + 1
      }
      rvals$saveOrRedraw  <- !rvals$saveOrRedraw
    }
    if (input$key_pressed[1] == 114) {
      #R
      rvals$redrawButton <- rvals$redrawButton + 1
    }
  })

  observeEvent(input$save, {
    if (rvals$saveOrRedraw){
        rvals$saveButton <- rvals$saveButton + 1
      } else {
        rvals$redrawButton <- rvals$redrawButton + 1
      }
    rvals$saveOrRedraw  <- !rvals$saveOrRedraw
  })

  observeEvent(rvals$saveButton, {
    pData <- pickolo_data()
    fixdata <- colfix_current()
    saveRDS(fixdata, "fixdataTest.rds")
    pData$growth[match(fixdata$well, pData$growth$well), ] <- select(fixdata,
      well, row, col, auto,  manual)
    write(toJSON(pData, pretty = TRUE, na = "string"),
      json_path())
    rvals$reloadData <- rvals$reloadData + 1
    rvals$colfix <- colfix_blank()
    rvals$colfixPlot <- colfix_blank()
  })

  observeEvent(rvals$redrawButton, {
    rvals$colfixPlot <- colfix_current()
  })

  ##############################################################################
  # Data                                                                       #
  ##############################################################################
  rvals <- reactiveValues(
    colfix = colfix_blank(),
    colfixPlot = colfix_blank(),
    reloadData = 0,
    saveButton = 0,
    redrawButton = 0,
    saveOrRedraw = FALSE)

  pickolo_path <- reactive({
    output_directory <- filter(experiments, id == input$exp)$output_directory
    pickolo_path <- str_c("../../experiments/", output_directory, "/pickolo")
    return(pickolo_path)
  })

  csv_path <- reactive({
    csv_path <- str_c(pickolo_path(), "/csv/", input$plate, ".csv")
    # if (str_detect(csv_path,"_A_")) {
    #   updateNumericInput(session, "row_start_coord", value = 180)
    #   updateNumericInput(session, "col_start_coord", value = 160)
    # } else {
    #   updateNumericInput(session, "row_start_coord", value = 190)
    #   updateNumericInput(session, "col_start_coord", value = 150)
    # }
    return(csv_path)
  })

  json_path <- reactive({
    json_path <- str_c(pickolo_path(), "/json/", input$plate, ".json")
    return(json_path)
  })

  pickolo_data <- reactive({
    shiny::validate(
      need(file.exists(csv_path()), "no plates left to correct!")
      )
    rvals$reloadData
    if (file.exists(json_path())) {
      pickolo_data <- fromJSON(json_path())
    } else {
      pickolo_data <- readColonies(csv_path(),
        row_start_coord = input$row_start_coord,
        col_start_coord = input$col_start_coord,
        row_d = input$row_d,
        col_d = input$col_d)
      pickolo_data$pars$csv_path <- csv_path()
      pickolo_data$pars$img_path <- str_c(pickolo_path(),
        "/img/", input$plate, ".jpg")
      return(pickolo_data)
    }
  })

  pickolo_img <- reactive({
    pickolo_img <- jpeg::readJPEG(pickolo_data()$pars$img_path)
  })

  colfix_current <- reactive({
    colfix_current <- filter(rvals$colfix, csv_path == csv_path())
    return(colfix_current)
  })

  ##############################################################################
  # Data Output                                                                #
  ##############################################################################

  output$test <- renderPrint({
    input$key_pressed[1]
  })

  output$colFix <- renderPrint({
    as.data.frame(select(colfix_current(), -x, -y, -csv_path))
  })

  ##############################################################################
  # Plot Output                                                                #
  ##############################################################################

  output$plot_agar <- renderPlot({
    plotColonies2(
      pickoloData = pickolo_data(),
      img_path = pickolo_data()$pars$img_path,
      img = pickolo_img(),
      plotLabels = input$show_labels, colChange = rvals$colfixPlot,
      img_opts$line_size,
      img_opts$title_size,
      img_opts$label_size,
      img_opts$axis_text_size)
  })

})
