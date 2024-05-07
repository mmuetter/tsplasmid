shinyUI(
  fluidPage(title = "plate correct",
    fluidRow(
      tags$script('
        $(document).on("keypress", function (e) {
          Shiny.onInputChange("key_pressed", [e.which,e.timeStamp]);
        });'),
      column(2,
        tags$b("Plate Corrector")),
      column(2,
        selectizeInput(
          "exp", "Experiment", choices = exp_choices, multiple = FALSE,
          selected = max(exp_choices), width = "100%")),
      column(2,
        selectizeInput(
          "transfer", "Transfer", choices = "", multiple = FALSE,
          selected = "", width = "100%")),
      column(3,
        selectizeInput(
          "plate", "Plate", choices = "", multiple = FALSE,
          selected = "", width = "100%")),
      column(3,
        checkboxInput("show_done", "Show allready corrected?", value = FALSE))
    ),
    fluidRow(
      tabsetPanel(type = "tabs", id = "tabs",
        tabPanel("Plot",
          column(9,
            plotOutput("plot_agar", width = "1200px", height = "900px",
              click = clickOpts(id = "click_agar", clip = TRUE))),
          column(3,
            # verbatimTextOutput("test"),
            verbatimTextOutput("colFix"),
            uiOutput("saveRedrawButton"),
            tags$hr(),
            numericInput("row_start_coord", "row_start_coord", 180,
              min = 0, max = 1200, step = 1),
            numericInput("col_start_coord", "col_start_coord", 160,
              min = 0, max = 1600, step = 1),
            numericInput("row_d", "row_d", 56.5,
              min = 0, step = 1),
            numericInput("col_d", "col_d", 56.5,
              min = 0, step = 1),
            checkboxInput("show_labels", "Labels?", value = FALSE))
        ),
        tabPanel("Data",
          column(12,
            )
        )
      )
    )
  )
)
