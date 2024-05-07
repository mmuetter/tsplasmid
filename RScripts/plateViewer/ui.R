shinyUI(fluidPage(title = "plate correct",
  fluidRow(
    tags$script('
        $(document).on("keypress", function (e) {
          Shiny.onInputChange("key_pressed", [e.which,e.timeStamp]);
        });'),
    column(1,
      tags$b("Plate Viewer")),
    column(2,
      selectizeInput(
        "exp", "Experiment", choices = exp_choices, multiple = FALSE,
        selected = max(exp_choices), width = "100%")),
    column(1,
      numericInput(
        "plate", "Plate", value = 1, min = 1, max = 6, width = "100%")),
    column(1,
      numericInput(
        "transfer", "Transfer", value = 0, min = 0, width = "100%"))
  ),
  fluidRow(
    column(12,
      plotOutput("chain_plot", width = "100%", height = "350px",
        click = "click_chain"))
  ),
  fluidRow(
    column(6,
      plotOutput("plot_plate", width = "825px", height = "550px",
        click = "click")),
    column(6,
      plotOutput("plot_pickolo", width = "825px", height = "550px",
        click = "click"))
  ),
  fluidRow(
    column(6,
      plotOutput("plot_od", width = "825px", height = "550px",
        click = "click")),
    column(6,
      plotOutput("plot_expectations", width = "875px", height = "550px",
        click = "click"))
  ),
  fluidRow(
    column(12,
      verbatimTextOutput("test"))
  )
  # fluidRow(
  #   column(12,
  #     dataTableOutput("selectedWell"))
  # )
))
