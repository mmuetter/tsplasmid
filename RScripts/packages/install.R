
install.packages("here")

packagePath <- here::here("RScripts/packages")

# Dependencies
install.packages(c(
  "languageserver",
  "xml2",
  "jsonlite",
  "tidyverse",
  "stringr",
  "jpeg",
  "rmarkdown",
  "knitr",
  "shiny",
  "TSP"
  ))

# custom packages
install.packages(paste0(packagePath, "/rrobot_1.2.tar.gz"))
install.packages(paste0(packagePath, "/treatR_0.2.2.tar.gz"))
