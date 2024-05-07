################################################################################
#                                                                              #
################################################################################
if (Sys.info()[1] == "Windows"){
  working_directory <- "D:/Polybox/robot-data/Robot-Shared/TreatmentStrategies"
} else {
  working_directory <-
    "/Users/Daniel/polybox/Shared/Robot-Shared/TreatmentStrategies"
}
setwd(working_directory)
library(jsonlite)
library(tidyverse)

strains <- data_frame(
  strain = c("A_r", "B_r", "AB_r", "wt", "UI", "bl"),
  phenotype = as.integer(c(3, 5, 15, 1, 0, 0)),
  color = c("#ebe507", "#ffa10a", "#acdf0c", "#0cb7eb", "gray", "gray")
  )

phenotypes <- data_frame(
  int = 0:15,
  let = c("E", "W", "F", "A", "F", "B", "F", "A/B", "F", "F", "F", "F", "F",
    "F", "F", "AB"),
  valid = let != "F",
  color = c(
    "gray",
    "#0cb7eb",
    "#e53b6a",
    "#ebe507",
    "#e53b6a",
    "#ffa10a",
    "#e53b6a",
    "#662f89",
    "#e53b6a",
    "#e53b6a",
    "#e53b6a",
    "#e53b6a",
    "#e53b6a",
    "#e53b6a",
    "#e53b6a",
    "#acdf0c"))

# add growth in plate to phenotype -> phenotypes from 0-31
ext_phenotypes <- read_csv("/Users/daniel/polybox/Shared/Robot-Shared/treatmentstrategies/RScripts/JSON_generator_fns/ext_phenotypes.csv")
# Definition of transition matrix
# Define expected and unexpected changes, etc
changes_legend <- data_frame(
  class_short = c("X","C","F","S","U","R","P","v","V"),
  class_long = c("expected", "contamination", "Error?", "Treatment successful",
    "Tx unsucessfull", "denovo resistance", "replacement", "single reversion",
    "double reversion"),
  color = c("green3", "red", "red3", "green1", "green2", "orange", "orange1",
    "orange2","orange3"))

agar_phenotype_changes <- tribble(
  ~pre,~post,~none,~A,~B,~AB,
  "E", "E"   , "X", "X", "X", "X", 
  "E", "W"   , "C", "C", "C", "C", 
  "E", "A"   , "C", "C", "C", "C", 
  "E", "B"   , "C", "C", "C", "C", 
  "E", "A/B" , "C", "C", "C", "C", 
  "E", "AB"  , "C", "C", "C", "C", 
  
  "W", "E"   , "F", "S", "S", "S", 
  "W", "W"   , "X", "U", "U", "U", 
  "W", "A"   , "F", "R", "F", "F", 
  "W", "B"   , "F", "F", "R", "F", 
  "W", "A/B" , "F", "F", "F", "F", 
  "W", "AB"  , "F", "F", "F", "R", 
  
  "A", "E"   , "F", "F", "S", "S", 
  "A", "W"   , "v", "F", "F", "F", 
  "A", "A"   , "X", "X", "U", "U", 
  "A", "B"   , "F", "F", "F", "F", 
  "A", "A/B" , "F", "F", "F", "F", 
  "A", "AB"  , "F", "F", "R", "R",  
  
  "B", "E"   , "F", "S", "F", "S", 
  "B", "W"   , "v", "F", "F", "F", 
  "B", "A"   , "F", "F", "F", "F", 
  "B", "B"   , "X", "U", "X", "U", 
  "B", "A/B" , "F", "F", "F", "F", 
  "B", "AB"  , "F", "R", "F", "R",  
  
  "A/B", "E"   , "F", "F", "F", "S", 
  "A/B", "W"   , "F", "F", "F", "F", 
  "A/B", "A"   , "P", "X", "F", "F", 
  "A/B", "B"   , "P", "F", "X", "F", 
  "A/B", "A/B" , "X", "U", "U", "U", 
  "A/B", "AB"  , "F", "F", "F", "R", 
  
  "AB", "E"   , "F", "F", "F", "F", 
  "AB", "W"   , "V", "F", "F", "F", 
  "AB", "A"   , "v", "v", "F", "F", 
  "AB", "B"   , "v", "F", "v", "F", 
  "AB", "A/B" , "F", "F", "F", "F", 
  "AB", "AB"  , "X", "X", "X", "X"
) %>% gather("treatment","class_short",3:6) %>%
  mutate(
    pre_int = phenotypes$int[match(pre,phenotypes$let)],
    post_int = phenotypes$int[match(post,phenotypes$let)])

# G = detectable OD, N = non-detectable OD
od_phenotype_changes <- tribble(
  ~last_phenotype,~od,~none,~A,~B,~AB,
  "E","G"      ,"C","C","C","C",
  "E","N"      ,"X","X","X","X",
  "W","G"      ,"X","U","U","U",
  "W","N"      ,"F","S","S","S",
  "A","G"      ,"X","X","U","U",
  "A","N"      ,"F","F","S","S",
  "B","G"      ,"X","U","X","U",
  "B","N"      ,"F","S","F","S",
  "A/B","G"    ,"X","X","X","U",
  "A/B","N"    ,"F","F","F","S",
  "AB","G"     ,"X","X","X","X",
  "AB","N"     ,"F","F","F","F"
) %>% gather("treatment","class_short",3:6) %>%
  mutate(
    last_phenotype_int = phenotypes$int[match(last_phenotype,phenotypes$let)])

phenotypes_json <- lst(
  strains,
  phenotypes,
  ext_phenotypes,
  changes_agar = agar_phenotype_changes,
  changes_od = od_phenotype_changes,
  changes_legend)

write(toJSON(phenotypes_json, pretty = TRUE),
  "experiments/templates/phenotypes.json")
