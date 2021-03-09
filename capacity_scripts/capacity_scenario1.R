# A script depicting the capacity calculation for a scenario 1 from a paper:

# "Phenotypic variability and not noise accounts for most of the cell-to-cell
# heterogeneity in cytokine signaling", Topolewski et al. 2021

## scenario 1: Population response distributions ##

#### loading libraries and functions ####
source("D:/Piotrek/scripts/basic_scripts/normalize_data.R")
library(SLEMI)

#### set-up of directories ####
base.path <- "D:/Piotrek/Experiments/ICF/"
folder.name <- "cellular_identity"
project <- "trajectories"
normaliz.chosen <- "pbs"
path <- paste(base.path, folder.name, sep='')

#### data loading ####
bio.traj <- read.csv(paste(path, "/", project, "/input/", 
                           normaliz.chosen,
                           "/trajectories_paper.csv", sep = ""),  ## if you want to plot TNF, remove "_paper"
                     header = TRUE, sep = ",")

bio.traj <- bio.traj %>%
  dplyr::filter(merged_nuclei %in% c(NA))

#### choosing the stimulants ####
chosen.stimulants <- c("IFNG", "OSM")

bio.capacity <- bio.traj %>%
  dplyr::filter((time == 15 & stim_type == "IFNG") |
                  (time == 30 & stim_type == "OSM")) %>%
  group_by(stim_type) %>%
  summarise(capacity = SLEMI::capacity_logreg_main(dataRaw = cur_data(),
                                                   signal = "stim_level",
                                                    response = "response")$cc)

print(bio.capacity)

