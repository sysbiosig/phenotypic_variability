# A script depicting the channel capacity calculation 
# for a scenario 1 from a paper:

# "Phenotypic variability and not noise accounts for most of the cell-to-cell
# heterogeneity in cytokine signaling", Topolewski et al. 2021

## scenario 1: Population response distributions ##

#### loading libraries and functions ####
library(SLEMI)
library(dplyr)

#### set-up of directories ####
base.path <- "D:/Piotrek/publications/syncytia_noise/phenotypic_variability"
input.path <- paste(base.path, "/data", sep = "")

#### data loading ####
bio.traj <- read.csv(paste(input.path,
                           "/trajectories_paper.csv", sep = ""),
                     header = TRUE, sep = ",")

bio.traj <- bio.traj %>%
  dplyr::filter(merged_nuclei %in% c(NA)) # single, unfused cells

#### calculating the channel capacity for scenario 1 ####
bio.capacity <- bio.traj %>%
  dplyr::filter((time == 15 & stim_type == "IFNG") |
                  (time == 30 & stim_type == "OSM")) %>%
  group_by(stim_type) %>%
  summarise(capacity = SLEMI::capacity_logreg_main(dataRaw = cur_data(),
                                                   signal = "stim_level",
                                                    response = "response")$cc)

print(bio.capacity)

