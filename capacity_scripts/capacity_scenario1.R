# A script depicting the channel capacity calculation 
# for a scenario 1

## scenario 1: Population response distributions ##

#### loading libraries and functions ####
library(dplyr)
#install the SLEMI package according to https://github.com/sysbiosig/SLEMI
library(SLEMI) 

#### set-up of directories ####
# either provide a path to the main directory, 
# or open the "phenotypic_variability" Rproject
base.path <- getwd()
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

