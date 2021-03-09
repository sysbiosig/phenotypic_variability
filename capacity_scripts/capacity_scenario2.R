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
project <- "fusion"
normaliz.chosen <- "pbs"
path <- paste(base.path, folder.name, sep='')

#### data loading ####
bio.AB <- read.csv(paste(path, "/", project, 
                         "/input/", normaliz.chosen,
                         "/fused_paper.csv", sep = ""), header = TRUE, sep = ",")

merge.chosen <- "no"
bio.AB <- bio.AB %>% 
  dplyr::filter(merged_nuclei == merge.chosen)

bio.traj <- read.csv(paste(path, "/trajectories", 
                           "/input/", normaliz.chosen,
                           "/trajectories_paper.csv", sep = ""), header = TRUE, sep = ",")
bio.traj <- bio.traj %>% 
  dplyr::filter(nuclearity == "single-nuclear")


#### choosing the stimulants ####
chosen.stimulants <- c("IFNG", "OSM")

bio.coefficients <- bio.AB %>%
  group_by(stim_type)%>%
  mutate(sd = sd.AB(A, B),
         mean = (A + B)/2) %>%
  summarise(a = summary(lm(sd ~ 0 + mean))$coefficients[[1]]) %>%
  ungroup()

gamma.capacity <- data.frame()
for(stim.type.chosen in chosen.stimulants){
  
  if(stim.type.chosen == "IFNG"){
    time.chosen <- 15
  } else if(stim.type.chosen == "OSM") { 
    time.chosen <- 30
  }
  
  sd.coefficient <- as.numeric(bio.coefficients %>%
                                 dplyr::filter(stim_type == stim.type.chosen) %>%
                                 select(a))
  
  CM.capacity <- confusion.matrix.gamma(data = bio.traj,
                                        stim_type_chosen = stim.type.chosen,
                                        time.chosen = time.chosen,
                                        sd.coefficient = sd.coefficient,
                                        percentiles = seq(0.05, 0.95, by = 0.01),
                                        n.sampled = 10)
  
  gamma.capacity <- rbind(gamma.capacity,
                          data.frame(stim_type = stim.type.chosen,
                                     gamma.capacity = CM.capacity$capacity)) 
}
print(gamma.capacity)
