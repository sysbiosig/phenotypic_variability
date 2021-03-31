# A script depicting the channel capacity calculation 
# for a scenario 2 from a paper:

# "Phenotypic variability and not noise accounts for most of the cell-to-cell
# heterogeneity in cytokine signaling", Topolewski et al. 2021

## scenario 1: Population response distributions ##

#### set-up of directories ####
base.path <- "D:/Piotrek/publications/syncytia_noise/phenotypic_variability"
input.path <- paste(base.path, "/data", sep = "")

#### loading libraries and functions ####
# first load all required custom functions: 
source(paste(base.path, "/Topolewski_auxillary_functions.R", sep = ""))

get.single.percentile.capacity <- function(distribution.parameters, 
                                           n.sampled,
                                           distribution.type){
  library(SLEMI)
  
  stim.level.n = nrow(distribution.parameters)
  all.stim.data <- data.frame()
  if(distribution.type == "gamma"){
    for(i in 1:stim.level.n){
      
      single.stim.data <- 
        data.frame(response = 
                     rgamma(n = n.sampled, 
                            shape = distribution.parameters[i, 1], 
                            rate = distribution.parameters[i,2]),
                   stim_level = i)
      
      all.stim.data <- rbind(all.stim.data,
                             single.stim.data)
    }
  } else if(distribution.type == "lognormal"){
    for(i in 1:stim.level.n){
      single.stim.data <-
        data.frame(response =
                     rlnorm(n = n.sampled,
                            meanlog = distribution.parameters[i, 1],
                            sdlog = distribution.parameters[i,2]),
                   stim_level = i)
      
      all.stim.data <- rbind(all.stim.data,
                             single.stim.data)
    }
  } else {
    return("wrong distribution type, either gamma or lognormal possible")
  }
  
  
  return(SLEMI::capacity_logreg_main(dataRaw = all.stim.data, 
                                     signal = "stim_level", 
                                     response = "response")$cc)
}

get.all.percentile.capacity <- function(data,
                                        time.chosen,
                                        stim_type_chosen,
                                        sd.coefficient,
                                        percentiles = seq(0.05,0.95,by=0.01),
                                        n.sampled = 1000,
                                        continuous = FALSE,
                                        distribution.type = "gamma"){
  
  
  D=list() # D is a list of densities correspoding to every stimulation level
  i=1
  
  Uss <- sort(unique(data$stim_level))
  for (ss in sort(Uss)){# for every stimulation level ss
    X_S=data %>% 
      dplyr::filter(stim_level == ss,
                    time == time.chosen,
                    stim_type == stim_type_chosen) 
    D[[i]] <- quantile(X_S$response,percentiles);
    i=i+1
  }
  
  # st.dev. as estimated from logistic regression
  if(distribution.type == "gamma"){
    ldlist=list() # parameters of gamma for different percentiles, shape and rate
    for(k in 1:length(percentiles)){
      ld = matrix(0, length(Uss), 2)
      for(l in 1:length(Uss)){
        ld[l,] = c(1/sd.coefficient^2,1 / (sd.coefficient^2*D[[l]][k]))
        ldlist[[k]]=ld
      }
      
    }
  } else if(distribution.type == "lognormal"){
    ldlist=list()
    for(k in 1:length(percentiles)){
      ld = matrix(0, length(Uss), 2)
      for(l in 1:length(Uss)){
        mean.log = log(D[[l]][k]^2 / (D[[l]][k]^2 + (sd.coefficient * D[[l]][k])^2)^0.5)
        sd.log = log(1+((sd.coefficient * D[[l]][k])^2 / D[[l]][k]^2))
        ld[l,] = c(mean.log, sd.log) # to be provided
        ldlist[[k]]=ld
      }
      
    }
  } else {
    return("wrong distribution type, either gamma or lognormal possible")
  }
  
  if(n.sampled >= 1){
    doses <- log10(c(0.01, 0.1, 1, 10))
    capacities <- c()
    
    for(k in 1:length(percentiles)){
      if(continuous == "FALSE"){
        
        capacities = c(capacities, 
                       get.single.percentile.capacity(distribution.parameters =ldlist[[k]], 
                                                      n.sampled = 
                                                        n.sampled,
                                                      distribution.type = 
                                                        distribution.type))
      } else if(continuous == "TRUE"){
        capacities <- c(capacities,
                        get.capacity.gamma.c(xvec = doses, 
                                             rvec = 1/(ldlist[[k]][, 2] * sd.coefficient^2), 
                                             sd.coefficient2 = sd.coefficient^2))
        
      } else {return("wrong continuous value, should be T/F")}
      
    }
    capacity <- log2(mean(2^capacities))
  } else {
    return("Too low sample number!")
  }
  return(list(capacity = capacity))
  
}


#### data loading ####
bio.AB <- read.csv(paste(input.path,
                         "/fused_paper.csv", sep = ""), header = TRUE, sep = ",")

merge.chosen <- "no"
bio.AB <- bio.AB %>% 
  dplyr::filter(merged_nuclei == merge.chosen)

bio.traj <- read.csv(paste(input.path,
                           "/trajectories_paper.csv", sep = ""), header = TRUE, sep = ",")
bio.traj <- bio.traj %>% 
  dplyr::filter(nuclearity == "single-nuclear")


#### choosing the stimulants ####
bio.coefficients <- bio.AB %>%
  group_by(stim_type)%>%
  mutate(sd = sd.AB(A, B),
         mean = (A + B)/2) %>%
  summarise(a = summary(lm(sd ~ 0 + mean))$coefficients[[1]]) %>%
  ungroup()

set.seed(12)
all.capacity <- data.frame()
chosen.stimulants <- unique(bio.traj$stim_type)
for(stim.type.chosen in chosen.stimulants){
  
  if(stim.type.chosen == "IFNG"){
    time.chosen <- 15
  } else if(stim.type.chosen == "OSM") { 
    time.chosen <- 30
  }
  
  sd.coefficient <- as.numeric(bio.coefficients %>%
                                 dplyr::filter(stim_type == stim.type.chosen) %>%
                                 select(a))
  
  
  stimulant.gamma.capacity <- 
    get.all.percentile.capacity(data = bio.traj,
                                stim_type_chosen = stim.type.chosen,
                                time.chosen = time.chosen,
                                sd.coefficient = sd.coefficient,
                                percentiles = seq(0.05, 0.95, by = 0.01),
                                n.sampled = 100,
                                continuous = FALSE,
                                distribution.type = "gamma")
  
  stimulant.log.capacity <-
    get.all.percentile.capacity(data = bio.traj,
                                stim_type_chosen = stim.type.chosen,
                                time.chosen = time.chosen,
                                sd.coefficient = sd.coefficient,
                                percentiles = seq(0.05, 0.95, by = 0.01),
                                n.sampled = 100,
                                continuous = FALSE,
                                distribution.type = "lognormal")
  
  all.capacity <- rbind(all.capacity,
                        data.frame(stim_type = stim.type.chosen,
                                   gamma.capacity = stimulant.gamma.capacity,
                                   lognormal.capacity = stimulant.log.capacity)) 
}
print(all.capacity)
