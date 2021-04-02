# A script depicting the channel capacity calculation 
# for a scenario 2 from a paper:

# "Phenotypic variability and not noise accounts for most of the cell-to-cell
# heterogeneity in cytokine signaling", Topolewski et al. 2021

## scenario 1: Population response distributions ##

#### loading libraries and functions ####
library(dplyr)

#install the SLEMI package according to https://github.com/sysbiosig/SLEMI
library(SLEMI) 

get.single.percentile.capacity <- function(distribution.type,
                                           distribution.parameters, 
                                           n.sampled){
  # a function to calculate channel capacity for n.sampled values 
  # sampled from distribution.type which has given distribution.parameters. 
  # The level of input (stimulations) is given by the distribution.parameters argument:
  # - each row of distribution.parameters correspond to differenet input level; 
  # - each colum to specific distribution parameter (has to be two columns in total);
  # only "gamma" or "lognormal" values allowed for distribution.type

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
      dplyr::filter(stim_level == ss) 
    D[[i]] <- quantile(X_S$response,percentiles);
    i=i+1
  }
  
  # st.dev. as estimated from logistic regression
  if(distribution.type == "gamma"){
    ldlist=list() # parameters of gamma for different percentiles, shape and rate
    for(k in 1:length(percentiles)){
      ld = matrix(0, length(Uss), 2)
      for(l in 1:length(Uss)){
        ld[l,] = c(1/sd.coefficient^2, 1 / (sd.coefficient^2*D[[l]][k]))
        ldlist[[k]]=ld
      }
      
    }
  } else if(distribution.type == "lognormal"){
    ldlist=list()
    for(k in 1:length(percentiles)){
      ld = matrix(0, length(Uss), 2)
      for(l in 1:length(Uss)){
        mean.log = log(D[[l]][k]^2 / (D[[l]][k]^2 + (sd.coefficient * D[[l]][k])^2)^0.5)
        sd.log = (log(1+((sd.coefficient * D[[l]][k])^2 / D[[l]][k]^2)))^0.5
        ld[l,] = c(mean.log, sd.log)
        ldlist[[k]]=ld
      }
    }
  } else {
    return("wrong distribution type, either gamma or lognormal possible")
  }
  
  if(n.sampled >= 1){
    doses <- log10(c(0.001, 0.1, 1, 10))
    capacities <- c()
    
    for(k in 1:length(percentiles)){
      if(continuous == "FALSE"){
        
        capacities = c(capacities, 
                       get.single.percentile.capacity(
                         distribution.parameters = ldlist[[k]], 
                         n.sampled = 
                           n.sampled,
                         distribution.type = 
                           distribution.type)
        )
      } else if(continuous == "TRUE"){
        rvec <- unlist(D)[seq(k, length(percentiles) * length(doses), by = length(percentiles))]
        capacities <- c(capacities,
                        get.single.percentile.capacity.cont(
                          xvec = doses, 
                          rvec = rvec, 
                          sd.coefficient2 = sd.coefficient^2,
                          distribution.type = distribution.type)
        )
        
      } else {return("wrong continuous value, should be T/F")}
      
    }
    capacity <- log2(mean(2^capacities))
  } else {
    return("Too low sample number!")
  }
  return(capacity)
  
}

sd.AB <- function(A, B){
  # function to calculate sd from two df.columns, 
  # representing signal in the nucleus A and B
  if(length(A) != length(B)){
    break()
  }
  sd <- c()
  for(i in 1:length(A)){
    mean <- (A[i] + B[i])/2 
    sd[i] <- (((A[i]-mean)^2 + (B[i]-mean)^2)/2)^0.5
  }
  return(sd)
}


#### set-up of directories ####
base.path <- "D:/Piotrek/publications/syncytia_noise/phenotypic_variability"
input.path <- paste(base.path, "/data", sep = "")

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
  
  bio.for.capacity <- bio.traj %>%
    dplyr::filter(time == time.chosen &
                    stim_type == stim.type.chosen)
  
  stimulant.gamma.capacity <-
    get.all.percentile.capacity(data = bio.for.capacity,
                                sd.coefficient = sd.coefficient,
                                percentiles = seq(0.05, 0.95, by = 0.01),
                                n.sampled = 1000,
                                continuous = FALSE,
                                distribution.type = "gamma")
  
  stimulant.log.capacity <-
    get.all.percentile.capacity(data = bio.for.capacity,
                                sd.coefficient = sd.coefficient,
                                percentiles = seq(0.05, 0.95, by = 0.01),
                                n.sampled = 1000,
                                continuous = FALSE,
                                distribution.type = "lognormal")
  
  all.capacity <- rbind(all.capacity,
                        data.frame(stim_type = stim.type.chosen,
                                   gamma.capacity = stimulant.gamma.capacity,
                                   lognormal.capacity = stimulant.log.capacity)) 
}
print(all.capacity)
