# A script depicting the channel capacity calculation 
# for a scenario 3 

## scenario 3: Single-cell noise response distributions with continuous input ##

#### set-up of directories ####
# either provide a path to the main directory, 
# or open the "phenotypic_variability" Rproject
base.path <- getwd()
input.path <- paste(base.path, "/data", sep = "")

#### loading libraries and functions ####
library(dplyr)

#install the SLEMI package according to https://github.com/sysbiosig/SLEMI
library(SLEMI) 

sqrtFIvec_gamma <- function(dose,
                         sd.coefficient2,
                         res){
  
  me=res(dose) # mean of the response for specifc dose 
  d_mean_d_dose=res(dose, deriv = 1) # derivative of the response with respect to the log of the dose
  shape=1/sd.coefficient2 # shape parameter of gamma distribution for mean mi and variance si2
  rate=1/(me*sd.coefficient2) # rate parameter of gamma distribution for mean mi and variance si2
  
  d_shape_d_dose=0 # derivative of the shape parameter with respect to the log of the dose
  d_rate_d_dose=-1/(sd.coefficient2*me*me)*d_mean_d_dose # derivative of the rate parameter with respect to the log of the dose
  
  FIGamma_11=trigamma(shape) # element 11 of the Fisher informtion matrix of the gamma distribution i.e. with respoect to shape
  FIGamma_22=shape*rate^(-2)  # element 22 of the Fisher informtion matrix of the gamma distribution i.e. with respoect to rate
  FIGamma_12=-rate^(-1) # element 12/21 of the Fisher informtion matrix of the gamma distribution i.e. with respoect to shape and rate
  v_1=d_shape_d_dose #element 1 the re-paramtrization vector
  v_2=d_rate_d_dose  #element 2 the re-paramtrization vector
  
  sqrt(FIGamma_11*v_1*v_1+FIGamma_22*v_2*v_2+2*FIGamma_12*v_1*v_2) # sqrt of the Fisher information with respect to log of the dose
}

sqrtFIvec_lognormal <- function(dose,
                             sd.coefficient2,
                             res){
  d_logmean_d_dose=res(dose, deriv = 1)/res(dose) # derivative of the log of the response with respect to the log of the dose
  FI=1/log(1+sd.coefficient2)*(d_logmean_d_dose)^2 # Fisher information of the log-normal distribution
  sqrt(FI)  
}

get.single.percentile.capacity.cont <- function(xvec,
                                             rvec,
                                             sd.coefficient2,
                                             distribution.type){

  # interpolation of the response with respect to log of the dose
  res = splinefun(xvec, rvec, method="hyman") 
  lower <- -3
  # integration of the square root of the Fisher information matrix
  if(distribution.type == "gamma"){
    resstates = integrate(sqrtFIvec_gamma, lower, 2, sd.coefficient2, res) 
  } else if(distribution.type == "lognormal"){
    resstates=integrate(sqrtFIvec_lognormal, lower, 2, sd.coefficient2, res)
  } else {
    return("wrong distribution type, either gamma or lognormal possible")
    }
  
  # Shannon capacity
  log2(resstates$value/sqrt(2*pi*exp(1)))   
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

#### data loading ####
bio.AB <- read.csv(paste(input.path,
                         "/fused_paper.csv", sep = ""), header = TRUE, sep = ",")

merge.chosen <- "no"
bio.AB <- bio.AB %>% 
  dplyr::filter(merged_nuclei == merge.chosen &
                  ((stim_type == "IFNG" &
                      time == 15) |
                     (stim_type == "OSM" &
                        time == 30)))

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
  
  stimulant.gamma.capacity.c <-
    get.all.percentile.capacity(data = bio.for.capacity,
                                sd.coefficient = sd.coefficient,
                                percentiles = seq(0.05, 0.95, by = 0.01),
                                continuous = TRUE,
                                distribution.type = "gamma")
  
  stimulant.log.capacity.c <-
    get.all.percentile.capacity(data = bio.for.capacity,
                                sd.coefficient = sd.coefficient,
                                percentiles = seq(0.05, 0.95, by = 0.01),
                                continuous = TRUE,
                                distribution.type = "lognormal")
  
  all.capacity <- rbind(all.capacity,
                        data.frame(stim_type = stim.type.chosen,
                                   gamma.capacity.cont = stimulant.gamma.capacity.c,
                                   lognormal.capacity.cont = stimulant.log.capacity.c))  
}
print(all.capacity)
