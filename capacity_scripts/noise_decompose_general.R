noise.decompose <- function(data,
                            colname.A,
                            colname.B,
                            round.digits = 2,
                            group.columns = c("stim_level",
                                              "time")){
  library(dplyr)
  library(rlist)
  
  mean_own=function(x){
    sum(x)/length(x)
  }
  
  var_own=function(x){
    sum((x-mean(x))^2)/length(x)
  }
  
  apply_paste <- function(single_row){
    paste(single_row, collapse = "_")
  }
  
  data$unique_indicator <- apply(data[, group.columns], MARGIN = 1, FUN = apply_paste)
  
  D=list()
  for(single_group in unique(data$unique_indicator)){
      X_S=data %>% 
        dplyr::filter(unique_indicator == single_group)
      X_S$AB <- (X_S[[colname.A]] + X_S[[colname.B]])/2
      X_S$VAB <- ((X_S[[colname.A]] - X_S$AB)^2 + (X_S[[colname.B]] - X_S$AB)^2)/2
      
      VV_S = c(var_own(c(X_S[[colname.A]], X_S[[colname.B]])), var_own(X_S$AB), mean_own(X_S$VAB))
      CV_S = round(VV_S/VV_S[1], digits = round.digits)
      unique_characteristics <- strsplit(single_group, "_")[[1]]
      D <- list.append(D, c(unique_characteristics, VV_S, CV_S))
  }
  
  decomposition <- data.frame(matrix(unlist(D), nrow=length(unique(data$unique_indicator)), byrow=TRUE ))
  variance.colnames <- c("Tvar","Ivar","Nvar", "Tcv","Icv","Ncv")
  colnames(decomposition) <- c(group.columns, variance.colnames)
  decomposition[, variance.colnames] <- apply(decomposition[variance.colnames], 
                                              MARGIN = 2,
                                               function(x){ as.numeric(as.character(x)) } )
  return(decomposition)
}


sd.AB <- function(A, B){
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

boot.AB <- function(data,
                    n,
                    boot.no,
                    colname.A = "A",
                    colname.B = "B",
                    estimator,
                    group.columns){
  
  
  decom.all <- data.frame()
  for(i in 1:boot.no){
    if(n <= 1){
      
      bio.sample <- data %>% 
        group_by_at(group.columns) %>%
        sample_n(size = length(A)*n, replace = TRUE)
    } else {
      
      bio.sample <- data %>% 
        group_by_at(group.columns) %>%
        sample_n(size = n, replace = TRUE)
    }
    decom.single <- noise.decompose(data = bio.sample,
                                    colname.A = colname.A,
                                    colname.B = colname.B,
                                    group.columns = group.columns)
    decom.single$iteration <- i
    decom.all <- rbind(decom.all, decom.single)
    
  }
  
  if(estimator == "CI"){
    CI.data <- decom.all %>% 
      group_by_at(group.columns) %>%
      summarise(CI_5_Icv = round(quantile(x = Icv, prob = 0.05), 2),
                CI_95_Icv = round(quantile(x = Icv, prob = 0.95), 2),
                CI_5_Ncv = round(quantile(x = Ncv, prob = 0.05), 2),
                CI_95_Ncv = round(quantile(x = Ncv, prob = 0.95), 2))
  } else if(estimator == "sd"){
    CI.data <- decom.all %>% 
      group_by_at(group.columns) %>%
      summarise(sd_Icv = round(sd(Icv), 2),
                sd_Ncv = round(sd(Ncv), 2),
                mean_Icv = round(mean(Icv), 2),
                mean_Ncv = round(mean(Ncv), 2))
    
    
  }
  return(ungroup(CI.data))
}

noise.decompose.boot <- function(data,
                                 n,
                                 boot.no,
                                 colname.A = "A",
                                 colname.B = "B",
                                 estimator,
                                 group.columns = c("stim_level", "time")){

  art.boot <- boot.AB(data = data,
                      boot.no = boot.no,
                      n = n,
                      colname.A = colname.A,
                      colname.B = colname.B,
                      estimator = estimator,
                      group.columns = group.columns)
  
  bio.decompose <- noise.decompose(data = data,
                                   colname.A = colname.A,
                                   colname.B = colname.B,
                                   group.columns = group.columns)
  
  return(left_join(bio.decompose, art.boot, group.columns))
  
}

