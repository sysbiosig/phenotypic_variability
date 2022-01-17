#### A script with auxiliary functions for a paper: ####

# "Phenotypic variability and not noise accounts for most of the cell-to-cell
# heterogeneity in cytokine signaling", Topolewski et al. 2021

# first install the SLEMI package according to https://github.com/sysbiosig/SLEMI #

package.list <- list("ggplot2", 
                     "gridExtra",
                     "data.table",
                     "reshape2",
                     "dplyr",
                     "scales",
                     "RColorBrewer",
                     "rlist",
                     "SLEMI")
try({package.list
  package.load <- sapply(package.list, function(package.name){
    package.exist <- require(package.name, character.only = TRUE)
    if(!package.exist){
      install.packages(package.name)
      return(library(package.name, character.only = TRUE))
    }
    return(package.exist)
  })
})

#### ggplot theme ####

theme_trajectories <- function(border.thickness = 0.5,
                               axis.num.size = 10, 
                               axis.name.size = 11,
                               aspect.ratio = FALSE){
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size=border.thickness, fill=NA),
        panel.background = element_blank(),
        axis.ticks = element_line(colour = "black", size = border.thickness),
        axis.text.x = element_text(colour = "black", size = axis.num.size),
        axis.text.y = element_text(colour = "black", size = axis.num.size),
        strip.background = element_blank(),
        axis.title = element_text(face="plain", size = axis.name.size),
        plot.title = element_text(hjust = 0.5),
        legend.key=element_blank())+
    if(aspect.ratio != FALSE){
      theme(aspect.ratio = aspect.ratio)
    } else {theme()}
}

theme.itrc <- function(){
  theme_itrc <-theme(axis.title = element_text(size= 11,
                                               face="plain",
                                               vjust = 0.5),
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     axis.line.x = element_blank(),
                     axis.line.y = element_blank(),
                     panel.background = element_blank(),
                     strip.background = element_blank(),
                     strip.text =
                       element_text(size= 10,
                                    face="plain",
                                    vjust = 0.5#,lineheight = theme.text_size*3
                       ))
}

add_letter <- function(letter = "A", size = 20, face = "bold", color = "black"){
  list(labs(tag = letter),
       theme(plot.tag = element_text(size = size, face = face, color = color)))
} 

#### palettes for IFNG and OSM ####
IFN.pal <- c(brewer.pal(8, "PuBu"))[c(3,5,7,8)]
OSM.pal <- c("#cbdeda", "#63bda3", "#218845", "#114320")

#### function to push the creation of a new folder if not yet existing ####
push.dir <- function(folder.name){
  if(!dir.exists(folder.name)){
    dir.create(folder.name, recursive = TRUE)
  }
}


#### noise_decompose ####
noise.decompose <- function(data,
                            colname.A,
                            colname.B,
                            round.digits = 2,
                            group.columns = c("stim_level",
                                              "time")){
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

boot.AB <- function(data,
                    n,
                    boot.no,
                    colname.A = "A",
                    colname.B = "B",
                    estimator,
                    group.columns,
                    round.digits = 2){
  
  
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
      summarise(CI_5_Icv = round(quantile(x = Icv, prob = 0.05), round.digits),
                CI_95_Icv = round(quantile(x = Icv, prob = 0.95), round.digits),
                CI_5_Ncv = round(quantile(x = Ncv, prob = 0.05), round.digits),
                CI_95_Ncv = round(quantile(x = Ncv, prob = 0.95), round.digits))
  } else if(estimator == "sd"){
    CI.data <- decom.all %>% 
      group_by_at(group.columns) %>%
      summarise(sd_Icv = round(sd(Icv), round.digits),
                sd_Ncv = round(sd(Ncv), round.digits),
                mean_Icv = round(mean(Icv), round.digits),
                mean_Ncv = round(mean(Ncv), round.digits))
    
    
  }
  return(ungroup(CI.data))
}

noise.decompose.boot <- function(data,
                                 n,
                                 boot.no,
                                 colname.A = "A",
                                 colname.B = "B",
                                 estimator,
                                 round.digits = 2,
                                 group.columns = c("stim_level", "time")){
  
  art.boot <- boot.AB(data = data,
                      boot.no = boot.no,
                      n = n,
                      colname.A = colname.A,
                      colname.B = colname.B,
                      estimator = estimator,
                      round.digits = round.digits,
                      group.columns = group.columns)
  
  bio.decompose <- noise.decompose(data = data,
                                   colname.A = colname.A,
                                   colname.B = colname.B,
                                   group.columns = group.columns)
  
  return(left_join(bio.decompose, art.boot, group.columns))
  
}


#### pie density ####
pie.density=function(ld){
  # all densities have to be on the same grid
  xx=ld[[1]]$x # for all desnities 
  n=length(ld) # number of concentration
  ## re-writing densities from list to matrix 
  ddM = matrix(0,length(xx),n)
  ddI = matrix(FALSE,length(xx),n)
  for(i in 1:n){
    ddM[,i]=ld[[i]]$y
  }
  
  MC=max.col(ddM) #coumn in which max is achieved
  for(i in 1:n){
    ddI[,i]=(MC==i) # True is max is column i
  }
  
  dx=diff(xx)[1] #dx for integration
  CM=matrix(0,n,n)
  for(i in 1:n) # trick
    for(j in 1:n)
      CM[i,j]=sum(ddM[ddI[,j],i])*dx
  return(CM)
}

confusion.matrix.bio <- function(data,
                                 time.chosen,
                                 xmax= 40,
                                 stim_type_chosen,
                                 matrix = FALSE){
  bio.pie <- data %>% dplyr::filter(stim_type == stim_type_chosen)
  D=list() # list of densities correspoding to every stimulation level
  i=1
  Uss <- sort(unique(bio.pie$stim_level))
  for (ss in sort(Uss)){# for every stimulation level ss
    X_S=bio.pie %>%
      dplyr::filter(stim_level == ss,
                    time == time.chosen)
    dens=density(X_S$response,from=0, to=xmax)# density estimation
    D[[i]] <- dens;
    i=i+1
  }
  CM=pie.density(D)
  colnames(CM) <- Uss
  rownames(CM) <- Uss
  if(matrix == TRUE){
    return(CM)
  } else{
    CM.pt <- melt(CM)
    colnames(CM.pt) <- c("stim", "push.stim", "prob")
    return(CM.pt)
  }
}


pie_gamma <- function(ld,xmax){
  
  xx=seq(0.01,xmax,by=0.01)
  n=nrow(ld)
  dd=list()
  ddM=matrix(0,length(xx),n)
  ddI=matrix(FALSE,length(xx),n)
  for(i in 1:n){
    dd[[i]]=dgamma(xx, shape =ld[i,1], rate = ld[i,2])
    ddM[,i]=dd[[i]]
  }
  
  MC=max.col(ddM)
  for(i in 1:n){
    ddI[,i]=(MC==i)
  }
  
  
  dx=diff(xx)[1]
  CM=matrix(0,n,n)
  for(i in 1:n)
    for(j in 1:n)
      CM[i,j]=sum(ddM[ddI[,j],i])*dx
  
  CM[CM>1]=1
  CM[CM<10^-8]=0
  CM
}



confusion.matrix.gamma <- function(data,
                                   time.chosen,
                                   xmax= 40,
                                   stim_type_chosen,
                                   sd.coefficient,
                                   percentiles = seq(0.05,0.95,by=0.01)){
  
  
  D=list() # D is a list of densities correspoding to every stimulation level
  i=1
  
  Uss <- sort(unique(data$stim_level))
  for (ss in sort(Uss)){# for every stimulation level ss
    X_S=data %>% 
      dplyr::filter(stim_level == ss,
                    time == time.chosen,
                    stim_type == stim_type_chosen) 
    dens=density(X_S$response,from=0, to=xmax)# density estimation
    D[[i]] <-dens;

    D[[i]]$q <-quantile(X_S$response,percentiles);
    D[[i]]$m <-mean(X_S$response);
    D[[i]]$v <-var(X_S$response);
    i=i+1
  }
  ldlist=list() # parameters of gamma for different percentiles
  for(k in 1:length(percentiles)){
    l=1
    ld = matrix(0,4,2)
    ld[l,]=c(1/sd.coefficient^2,1 / (sd.coefficient^2*D[[l]]$q[k]))
    l=2
    ld[l,]=c(1/sd.coefficient^2,1 / (sd.coefficient^2*D[[l]]$q[k]))
    l=3
    ld[l,]=c(1/sd.coefficient^2,1 / (sd.coefficient^2*D[[l]]$q[k]))
    l=4
    ld[l,]=c(1/sd.coefficient^2,1 / (sd.coefficient^2*D[[l]]$q[k]))
    ldlist[[k]]=ld
  }
  
  

    CMloop=matrix(0,4,4)
    for(k in 1:length(percentiles)){
      CMloop = CMloop+pie_gamma(ldlist[[k]],xmax)
      
      
    }

  
  CM=CMloop/length(percentiles) 
  colnames(CM) <- Uss
  rownames(CM) <- Uss
  CM.pt <- melt(CM)
  colnames(CM.pt) <- c("stim", "push.stim", "prob")
  return(CM.pt)
  
}
