# A script for reproducing the Figure S4

#### set-up of directories ####
# either provide a path to the main directory, 
# or open the "phenotypic_variability" Rproject
base.path <- getwd()
input.path <- paste(base.path, "/data", sep = "")

#### loading libraries and functions ####
source(paste(base.path, "/Topolewski_auxiliary_functions.R", sep = ""))
# all code noted as "aux code" is first introduced in the above auxiliary file

#### data loading ####
bio.AB <- read.csv(paste(input.path,
                         "/fused_paper.csv", sep = ""), 
                   header = TRUE, sep = ",")

bio.AB <- bio.AB %>%
  dplyr::filter(!((time %in% c(15) & stim_type == "IFNG" & stim_level %in% c(0.1, 1)) |
                    (time %in% c(30) & stim_type == "OSM" & stim_level %in% c(0.1, 1))))

# to present 0 on the facet as time point, annotate the data accordingly:
bio.AB[bio.AB$stim_level == 0, "time"] <- 0
bio.AB[bio.AB$stim_level == 0, "stim_level"] <- 10

merge.chosen <- "no"

#### grid ####
set.seed(122)
plots <- list()
for(stim_type_chosen in c("IFNG", "OSM")){
  
  bio.plot <- bio.AB[bio.AB$merged_nuclei == merge.chosen &
                       bio.AB$stim_type == stim_type_chosen, ]
  
  if(stim_type_chosen == "IFNG"){
    limits <- c(0, 60)
    x.text <- limits[1] + (limits[2] - limits[1])/30
    y.text <- limits[1] + (limits[2] - limits[1])*5/6
    breaks <- seq(limits[1], limits[2]-1, 20)
    my.pal <- IFN.pal[4]
    
    # to keep the facet size consistent, create an artificial facet, 
    # which will be removed in further processing of a graph
    art <- bio.plot[1,]
    art$time <- NA
    art$A <- NA
    art$B <- NA
    bio.plot <- rbind(bio.plot,
                      art)
    
  } else if(stim_type_chosen == "OSM") {
    limits <- c(0, 40)
    x.text <- limits[1] + (limits[2] - limits[1])/30
    y.text <- limits[1] + (limits[2] - limits[1])*5/6
    breaks <- seq(limits[1], limits[2]-1, 15)
    my.pal <- OSM.pal[4]
  }
  
  decom.fold <- noise.decompose.boot(data = bio.plot,
                                     n = 1,
                                     round.digits = 2,
                                     boot.no = 1000,
                                     estimator = "sd",
                                     colname.A = "A",
                                     colname.B = "B",
                                     group.columns = c("stim_level", 
                                                       "time"))
  
  decom.fold$time <- as.numeric(as.character(decom.fold$time))
  decom.fold$stim_level <- as.numeric(as.character(decom.fold$stim_level))
  
  alpha <- 1
  size <- 1.5
  shape <- 16
  size.text <- 4
  
  
  plots[[paste(stim_type_chosen, "_scatter", sep = "")]] <- ggplot()+
    geom_point(data = bio.plot,
               aes(x = A,
                   y = B,
                   color = factor(stim_level),
                   shape = replicate),
               alpha = alpha,
               size = size)+
    geom_text(data = decom.fold,
              parse = TRUE,
              x = x.text,
              y = y.text,
              aes(label = paste("atop(''* rho[italic(phen.)] * ''== ", Icv,
                                ",''* rho[italic(noise)] * ''*''== ", Ncv, ")",
                                sep = "")),
              hjust = 0,
              size = 4)+
    geom_text(data = decom.fold,
              parse = TRUE,
              x = x.text+limits[2]*0.45,
              y = y.text,
              aes(label = paste("phantom(0)%+-%", sd_Ncv, sep = "")),
              hjust = 0,
              size = 3.25)+
    geom_abline(slope = 1, linetype = 2, color = "grey40")+
    theme_trajectories(aspect.ratio = 1)+
    theme(axis.text=element_text(size=6))+
    facet_grid(stim_level~time)+
    scale_shape_manual(values = c(16, 6),
                       guide = NULL)+
    scale_color_manual(values = my.pal,
                       guide = "none")+
    scale_x_continuous(expand = c(0, 0),
                       breaks = breaks,
                       limits = limits,
                       name = "nucleus A response [a.u.]")+
    scale_y_continuous(expand = c(0, 0),
                       breaks = breaks,
                       limits = limits,
                       name = "nucleus B response [a.u.]")
  plots[[paste(stim_type_chosen, "_scatter", sep = "")]]
}

arranged.plots <- arrangeGrob(grobs = plots, layout_matrix = rbind(c(1),
                                                                   c(2)))

path.to.save <- paste(base.path, "/figures", sep = "")
push.dir(path.to.save) # aux code


width <- 7.5
height <- 6
filename_suffix <- "Figure_S5.pdf"

cairo_pdf(paste(path.to.save, "/",
                filename_suffix, sep = ""),
          width = width, height = height)
grid.arrange(arranged.plots)
dev.off()


