# A script for reproducing the Figure 2 

#### set-up of directories ####
# either provide a path to the main directory, 
# or open the "phenotypic_variability" Rproject
base.path <- getwd()
input.path <- paste(base.path, "/data", sep = "")

#### loading libraries and functions ####
source(paste(base.path, "/Topolewski_auxiliary_functions.R", sep = ""))
# all code noted as "aux code" is first introduced in the above auxiliary file

#### data loading ####
merge.chosen <- "no"

bio.AB <- read.csv(paste(input.path,
                         "/fused_paper.csv", sep = ""), 
                   header = TRUE, sep = ",")

#### plotting Figure 2 ####
chosen.stimulants.AB <- c("IFNG", "OSM")
plots <- list()
set.seed(12)
for(stim.type.chosen in chosen.stimulants.AB){
  
  if(stim.type.chosen == "IFNG"){
    limits <- c(0, 60)
    x.text <- limits[1] + (limits[2] - limits[1])/30
    y.text <- limits[1] + (limits[2] - limits[1])*5/6
    breaks <- seq(limits[1], limits[2]-1, 20)
    my.pal <- IFN.pal # aux code
    time.chosen <- 15
    ggtitle <- bquote("IFN-"*gamma ~ .(time.chosen)~"min")
  } else if(stim.type.chosen == "OSM") {
    limits <- c(0, 40)
    x.text <- limits[1] + (limits[2] - limits[1])/30
    y.text <- limits[1] + (limits[2] - limits[1])*5/6
    breaks <- seq(limits[1], limits[2]-1, 15)
    my.pal <- OSM.pal # aux code
    time.chosen <- 30
    ggtitle <- bquote("OSM" ~ .(time.chosen)~"ng/mL")
  } 
  bio.plot <- bio.AB %>%
    dplyr::filter(merged_nuclei == merge.chosen &
                    stim_type == stim.type.chosen &
                    time == time.chosen)
    
  
  decom.fold <- noise.decompose.boot(data = bio.plot,
                                     n = 1,
                                     boot.no = 1000,
                                     estimator = "sd",
                                     colname.A = "A",
                                     colname.B = "B",
                                     group.columns = c("stim_level",
                                                       "time")) # aux code
  
  alpha <- 1
  size <- 1.5
  shape <- 16
  
  
  plots[[paste(stim.type.chosen, "_scatter", sep = "")]] <- ggplot()+
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
    theme(axis.text=element_text(size=6),
          legend.position = "bottom")+
    facet_grid(.~stim_level)+
    scale_color_manual(values = my.pal,
                       guide = "none")+
    scale_x_continuous(expand = c(0, 0),
                       breaks = breaks,
                       limits = limits,
                       name = "nucleus A response [a.u.]")+
    scale_shape_manual(values = c(16, 6),
                       guide = NULL)+
    scale_y_continuous(expand = c(0, 0),
                       breaks = breaks,
                       limits = limits,
                       name = "nucleus B response [a.u.]")
  plots[[paste(stim.type.chosen, "_scatter", sep = "")]] 
}


arranged.plots <- arrangeGrob(grobs = plots, layout_matrix = rbind(c(1, 1, 1),
                                                                   c(2, 2, 2)))

width <- 7.5
height <- 7
filename_suffix <- "Figure_2.pdf"


path.to.save <- paste(base.path, "/figures", sep = "")
push.dir(path.to.save) # aux code

cairo_pdf(paste(path.to.save, "/",
                filename_suffix, sep = ""),
          width = width, height = height)
grid.arrange(arranged.plots)
dev.off()

