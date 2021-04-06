# A script for reproducing the Figure S5 of a paper:

# "Phenotypic variability and not noise accounts for most of the cell-to-cell
# heterogeneity in cytokine signaling", Topolewski et al. 2021

#### set-up of directories ####
base.path <- "D:/Piotrek/publications/syncytia_noise/phenotypic_variability"
input.path <- paste(base.path, "/data", sep = "")

#### loading libraries and functions ####
source(paste(base.path, "/Topolewski_auxillary_functions.R", sep = ""))
# all code noted as "aux code" is first introduced in the above auxillary file

#### data loading ####
bio.AB <- read.csv(paste(input.path,
                         "/fused_paper.csv", sep = ""), 
                   header = TRUE, sep = ",")
merge.chosen <- "no"
#### plotting the Figure S5 ####
plots <- list()
chosen.stimulants <- unique(bio.AB$stim_type)
for(stim_type in chosen.stimulants){
  
  if(stim_type == "IFNG"){
    limits <- c(0, 60)
    x.text <- limits[1] + (limits[2] - limits[1])/30
    y.text <- limits[1] + (limits[2] - limits[1])*5/6
    breaks <- seq(limits[1], limits[2]-1, 20)
    my.pal <- IFN.pal # aux code
    time.chosen <- 15
    ggtitle <- bquote("IFN-"*gamma ~ .(time.chosen)~"min")
  } else if(stim_type == "OSM") { 
    limits <- c(0, 40)
    x.text <- limits[1] + (limits[2] - limits[1])/30
    y.text <- limits[1] + (limits[2] - limits[1])*5/6
    breaks <- seq(limits[1], limits[2]-1, 15)
    my.pal <- OSM.pal # aux code
    time.chosen <- 30
    ggtitle <- bquote("OSM" ~ .(time.chosen)~"ng/mL")
  } 
  
  bio.plot <- bio.AB[bio.AB$merged_nuclei == merge.chosen &
                       bio.AB$stim_type == stim_type, ]
  
  count.experiment <- bio.plot %>% 
    group_by(time, stim_level) %>%
    summarise(n = length(replicate)) 
  
  size.text <- 4
  alpha <- 1
  size <- 1.5
  shape <- 16
  
  plots[[paste(stim_type, "_distinct", sep = "")]] <- ggplot()+
    geom_point(data = bio.plot,
               aes(x = A,
                   y = B,
                   color = factor(replicate)),
               alpha = alpha,
               size = size,
               shape = shape)+
    geom_text(data = count.experiment,
              x = x.text,
              y = y.text,
              aes(label = paste("N = ", n, "     ", sep = "")),
              hjust = 0,
              size = size.text)+
    geom_abline(slope = 1, linetype = 2, color = "grey40")+
    theme_trajectories(aspect.ratio = 1)+
    theme(axis.text=element_text(size=6))+
    facet_grid(.~stim_level)+
    scale_color_manual(values = c("black", "grey60"),
                       guide = "none")+
    scale_x_continuous(expand = c(0, 0),
                       breaks = breaks, 
                       limits = limits,
                       name = "nucleus A response [a.u.]")+
    scale_y_continuous(expand = c(0, 0),
                       breaks = breaks, 
                       limits = limits,
                       name = "nucleus B response [a.u.]")
}

arranged.plots <- arrangeGrob(grobs = plots, nrow = 2, ncol = 1)

path.to.save <- paste(base.path, "/figures/", sep = "")
push.dir(path.to.save)
cairo_pdf(paste(path.to.save,
                "/Figure_S5.pdf", sep = ""),
          width = 7.5, height = 5)
grid.arrange(arranged.plots)
dev.off()