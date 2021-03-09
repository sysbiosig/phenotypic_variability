# A script for reproducing the Figure ___ of a paper:

# "Phenotypic variability and not noise accounts for most of the cell-to-cell
# heterogeneity in cytokine signaling", Topolewski et al. 2021

#### loading libraries and functions ####
source("D:/Piotrek/scripts/basic_scripts/normalize_data.R")
source("D:/Piotrek/scripts/identity_analysis/pie_density.R")
library(SLEMI)

#### set-up of directories ####
base.path <- "D:/Piotrek/Experiments/ICF/"
folder.name <- "cellular_identity"
project <- "trajectories"
normaliz.chosen <- "pbs"
path <- paste(base.path, folder.name, sep='')

#### data loading ####
bio.traj <- read.csv(paste(path, "/", project, "/input/", 
                           normaliz.chosen,
                           "/trajectories_paper.csv", sep = ""),  ## if you want to plot TNF, remove "_paper"
                     header = TRUE, sep = ",")

bio.traj <- bio.traj %>%
  dplyr::filter(merged_nuclei %in% c(NA))

#### choosing the stimulants ####
chosen.stimulants <- c("IFNG", "OSM")
# chosen.stimulants <- c("TNF")

#### plotting the time trajectory for single-nuclear cells ####
plots <- list()
for(stim.type.chosen in chosen.stimulants){
  
  ## providing replicate specific details for plotting ##
  if(stim.type.chosen == "IFNG"){
    x.limits <- c(0, 17)
    y.limits <- c(0, 1.1)
    my.pal <- IFN.pal
    time.chosen <- 15
    
    x.text <- x.limits[1] + (x.limits[2] - x.limits[1])/35
    y.text <- y.limits[1] + (y.limits[2] - y.limits[1])*4/6*0.9
    
    x.breaks <- seq(x.limits[1], x.limits[2], 6)
    y.breaks <- seq(y.limits[1], y.limits[2], 0.5)
    
    ggtitle <- bquote("IFN-"*gamma ~ .(time.chosen)~"min")
  } else if(stim.type.chosen == "OSM") { 
    x.limits <- c(0, 25)
    y.limits <- c(0, 1.3)
    my.pal <- OSM.pal
    time.chosen <- 30
    x.text <- x.limits[1] + (x.limits[2] - x.limits[1])/50
    y.text <- y.limits[1] + (y.limits[2] - y.limits[1])*3/6
    
    x.breaks <- seq(x.limits[1], x.limits[2], 10)
    y.breaks <- seq(y.limits[1], y.limits[2], 0.6)
    ggtitle <- bquote("OSM" ~ .(time.chosen)~"min")
  } else if(stim.type.chosen == "TNF"){
    x.limits <- c(0, 8)
    y.limits <- c(0, 1.1)
    my.pal <- TNF.pal
    time.chosen <- 30
    
    x.text <- x.limits[1] + (x.limits[2] - x.limits[1])/12
    y.text <- y.limits[1] + (y.limits[2] - y.limits[1])*4/6*0.9
    
    x.breaks <- seq(x.limits[1], x.limits[2], 6)
    y.breaks <- seq(y.limits[1], y.limits[2], 0.5)
    
    ggtitle <- bquote("TNF-"*alpha ~ .(time.chosen)~"min")
  } else if(stim.type.chosen == "TNF_h"){
    x.limits <- c(0, 8)
    y.limits <- c(0, 1.1)
    my.pal <- TNF.pal
    time.chosen <- 30
    
    x.text <- x.limits[1] + (x.limits[2] - x.limits[1])/12
    y.text <- y.limits[1] + (y.limits[2] - y.limits[1])*4/6*0.9
    
    x.breaks <- seq(x.limits[1], x.limits[2], 6)
    y.breaks <- seq(y.limits[1], y.limits[2], 0.5)
    
    ggtitle <- bquote("TNF-h"*alpha ~ .(time.chosen)~"min")
  }
  
  bio.plot <- bio.traj %>%
    dplyr::filter(time == time.chosen &
                         stim_type == stim.type.chosen)

  
  plots[[stim.type.chosen]] <- ggplot()+
    geom_line(stat="density",
              data = bio.plot,
              aes(x = (response),
                  group = stim_level,
                  color = factor(stim_level)),
              lwd=1,
              position = "identity")+
    scale_color_manual(values = my.pal,
                       name = "dose [ng/mL]")+
    scale_fill_manual(values = my.pal,
                      guide = FALSE)+
    theme_trajectories(aspect.ratio = 1)+
    theme(legend.position = c(x.text, y.text),
          legend.key.size = unit(4, "mm"))+
    coord_cartesian(xlim = x.limits,
                    ylim = y.limits)+
    scale_x_continuous(expand = c(0, 0),
                       breaks = x.breaks,
                       name = "relative fluorescence intensity")+
    scale_y_continuous(expand = c(0, 0),
                       breaks = y.breaks,
                       name = "probabilty density")+
    ggtitle(ggtitle)
  plots[[stim.type.chosen]]
  
}
#### concentraion trajectory in wide version, Fig. 1C ####

## providing replicate specific details for plotting ##
x.limits <- c(0, 17)
y.limits <- c(0, 1.1)
my.pal <- IFN.pal
time.chosen <- 15
stim.type.chosen <- "IFNG"
x.text <- x.limits[1] + (x.limits[2] - x.limits[1])/50
y.text <- y.limits[1] + (y.limits[2] - y.limits[1])*4/6

x.breaks <- seq(x.limits[1], x.limits[2], 6)
y.breaks <- seq(y.limits[1], y.limits[2], 0.5)

ggtitle <- bquote("IFN-"*gamma ~ .(time.chosen)~"min")

bio.plot <- bio.traj %>%
  dplyr::filter(time == time.chosen &
                  stim_type == stim.type.chosen)

if(stim.type.chosen != "TNF"){
plots[["long"]] <- ggplot()+
  geom_line(stat="density",
            data = bio.plot,
            aes(x = (response),
                group = stim_level,
                color = factor(stim_level)),
            lwd=1,
            position = "identity",
            kernel = "e")+
  geom_density(data = bio.plot[bio.plot$stim_level == 1,],
               fill = "red",
               color = NA,
               aes(x = response),
               alpha = 1)+
  scale_color_manual(values = my.pal,
                     guide = FALSE)+
  scale_fill_manual(values = my.pal,
                    guide = FALSE)+
  theme_trajectories(aspect.ratio = 1/2.5)+
  theme(legend.position = c(x.text, y.text))+
  coord_cartesian(xlim = x.limits,
                  ylim = y.limits)+
  scale_x_continuous(expand = c(0, 0),
                     breaks = x.breaks,
                     name = "relative fluorescence intensity")+
  scale_y_continuous(expand = c(0, 0),
                     breaks = y.breaks,
                     name = "probabilty density")+
  ggtitle(ggtitle)
# plots[["long"]]
}

#### plotting the pie charts for single-nuclear cells ####
for(stim.type.chosen in chosen.stimulants){
  
  if(stim.type.chosen == "IFNG"){
    my.pal <- IFN.pal
    time.chosen <- 15
    paste("dose", stim.type.chosen, " [ng/mL]")
    xtitle <- bquote("IFN-"*gamma ~ "dose [ng/mL]")
  } else if(stim.type.chosen == "OSM") { 
    my.pal <- OSM.pal
    time.chosen <- 30
    xtitle <- bquote("OSM dose [ng/mL]")
  } else if(stim.type.chosen == "TNF"){
    my.pal <- TNF.pal
    time.chosen <- 30
    paste("dose", stim.type.chosen, " [ng/mL]")
    xtitle <- bquote("TNF-"*alpha ~ "dose [ng/mL]")
  } else if(stim.type.chosen == "TNF_h"){
    my.pal <- TNF.pal
    time.chosen <- 30
    paste("dose", stim.type.chosen, " [ng/mL]")
    xtitle <- bquote("TNF_h-"*alpha ~ "dose [ng/mL]")
  }
  bio.CM <- confusion.matrix.bio(data = bio.traj,
                                 stim_type_chosen = stim.type.chosen,
                                 time.chosen = time.chosen,
                                 matrix = FALSE)
  
  bio.capacity <- bio.traj %>%
    dplyr::filter(time == time.chosen &
                    stim_type == stim.type.chosen)
  
  capacity <- SLEMI::capacity_logreg_main(dataRaw = bio.capacity, 
                                          signal = "stim_level", 
                                          response = "response")$cc
  
  if(0 %in% bio.CM$prob){
    # for clarity on the plot, all probabilities = 0 are set to NA
    bio.CM[bio.CM$prob == 0, ]$prob <- NA 
}
  plots[[paste(stim.type.chosen, "CM")]] <- ggplot()+
    geom_bar(data = bio.CM,
              aes(fill = factor(push.stim),
                  x = 1,
                  y = prob),
              stat = "identity",
             color = "black") +
    labs(caption = bquote(C*"*" * "=" * .(round(capacity, 2))))+
    facet_grid(stim~push.stim, switch = "y")+
    coord_polar(theta = "y", start = 0) +
    scale_y_continuous(breaks = NULL,
                       position = "right",
                       limits = c(0, 1),
                       name = "dose for which response \nis typical [ng/mL]")+
    scale_x_continuous(breaks = NULL,
                       name = xtitle)+
    geom_vline(xintercept = 1.45,
               colour = "black") +
    scale_fill_manual(values = my.pal,
                      guide = "none")+
    theme_trajectories(axis.name.size = 11)+
    theme(panel.border = element_blank(),
          panel.spacing.x=unit(0, "mm"),
          panel.spacing.y=unit(0, "mm"),
          plot.caption = element_text(size = 25,
                                      hjust = 0.5))
  plots[[paste(stim.type.chosen, "CM")]]
  
}
arranged.plots <- arrangeGrob(grobs = plots, nrow = 3, ncol = 2,
                              layout_matrix = rbind(c(1,2), c(3,3), c(4, 5)))
width <- 5
height <- 9
filename_suffix <- "Figure_1_capacity.pdf"

# arranged.plots <- arrangeGrob(grobs = plots, nrow = 3, ncol = 3, 
#                               layout_matrix = rbind(c(1,2,3), c(4,4,4), c(5, 6, 7)))
# 
# arranged.plots <- arrangeGrob(grobs = plots, nrow = 3, ncol = 4, 
#                               layout_matrix = rbind(c(1,2,3,4), c(5,5,5,5), c(6, 7, 8, 9)))

## for TNF mouse only ##
# arranged.plots <- arrangeGrob(grobs = plots[-2],
#                               nrow = 2, ncol = 1)
# width <- 10/4
# height <- 7
# filename_suffix <- "_conc_TNF_only.pdf"

## plot saving ##
path.to.save <- paste(working.path, "/output/", normaliz.chosen, sep = "")
push.dir(path.to.save)

pdf(paste(path.to.save, "/", 
          filename_suffix, sep = ""), 
    width = width, height = height)
grid.arrange(arranged.plots)
dev.off()

# jpeg(paste(path.to.save, "/", project, "_", merge.chosen,
#            "_conc_TNF2.jpeg", sep = ""),
#      width = 10*300, height = 8.5*300,
#      quality = 100,
#      res = 300,
#      type = "cairo")
# grid.arrange(arranged.plots)
# dev.off()


#### saving confusion matrixes ####
# path.to.save <- paste(working.path, "/output/", normaliz.chosen, sep = "")
# for(stim.type.chosen in c("IFNG", "OSM")){
#   
#   if(stim.type.chosen == "IFNG"){
#     time.chosen <- 15
#   } else if(stim.type.chosen == "OSM") { 
#     time.chosen <- 30
#     }
# bio.CM <- confusion.matrix.bio(data = bio.traj,
#                                stim.type.chosen = stim.type.chosen,
#                                time.chosen = time.chosen,
#                                matrix = TRUE)
# write.csv(bio.CM, paste(path.to.save, "/", stim.type.chosen, "_CM.csv", sep = ""))
# }
