# A script for reproducing the Figure ___ of a paper:

# "Phenotypic variability and not noise accounts for most of the cell-to-cell
# heterogeneity in cytokine signaling", Topolewski et al. 2021

#### loading libraries and functions ####
source("D:/Piotrek/scripts/basic_scripts/normalize_data.R")
source("D:/Piotrek/scripts/identity_analysis/pie_density.R")

#### set-up of directories ####
base.path <- "D:/Piotrek/Experiments/ICF/"
folder.name <- "cellular_identity"
project <- "trajectories"
normaliz.chosen <- "pbs"
path <- paste(base.path, folder.name, sep='')

#### data loading ####
bio.traj <- read.csv(paste(path, "/", project, "/input/",
                          normaliz.chosen, "/",
                          project,
                          "_papeR.csv", sep = ""))
bio.traj <- bio.traj %>%
  dplyr::filter(nuclearity == "single-nuclear")

#### choosing the stimulants ####
chosen.stimulants <- c("IFNG", "OSM")

#### plotting the concentration trajectory with piecharts ####
plots <- list()
letter.i <- 0
for(stim.type.chosen in chosen.stimulants){
  
  for(time.point in 1:5){
    letter.i <- letter.i + 1
    #### plotting the concentration trajectories ####
    
    ## providing details about the repetitions ##
    if(stim.type.chosen == "IFNG"){
      x.limits <- c(0, 17)
      y.limits <- c(0, 1.1)
      my.pal <- IFN.pal
      all.time.points <- c(5, 15, 30, 45, 90)
      time.chosen <- all.time.points[time.point]
      
      x.text <- x.limits[1] + (x.limits[2] - x.limits[1])/35
      y.text <- y.limits[1] + (y.limits[2] - y.limits[1])*4/6*0.9
      
      x.breaks <- seq(x.limits[1], x.limits[2], 6)
      y.breaks <- seq(y.limits[1], y.limits[2], 0.5)
      
      ggtitle <- bquote("IFN-"*gamma ~ .(time.chosen)~"min")
    } else if(stim.type.chosen == "OSM") { 
      x.limits <- c(0, 25)
      y.limits <- c(0, 1.3)
      my.pal <- OSM.pal
      all.time.points <- c(5, 15, 30, 60, 90)
      time.chosen <- all.time.points[time.point]
      x.text <- x.limits[1] + (x.limits[2] - x.limits[1])/50
      y.text <- y.limits[1] + (y.limits[2] - y.limits[1])*3/6
      
      x.breaks <- seq(x.limits[1], x.limits[2], 10)
      y.breaks <- seq(y.limits[1], y.limits[2], 0.6)
      ggtitle <- bquote("OSM" ~ .(time.chosen)~"min")
    } 
    ## subsetting the data accordingly ##
    bio.plot <- bio.traj[bio.traj$time == time.chosen &
                           bio.traj$stim_type == stim.type.chosen, ]
    letter = LETTERS[letter.i] # choosing the letter to annotate the figure 
    

    plots[[paste(stim.type.chosen, time.point)]] <- ggplot()+
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
      add_letter(letter)+
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
    
  
  
  #### plotting the pie charts ####
    
    ## providing replicate details ##
  if(stim.type.chosen == "IFNG"){
    my.pal <- IFN.pal
    paste("dose", stim.type.chosen, " [ng/mL]")
    xtitle <- bquote("IFN-"*gamma ~ "dose [ng/mL]")
  } else if(stim.type.chosen == "OSM") { 
    my.pal <- OSM.pal
    xtitle <- bquote("OSM dose [ng/mL]")
  }
    ## calculating a confusion matrix ## 
  bio.CM <- confusion.matrix.bio(data = bio.traj,
                                 stim_type_chosen = stim.type.chosen,
                                 time.chosen = time.chosen,
                                 matrix = FALSE)
  
  if(0 %in% bio.CM$prob){
    # for clarity on the plot, all probabilities = 0 are set to NA
    bio.CM[bio.CM$prob == 0, ]$prob <- NA 
  }
  
  
  plots[[paste(stim.type.chosen, time.point, "CM")]] <- ggplot(bio.CM,
                                                               aes(fill = factor(push.stim),
                                                                   x = 1,
                                                                   y = prob))+
    geom_bar(stat = "identity",
             color = "black") +
    facet_grid(stim~push.stim, switch = "y")+
    coord_polar(theta = "y", start = 0) +
    scale_y_continuous(breaks = NULL,
                       position = "right",
                       limits = c(0, 1),
                       name = "dose for which response \nis typical [ng/mL]")+
    scale_x_continuous(breaks = NULL,
                       name = xtitle)+
    add_letter(" ", color = "white")+
    geom_vline(xintercept = 1.45,
               colour = "black") +
    scale_fill_manual(values = my.pal,
                      guide = "none")+
    theme_trajectories(axis.name.size = 11)+
    theme(panel.border = element_blank(),
          panel.spacing.x=unit(0, "mm"),
          panel.spacing.y=unit(0, "mm"))
  plots[[paste(stim.type.chosen, time.point, "CM")]]
  }
}

arranged.plots <- arrangeGrob(grobs = plots,
                              layout_matrix = rbind(seq(1, 9, by = 2), 
                                                    seq(2, 10, by = 2), 
                                                    seq(11, 19, by = 2),
                                                    seq(12, 20, by = 2)))
width <- 10/4*5
height <- 8.5/3*4.5
filename_suffix <- "Fig1_allTimes.pdf"

## plot saving ##
path.to.save <- paste(working.path, "/output/", normaliz.chosen, sep = "")
push.dir(path.to.save)

pdf(paste(path.to.save, "/", 
          filename_suffix, sep = ""), 
    width = width, height = height)
grid.arrange(arranged.plots)
dev.off()

