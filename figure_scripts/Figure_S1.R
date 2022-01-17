# A script for reproducing the Figure S1 of a paper:

# "Phenotypic variability, not noise, accounts for most of the cell-to-cell
# heterogeneity of selected cytokine-induced JAK-STAT signaling responses", 
# Topolewski et al. 2021

#### set-up of directories ####
base.path <- "D:/Piotrek/publications/syncytia_noise/phenotypic_variability"
input.path <- paste(base.path, "/data", sep = "")

#### loading libraries and functions ####
source(paste(base.path, "/Topolewski_auxiliary_functions.R", sep = ""))
# all code noted as "aux code" is first introduced in the above auxiliary file

#### data loading ####
merge.chosen <- "no"

bio.traj <- read.csv(paste(input.path,
                           "/trajectories_paper.csv", sep = ""), 
                     header = TRUE, sep = ",")

bio.traj <- bio.traj %>%
  dplyr::filter(merged_nuclei %in% c("no", NA))

#### plotting the time trajectory for single-nuclear cells ####
plots <- list()
adjust.kernel <- 2 # adjusting the kernel density
aspect.chosen <- 0.75 # 3:4 ratio of a plot area
letter.i <- 1 # for figure annotation

for(stim.type.chosen in c("IFNG", "OSM")){
  
  ## providing replicate specific details for plotting ##
  stim.level.chosen = 10
  if(stim.type.chosen == "IFNG"){
    y.limits <- c(0, 40)
    my.pal <- IFN.pal
    y.breaks <- seq(y.limits[1], y.limits[2], 15)
    ggtitle <- bquote("IFN-"*gamma ~ .(stim.level.chosen)~"ng/mL, single-cells")
  } else if(stim.type.chosen == "OSM") {
    y.limits <- c(0, 25)
    my.pal <- OSM.pal
    y.breaks <- seq(y.limits[1], y.limits[2], 10)
    ggtitle <- bquote("OSM" ~ .(stim.level.chosen)~"ng/mL, single-cells")
  }
  
  bio.plot <- bio.traj %>%
    dplyr::filter(stim_level == stim.level.chosen &
                    stim_type == stim.type.chosen &
                    nuclearity == "single-nuclear")
  x.breaks <- unique(bio.plot$time)
  
  width <- 9 * (y.limits[2] / 25)
  print(paste(stim.type.chosen, width))
  plots[[paste("single.nuc", stim.type.chosen, sep = "")]] <- ggplot()+
    geom_violin(data = bio.plot,
                aes(x = time,
                    y = (response),
                    group = time),
                color = NA,
                trim = FALSE,
                fill = my.pal[4],
                width = width,
                scale = "area",
                adjust = adjust.kernel)+
    add_letter(letter = LETTERS[letter.i])+ # aux code
    theme_trajectories(aspect.ratio = aspect.chosen)+ # aux code
    scale_y_continuous(expand = c(0, 0),
                       breaks = y.breaks)+
    scale_x_continuous(breaks = x.breaks)+
    ylab("relative fluorescence intensity")+
    xlab("time [min]")+
    coord_cartesian(ylim = y.limits,
                    xlim = c(-2, 92))+
    ggtitle(ggtitle)
  plots[[paste("single.nuc", stim.type.chosen, sep = "")]]
  
  letter.i <- letter.i + 1
}

#### plotting the time trajectory for bi-nuclear cells ####

for(stim.type.chosen in c("IFNG", "OSM")){
  
  if(stim.type.chosen == "IFNG"){
    y.limits <- c(0, 40)
    my.pal <- IFN.pal # aux code
    y.breaks <- seq(y.limits[1], y.limits[2], 15)
    ggtitle <- bquote("IFN-"*gamma ~ .(stim.level.chosen)~"ng/mL, syncytia")
    width <- 9 * (y.limits[2] / 25)
  } else if(stim.type.chosen == "OSM") {
    y.limits <- c(0, 25)
    my.pal <- OSM.pal # aux code
    y.breaks <- seq(y.limits[1], y.limits[2], 10)
    ggtitle <- bquote("OSM" ~ .(stim.level.chosen)~"ng/mL, syncytia")
    width <- 9 * (y.limits[2] / 50)
  }
  bio.plot <- bio.traj %>%
    dplyr::filter(stim_level == stim.level.chosen &
                    stim_type == stim.type.chosen &
                    nuclearity == "bi-nuclear")
  x.breaks <- unique(bio.plot$time)
  
  
  plots[[paste("bi.nuc", stim.type.chosen, sep = "")]] <- ggplot()+
    geom_violin(data = bio.plot,
                aes(x = time,
                    y = (response),
                    group = time),
                color = NA,
                trim = FALSE,
                fill = my.pal[4],
                width = width,
                scale = "area",
                adjust = adjust.kernel)+
    add_letter(letter = LETTERS[letter.i])+ # aux code
    theme_trajectories(aspect.ratio = aspect.chosen)+ # aux code
    scale_y_continuous(expand = c(0, 0),
                       breaks = y.breaks)+
    scale_x_continuous(breaks = x.breaks)+
    ylab("relative fluorescence intensity")+
    xlab("time [min]")+
    coord_cartesian(ylim = y.limits,
                    xlim = c(-2, 92))+
    ggtitle(ggtitle)
  plots[[paste("bi.nuc", stim.type.chosen, sep = "")]]
  
  letter.i <- letter.i + 1
  
}

arranged.plots <- arrangeGrob(grobs = plots, 
                              layout_matrix = matrix(1:4,
                                                     nrow = 2,
                                                     ncol = 2,
                                                     byrow = TRUE))

path.to.save <- paste(base.path, "/figures", sep = "")
push.dir(path.to.save) # aux code

pdf.suffix <- "Figure_S1.pdf"

pdf(paste(path.to.save, "/",
          pdf.suffix, sep = ""),
    width = 7.3, height = 7)
grid.arrange(arranged.plots)
dev.off()
