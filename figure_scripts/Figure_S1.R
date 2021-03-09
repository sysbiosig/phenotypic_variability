# A script for reproducing the Figure ___ of a paper:

# "Phenotypic variability and not noise accounts for most of the cell-to-cell
# heterogeneity in cytokine signaling", Topolewski et al. 2021

#### loading libraries and functions ####
source("D:/Piotrek/scripts/basic_scripts/normalize_data.R")
source("D:/Piotrek/scripts/identity_analysis/pie_density.R")
source("D:/Piotrek/scripts/identity_analysis/noise.decompose.R")

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
                           "_paper.csv", sep = ""))
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
    y.limits <- c(0, 17)
    my.pal <- IFN.pal
    y.breaks <- seq(y.limits[1], y.limits[2], 6)
    ggtitle <- bquote("IFN-"*gamma ~ .(stim.level.chosen)~"ng/mL, single cells")
  } else if(stim.type.chosen == "OSM") {
    y.limits <- c(0, 25)
    my.pal <- OSM.pal
    y.breaks <- seq(y.limits[1], y.limits[2], 10)
    ggtitle <- bquote("OSM" ~ .(stim.level.chosen)~"ng/mL, single cells")
  }

  bio.plot <- bio.traj %>%
    dplyr:filter(stim_level == stim.level.chosen &
                            stim_type == stim.type.chosen)
  x.breaks <- unique(bio.plot$time)

  width <- 9 * (y.limits[2] / 25)

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
    theme_trajectories(aspect.ratio = aspect.chosen)+
    scale_y_continuous(expand = c(0, 0),
                       breaks = y.breaks)+
    scale_x_continuous(breaks = x.breaks)+
    ylab("relative fluorescence intensity")+
    xlab("time [min]")+
    add_letter(LETTERS[letter.i])+
    coord_cartesian(ylim = y.limits,
                    xlim = c(-2, 92))+
    ggtitle(ggtitle)
  plots[[paste("single.nuc", stim.type.chosen, sep = "")]]

  letter.i <- letter.i + 2
}

#### plotting the time trajectory for bi-nuclear cells ####
letter.i <- 2
for(stim.type.chosen in c("IFNG", "OSM")){
  
  if(stim.type.chosen == "IFNG"){
    y.limits <- c(0, 40)
    my.pal <- IFN.pal
    y.breaks <- seq(y.limits[1], y.limits[2], 15)
    ggtitle <- bquote("IFN-"*gamma ~ .(stim.level.chosen)~"ng/mL, syncytia")
    width <- 9 * (y.limits[2] / 25)
  } else if(stim.type.chosen == "OSM") {
    y.limits <- c(0, 25)
    my.pal <- OSM.pal
    y.breaks <- seq(y.limits[1], y.limits[2], 10)
    ggtitle <- bquote("OSM" ~ .(stim.level.chosen)~"ng/mL, syncytia")
    width <- 9 * (y.limits[2] / 50)
  }
  
  bio.plot <- bio.traj %>%
    dplyr::filter(stim_level == stim.level.chosen &
                            stim_type == stim.type.chosen)
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
    labs(tag = LETTERS[letter.i])+
    theme_trajectories(aspect.ratio = aspect.chosen)+
    theme(plot.tag = element_text(size = 20, face = "bold"))+
    scale_y_continuous(expand = c(0, 0),
                       breaks = y.breaks)+
    scale_x_continuous(breaks = x.breaks)+
    ylab("relative fluorescence intensity")+
    xlab("time [min]")+
    coord_cartesian(ylim = y.limits,
                    xlim = c(-2, 92))+
    ggtitle(ggtitle)
  plots[[paste("bi.nuc", stim.type.chosen, sep = "")]]
  
  letter.i <- letter.i + 2
  
}

arranged.plots <- arrangeGrob(grobs = plots, layout_matrix = matrix(1:4,
                                                                    nrow = 2,
                                                                    ncol = 2,
                                                                    byrow = FALSE))
path.to.save <- paste(path, "/", project, "/output/", normaliz.chosen, sep = "")
push.dir(path.to.save)

pdf.suffix <- "FIgure_S1.pdf"

pdf(paste(path.to.save, "/",
          pdf.suffix, sep = ""),
    width = 7.3, height = 7)
grid.arrange(arranged.plots)
dev.off()
