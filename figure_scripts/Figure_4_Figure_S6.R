# A script for reproducing the Figure 4 and S6 of a paper:

# "Phenotypic variability and not noise accounts for most of the cell-to-cell
# heterogeneity in cytokine signaling", Topolewski et al. 2021

#### loading libraries and functions ####
source("D:/Piotrek/scripts/basic_scripts/normalize_data.R")

#### set-up of directories ####
base.path <- "D:/Piotrek/Experiments/cytometry/"
folder.name <- "molecular_phenotype"
short.folder.name <- "molphen"
path <- paste(base.path, folder.name, sep='')
#### data loading ####
bio.all <- read.csv(paste(path, "/input/", 
                          short.folder.name,
                          "_papeR.csv", sep = ""))

#### providing replicate details ####
## for each cell line there is separate limits to each axis ##
  y.BJ.lim <- c(-1, 2)
  x.BJ.lim <- c(-1.5, 1.5)
  
  y.MEF.lim <- c(-1, 3)
  x.MEF.lim <- x.BJ.lim
  receptor.range <- c(-1, 1)

rep.details <- list(replicate = c("FC-PT63", 
                                   "FC-PT74", 
                                   "FC-PT75",
                                   "FC-PT79"),
                    cell_type = c("MEF",
                                  "BJ",
                                  "MEF",
                                  "BJ"),
                    stim_type = c("IFN-gamma",
                                  "OSM",
                                  "IFN-gamma",
                                  "OSM"),
                    receptor_name = c("IFNGR",
                                      "OSMR",
                                      "IFNGR",
                                      "OSMR"),
                    STAT_name = c("STAT1",
                                  "STAT3",
                                  "STAT1",
                                  "STAT3"),
                    pSTAT_name = c("phospho-STAT1",
                                   "phospho-STAT3",
                                   "phospho-STAT1",
                                   "phospho-STAT3"),
                    x_lim = list(x.MEF.lim,
                                 x.BJ.lim,
                                 x.MEF.lim,
                                 x.BJ.lim),
                    y_lim = list(y.MEF.lim,
                                 y.BJ.lim,
                                 y.MEF.lim,
                                 y.BJ.lim))


#### choosing the stim_level levels - remove to paper ####
bio.relevant <- bio.all %>% 
  dplyr::filter(stim_level %in% c(0, 10))

#### calculating R2 values ####
bio.correlations <- bio.relevant %>%
  group_by(stim_level, replicate) %>% 
  summarize(multilinear_R2 = summary(lm(log10(phospho_STAT) ~
                                          log10(total_STAT) +
                                          log10(receptor)))$r.squared,
            phospho_vs_total = summary(lm(log10(phospho_STAT)~
                                   log10(total_STAT)))$r.squared,
            phospho_vs_receptor = summary(lm(log10(phospho_STAT)~
                                               log10(receptor)))$r.squared,
            total_vs_receptor = summary(lm(log10(total_STAT)~
                                             log10(receptor)))$r.squared) %>%
  left_join(as.data.frame(rep.details)[, 1:6], by = "replicate") # providing some details about replicates

#### building the convenient data frame with details of replicates ####:
bio.names <- bio.relevant %>%
  group_by(replicate, stim_level) %>%
  summarise() %>%
  left_join(as.data.frame(rep.details)[, 1:6], by = "replicate") %>%
  ungroup()

#### preparing for plotting plotting ####
## all plots will be saved in a list ##
# all.exp.list <- list() 
chosen.replicates <- c("FC-PT63",
                        "FC-PT79")
## size of a point and iterating variable for figure's letter
point.size <- 0.5
letter.i <- 1
plots <- list()
for(replicate.name in chosen.replicates){
  ## iterating variable for replicate name
  i <- which(rep.details$replicate == replicate.name)
  cell.type.chosen <- rep.details$cell_type[i]
  
  ## subsetting the data ## 
  bio.plot <- bio.relevant %>%
    dplyr::filter(replicate == replicate.name)
  bio.names.plot <- bio.names[bio.names$replicate == replicate.name, ]
  
  ## getting the xy limits from experiment details ##
  x.limits <- rep.details$x_lim[[i]]
  y.limits <- rep.details$y_lim[[i]]
  
  ## providing the position of text annotation on a plot ##
  x.text <- x.limits[1] + (x.limits[2] - x.limits[1])/12
  y.text <- y.limits[1] + (y.limits[2] - y.limits[1])*6/10
  
  ## the specific name of a receptor ##
  receptor.name <-  as.character(unique(bio.names.plot$receptor_name))
  
  ## exact plotting ##
  plots[[replicate.name]] <- ggplot()+
    geom_point(data = bio.plot,
               aes(y = phospho_STAT,
                   x = total_STAT,
                   color = receptor),
               size = point.size,
               alpha = 1)+
    geom_text(data = bio.correlations[bio.correlations$replicate == replicate.name, ],
              x = x.limits[1]+ (x.limits[2] - x.limits[1])/4*2,
              y = y.limits[1]+ (y.limits[2] - y.limits[1])/5*4.75,
              parse = TRUE,
              aes(label = paste("R[", pSTAT_name, "%~%", receptor_name, "+", STAT_name, "]^2 ==",
                                round(multilinear_R2, 2), sep = "")),
              size = 5,
              color = "black")+
    geom_text(data = bio.correlations[bio.correlations$replicate == replicate.name, ],
              x = x.limits[1]+ (x.limits[2] - x.limits[1])/4*0.9,
              y = y.limits[1]+ (y.limits[2] - y.limits[1])/5*3.8,
              parse = TRUE,
              aes(label = paste("atop(atop(textstyle(R[", pSTAT_name, "%~%", receptor_name, "]^2 ==",
                                round(phospho_vs_receptor, 2), "),textstyle(",  
                                "R[", pSTAT_name, "%~%", STAT_name, "]^2 ==",
                                round(phospho_vs_total, 2),
                                ")), ",
                                "textstyle(R[", STAT_name, "%~%", receptor_name, "]^2 ==",
                                round(total_vs_receptor, 2),
                                ")~phantom(0)~phantom(0)~phantom(0)~phantom(0))", sep = "")),
              size = 3,
              color = "black")+
    geom_text(data = bio.names.plot,
              x = x.limits[1]+ (x.limits[2] - x.limits[1])/4*2.5,
              y = y.limits[1]+ (y.limits[2] - y.limits[1])/20,
              parse = TRUE,
              aes(label = paste(cell_type, "+", stim_level, " ~ ng/mL ~ ", stim_type, sep = "")),
              size = 5,
              fontface = "bold",
              color = "black")+
    scale_x_log10(labels = trans_format("log10", math_format(10^.x)),
                  limits = 10^x.limits,
                  name = paste(unique(bio.names.plot$STAT_name), 
                               " relative fluorescence intensity", sep = "")) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                  limits = 10^y.limits,
                  name = paste(unique(bio.names.plot$pSTAT_name), 
                               "\nrelative fluorescence intensity", sep = "")) +
    labs(tag = LETTERS[letter.i])+
    theme_trajectories(aspect.ratio = 1)+
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.title = element_text(size = 9),
          legend.position = c(0.93, 0.4),
          legend.key.size = unit(4, "mm"),
          legend.background = element_blank(),
          plot.tag = element_text(size = 20, face = "bold"))+
    facet_grid(.~stim_level)+
    scale_color_distiller(name = paste(unique(bio.names.plot$receptor_name), 
                                       "\nrelative\nfluorescence\nintensity", sep = ""),
                          palette = "YlOrRd",
                          limits = 10^receptor.range,
                          direction = 1,
                          trans = "log10",
                          labels = trans_format("log10", math_format(10^.x)),
                          na.value = element_blank())
  letter.i <- letter.i + 1
}

arranged.plots <- marrangeGrob(grobs = plots, 
                               nrow = 2, ncol = 1,
                               top = element_blank())
path.to.save <- paste(path, "/output/", sep = "")
push.dir(path.to.save)
pdf(paste(path.to.save, "/", "Figure_4.pdf", sep = ""),
    width = 7.3, height = 8)
arranged.plots
dev.off()

#### distributions of protein in 0 vs 10 ng ####
plots <- list()
letter.i <- 1
for(replicate.name in chosen.replicates){
  
  bio.plot <- bio.relevant %>%
    dplyr::filter(replicate %in% replicate.name) %>%
    mutate(cell_and_stim_type = paste(cell_type, stim_type, sep = "+")) %>%
    rename(pSTAT = phospho_STAT,
           STAT = total_STAT,
           receptors = receptor) %>%
    melt(measure.vars = c("pSTAT",
                          "STAT",
                          "receptors"), 
         variable.name = "fluorophore", 
         value.name = "fluorescence")
  bio.plot$cell_and_stim_type <- factor(bio.plot$cell_and_stim_type,
                                        levels = c("MEF+IFNG",
                                                   "BJ+OSM"),
                                        labels = c("MEF+IFN-gamma",
                                                   "BJ+OSM"))
  
  if(unique(bio.plot$stim_type) == "IFNG"){
    levels(bio.plot$fluorophore) <- c("phospho-STAT1",
                                      "STAT1",
                                      "IFNGR")
  } else {
    levels(bio.plot$fluorophore) <- c("phospho-STAT3",
                                      "STAT3",
                                      "OSMR")
  }
  
  plots[[replicate.name]] <- ggplot()+
    geom_line(stat = "density",
              data = bio.plot,
              aes(x = fluorescence,
                  linetype = factor(stim_level)),
              adjust = 2)+
    facet_grid(cell_and_stim_type~fluorophore,
               labeller = label_parsed)+
    scale_x_log10(labels = trans_format("log10", math_format(10^.x)),
                  limits = 10^c(-1.5, 3),
                  name = "relative fluorescence intensity") +
    labs(tag = LETTERS[letter.i])+
    theme_trajectories(aspect.ratio = 1)+
    theme(legend.position = c(0.93, 0.5),
          legend.key.size = unit(4, "mm"),
          legend.background = element_blank(),
          strip.text.y = element_text(size = 12),
          plot.tag = element_text(size = 20, face = "bold"))+
    scale_linetype_discrete(name = "dose\n[ng/mL]")
  plots[[replicate.name]]
  
  letter.i <- letter.i + 1
}

arranged.plots <- marrangeGrob(grobs = plots, 
                               nrow = 2, ncol = 1,
                               top = element_blank())
path.to.save <- paste(path, "/output/", sep = "")
push.dir(path.to.save)
pdf(paste(path.to.save, "/", "FlowCytometry_distributions_0vs10.pdf", sep = ""),
    width = 7.3, height = 6)
arranged.plots
dev.off()
