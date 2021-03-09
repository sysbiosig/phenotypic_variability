# A script for reproducing the Figure ___ of a paper:

# "Phenotypic variability and not noise accounts for most of the cell-to-cell
# heterogeneity in cytokine signaling", Topolewski et al. 2021

#### loading libraries and functions ####
source("D:/Piotrek/scripts/basic_scripts/normalize_data.R")
source("D:/Piotrek/scripts/identity_analysis/pie_density.R")
source("D:/Piotrek/scripts/identity_analysis/noise_decompose_general.R")

#### set-up of directories ####
base.path <- "D:/Piotrek/Experiments/ICF/"
folder.name <- "cellular_identity"
project <- "fusion"
normaliz.chosen <- "pbs"
path <- paste(base.path, folder.name, sep='')

#### data loading ####
bio.AB <- read.csv(paste(path, "/", project, 
                         "/input/", normaliz.chosen,
                         "/fused_rename_paper.csv", sep = ""), header = TRUE, sep = ",")


bio.traj <- read.csv(paste(path, "/trajectories", 
                           "/input/", normaliz.chosen,
                           "/trajectories_paper.csv", sep = ""), header = TRUE, sep = ",")
bio.traj <- bio.traj %>% 
  dplyr::filter(nuclearity == "single-nuclear")
merge.chosen <- "no"

#### preparing the Figure S3 ####
ratio.limits <- c(0, 4)
differ.factor <- 0.5
bio.plot <- bio.AB
alpha <- 1
size <- 0.5
plots <- list()
plots[["mixing"]] <- ggplot()+
  geom_point(data = bio.plot,
             aes(y = (ratio_dye_B),
                 x = (ratio_dye_A),
                 color = factor(merged_nuclei)),
             alpha = alpha,
             size = size)+
  theme_trajectories(axis.num.size = 8,
                     aspect.ratio = 1)+
  theme(axis.title=element_text(size=8, face = "plain"))+
  scale_color_manual(values = c("lightgoldenrod2",
                                "mediumorchid4"),
                     name = "merged nuclei")+
  geom_abline(slope = 1, intercept = -differ.factor,
              linetype = 2, color = "grey40")+
  geom_abline(slope = 1, intercept = differ.factor,
              linetype = 2, color = "grey40")+
  scale_x_continuous(breaks = c(differ.factor, 1, seq(0, 4, 2)), limits = ratio.limits, expand = c(0, 0),
                     name = expression(paste(frac("nucleus A", "nucleus B"), " fluorescence of Dye A", sep = "")))+
  scale_y_continuous(breaks = c(differ.factor, 1, seq(0, 4, 2)), limits = ratio.limits, expand = c(0, 0),
                     name = expression(paste(frac("nucleus A", "nucleus B"), " fluorescence of Dye B", sep = "")))+
  guides(color = guide_legend(override.aes = list(size=8)))
plots[["mixing"]]


pdf.push(file = paste(path, "/", project, "/output/", normaliz.chosen, "/fused_ratio.pdf", sep = ""),
         width = 4, height = 3)
plots
dev.off()

#### Figure __ ####
plots <- list()
for(stim_type in c("IFNG", "OSM")){
  
  if(stim_type == "IFNG"){
    limits <- c(0, 60)
    x.text <- limits[1] + (limits[2] - limits[1])/30
    y.text <- limits[1] + (limits[2] - limits[1])*5/6
    breaks <- seq(limits[1], limits[2]-1, 20)
    my.pal <- IFN.pal
    time.chosen <- 15
    ggtitle <- bquote("IFN-"*gamma ~ .(time.chosen)~"min")
  } else if(stim_type == "OSM") { 
    limits <- c(0, 40)
    x.text <- limits[1] + (limits[2] - limits[1])/30
    y.text <- limits[1] + (limits[2] - limits[1])*5/6
    breaks <- seq(limits[1], limits[2]-1, 15)
    my.pal <- OSM.pal
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
              aes(label = paste("n = ", n, "     ", sep = "")),
              hjust = 0,
              size = size.text,
              family = "LM Roman 10")+
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
  # plots[[paste(stim_type, "_distinct", sep = "")]]
}

arranged.plots <- arrangeGrob(grobs = plots, nrow = 2, ncol = 1)
path.to.save <- paste(path, "/", project, "/output/", normaliz.chosen, sep = "")
push.dir(path.to.save)

cairo_pdf(paste(path.to.save, "/", project, 
                "_grid_distinct.pdf", sep = ""),
          width = 7.5, height = 5)
grid.arrange(arranged.plots)
dev.off()

#### Figure 2 ####
chosen.stimulants.AB <- c("IFNG", "OSM")
plots <- list()
for(stim.type.chosen in chosen.stimulants.AB){

  if(stim.type.chosen == "IFNG"){
    limits <- c(0, 60)
    x.text <- limits[1] + (limits[2] - limits[1])/30
    y.text <- limits[1] + (limits[2] - limits[1])*5/6
    breaks <- seq(limits[1], limits[2]-1, 20)
    my.pal <- IFN.pal
    time.chosen <- 15
    ggtitle <- bquote("IFN-"*gamma ~ .(time.chosen)~"min")
  } else if(stim.type.chosen == "OSM") {
    limits <- c(0, 40)
    x.text <- limits[1] + (limits[2] - limits[1])/30
    y.text <- limits[1] + (limits[2] - limits[1])*5/6
    breaks <- seq(limits[1], limits[2]-1, 15)
    my.pal <- OSM.pal
    time.chosen <- 30
    ggtitle <- bquote("OSM" ~ .(time.chosen)~"ng/mL")
  } else if(stim.type.chosen == "TNF") { 
    limits <- c(0, 7)
    x.text <- limits[1] + (limits[2] - limits[1])/30
    y.text <- limits[1] + (limits[2] - limits[1])*5/6
    breaks <- seq(limits[1], limits[2]-1, 3)
    my.pal <- TNF.pal
    time.chosen <- 30
    ggtitle <- bquote("TNF-"*alpha ~ .(time.chosen)~"ng/mL")
  }
  bio.plot <- bio.AB %>%
    dplyr::filter(merged_nuclei == merge.chosen &
                       stim_type == stim.type.chosen)

  decom.fold <- noise.decompose.boot(data = bio.plot,
                                     n = 1,
                                     boot.no = 1000,
                                     estimator = "sd",
                                colname.A = "A",
                                colname.B = "B",
                                group.columns = c("stim_level",
                                                  "time"))
  decom.fold
  alpha <- 1
  size <- 1.5
  shape <- 16


  plots[[paste(stim.type.chosen, "_scatter", sep = "")]] <- ggplot()+
    geom_point(data = bio.plot,
               aes(x = A,
                   y = B,
                   color = factor(stim_level)),
               alpha = alpha,
               size = size,
               shape = shape)+
    geom_text(data = decom.fold,
              parse = TRUE,
              x = x.text,
              y = y.text,
              aes(label = paste("atop(''* rho[italic(phen.)] * ''== ", Icv,
                                ",''* rho[italic(noise)] * ''*''== ", Ncv, ")",
                                sep = "")),
              hjust = 0,
              size = 4,
              family = "LM Roman 10")+
    geom_text(data = decom.fold,
              parse = TRUE,
              x = x.text+limits[2]*0.45,
              y = y.text,
              aes(label = paste("phantom(0)%+-%", sd_Ncv, sep = "")),
              hjust = 0,
              size = 3.25,
              family = "LM Roman 10")+
    geom_abline(slope = 1, linetype = 2, color = "grey40")+
    theme_trajectories(aspect.ratio = 1)+
    theme(axis.text=element_text(size=6))+
    facet_grid(.~stim_level)+
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
  # plots[[paste(stim.type.chosen, "_scatter", sep = "")]]
}


# arranged.plots <- arrangeGrob(grobs = plots, layout_matrix = rbind(c(1, 1, 1),
#                                                                    c(2, 2, 2),
#                                                                    c(3, 3, 3)))
# 
# width <- 7.5
# height <- 5
# filename_suffix <- "_AB.pdf"

# arranged.plots <- arrangeGrob(grobs = plots, layout_matrix = rbind(c(1, 1, 1),
#                                                                    c(2, 2, 2),
#                                                                    c(3, 3, 3)))
# 
# width <- 7.5
# height <- 7.5
# filename_suffix <- "_AB_TNF.pdf"

## for TNF mouse only ##
# arranged.plots <- arrangeGrob(grobs = plots, 
#                               nrow = 1, ncol = 1)
# width <- 7.25
# height <- 2.5
# filename_suffix <- "_AB_TNF_only.pdf"

path.to.save <- paste(path, "/", project, "/output/", normaliz.chosen, sep = "")
push.dir(path.to.save)

cairo_pdf(paste(path.to.save, "/", project,
                filename_suffix, sep = ""),
    width = width, height = height)
grid.arrange(arranged.plots)
dev.off()

#### Plotting the Figure 3 ####
sd.coefficients <- list()
i = 1
chosen.stimulants.AB <- c("IFNG", "OSM")
# chosen.stimulants.AB <- c("TNF")
plots <- list()
for(stim.type.chosen in chosen.stimulants.AB){
  
  if(stim.type.chosen == "IFNG"){
    x.limits <- c(0, 60)
    y.limits <- c(0, 60)*0.5
    x.text <- x.limits[1] + (x.limits[2] - x.limits[1])/10
    y.text <- y.limits[1] + (y.limits[2] - y.limits[1])*4/6
    x.breaks <- seq(x.limits[1], x.limits[2]-1, 20)
    y.breaks <- seq(y.limits[1], y.limits[2]-1, 20)
    my.pal <- IFN.pal
    time.chosen <- 15
    ggtitle <- bquote("IFN-"*gamma ~ .(time.chosen)~"min")
  } else if(stim.type.chosen == "OSM") { 
    x.limits <- c(0, 40)
    y.limits <- c(0, 40)*0.5
    x.text <- x.limits[1] + (x.limits[2] - x.limits[1])/10
    y.text <- y.limits[1] + (y.limits[2] - y.limits[1])*4/6
    
    x.breaks <- seq(x.limits[1], x.limits[2]-1, 15)
    y.breaks <- seq(y.limits[1], y.limits[2]-1, 15)
    my.pal <- OSM.pal
    time.chosen <- 30
    ggtitle <- bquote("OSM" ~ .(time.chosen)~"ng/mL")
  } else if(stim.type.chosen == "TNF") { 
    x.limits <- c(0, 7)
    y.limits <- x.limits*0.5
    x.text <- x.limits[1] + (x.limits[2] - x.limits[1])/30
    y.text <- y.limits[1] + (y.limits[2] - y.limits[1])*4/6
    
    x.breaks <- seq(x.limits[1], x.limits[2], 3)
    y.breaks <- seq(y.limits[1], y.limits[2], 3)
    my.pal <- TNF.pal
    time.chosen <- 30
    ggtitle <- bquote("TNF-"*alpha ~ .(time.chosen)~"ng/mL")
  } else if(stim.type.chosen == "TNF_h") { 
    stim.type.chosen <- "TNF"
    x.limits <- c(0, 7)
    y.limits <- x.limits*0.5
    x.text <- x.limits[1] + (x.limits[2] - x.limits[1])/30
    y.text <- y.limits[1] + (y.limits[2] - y.limits[1])*4/6
    
    x.breaks <- seq(x.limits[1], x.limits[2], 3)
    y.breaks <- seq(y.limits[1], y.limits[2], 3)
    my.pal <- TNF.pal
    time.chosen <- 30
    ggtitle <- bquote("TNF-"*alpha ~ .(time.chosen)~"ng/mL")
  }
  
  bio.sd <- bio.AB %>%
    dplyr::filter(merged_nuclei == merge.chosen,
           stim_type == stim.type.chosen) %>%
    mutate(sd = sd.AB(A, B),
           mean = (A + B)/2)
  
  bio.lm <- bio.sd %>%
    summarise(a = summary(lm(sd ~ 0 + mean))$coefficients[[1]],
              sd_err = summary(lm(sd ~ 0 + mean))$coefficients[[2]])
  sd.coefficients[[stim.type.chosen]] <- bio.lm$a
  bio.plot <- bio.sd
  
  alpha <- 1
  size <- 1.5
  shape <- 16
  if(i == 4) {stim.type.chosen = "TNF_h" }
  plots[[paste(stim.type.chosen, "_sd", sep = "")]] <- ggplot()+
    geom_point(data = bio.plot,
               aes(x = mean,
                   y = sd,
                   color = factor(stim_level)),
               alpha = alpha,
               size = size,
               shape = shape)+
    geom_text(data = bio.lm,
              parse = TRUE,
              x = x.text,
              y = y.text,
              aes(label = paste("atop(''* sigma[italic(noise)]^(i) * ''%~~% ", 
                                round(a, 3),  " ~~ italic(bar(y)[i]))", sep = "")),
              hjust = 0,
              size = 4,
              family="LM Roman 10")+
    geom_text(data = bio.lm,
              parse = TRUE,
              x = x.text + x.limits[2]*0.25,
              y = y.text - y.limits[2]*0.04,
              aes(label = paste("phantom(0)%+-%", round(sd_err, 4), sep = "")),
              hjust = 0,
              size = 2.5,
              family="LM Roman 10")+
    geom_abline(data = bio.lm, aes(slope = a, intercept = 0), 
                linetype = 1, 
                color = "grey70",
                size = 0.75)+
    theme_trajectories(aspect.ratio = 1)+
    theme(axis.title.y = element_text(family = "LM Roman 10", 
                                      margin = margin(0,-5,0,0), size = 12),
          axis.title.x = element_text(family = "LM Roman 10"))+
    scale_color_manual(values = my.pal,
                       guide = "none")+
    scale_x_continuous(expand = c(0, 0),
                       breaks = x.breaks, 
                       limits = x.limits,
                       name = expression(italic(bar(y)[i])))+
    scale_y_continuous(expand = c(0, 0),
                       breaks = y.breaks, 
                       limits = y.limits,
                       name = expression(sigma[italic(noise)]^(i)))
  # plots[[paste(stim.type.chosen, "_sd", sep = "")]]
  i = i + 1
}

#### percentile responses ####
i = 1
chosen.stimulants.trajectory <- c("IFNG", "OSM")
# chosen.stimulants.trajectory <- c("TNF")
for(stim.type.chosen in chosen.stimulants.trajectory){
  if(stim.type.chosen == "IFNG"){
    y.limits <- c(0, 17)
    my.pal <- IFN.pal
    time.chosen <- 15
    
    x.text <- 0.4
    y.text <- 0.85
    
    y.breaks <- seq(x.limits[1], x.limits[2], 6)
    x.breaks <- seq(y.limits[1], y.limits[2], 0.5)
    
    ggtitle <- bquote("IFN-"*gamma ~ .(time.chosen)~"min")
  } else if(stim.type.chosen == "OSM") { 
    y.limits <- c(0, 25)
    my.pal <- OSM.pal
    time.chosen <- 30
    x.text <- 0.4
    y.text <- 0.85
    
    y.breaks <- seq(x.limits[1], x.limits[2], 10)
    x.breaks <- seq(y.limits[1], y.limits[2], 0.6)
    ggtitle <- bquote("OSM" ~ .(time.chosen)~"min")
  } else if(stim.type.chosen == "TNF") { 
    y.limits <- c(0, 7)
    my.pal <- TNF.pal
    time.chosen <- 30
    x.text <- 0.4
    y.text <- 0.85
    
    y.breaks <- seq(x.limits[1], x.limits[2], 3)
    x.breaks <- seq(y.limits[1], y.limits[2], 0.6)
    ggtitle <- bquote("TNF"*alpha ~ .(time.chosen)~"min")
  } else if(stim.type.chosen == "TNF_h") { 
    y.limits <- c(0, 7)
    my.pal <- TNF.pal
    time.chosen <- 30
    x.text <- 0.4
    y.text <- 0.85
    
    y.breaks <- seq(x.limits[1], x.limits[2], 3)
    x.breaks <- seq(y.limits[1], y.limits[2], 0.6)
    ggtitle <- bquote("TNF"*alpha ~ .(time.chosen)~"min")
  }
  
  bio.percentile <- bio.traj %>% 
    dplyr::filter(time == time.chosen,
           stim_type == stim.type.chosen) %>%
    group_by(stim_level)%>%
    summarise(th5 = quantile(x = response, probs = 0.05),
              th50 = quantile(x = response, probs = 0.5),
              th95 = quantile(x = response, probs = 0.95)) %>%
    melt(id.vars = c("stim_level"), 
         value.name = "response",
         variable.name = "percentile")
  
  legend.size <- 8
  # if(i == 4) {stim.type.chosen = "TNF_h" }
  plots[[paste(stim.type.chosen, "_percentile", sep = "")]] <- 
    ggplot(data = bio.percentile,
           aes(y = response,
               x = factor(stim_level)))+
    geom_line(aes(group = percentile),
              color = "black",
              lwd = 0.5)+
    geom_point(aes(shape = percentile,
                   fill = factor(stim_level)),
               size = 3)+
    scale_fill_manual(values = my.pal,
                      name = "dose [ng/mL]",
                      guide = FALSE)+
    scale_shape_manual(values = c(24, 23, 21),
                       limits = c("th95", "th50", "th5"),
                       name = "percentile",
                       labels = c("95% (strong response)",
                                  "50% (moderate response)", 
                                  "5% (weak response)"))+
    theme_trajectories(aspect.ratio = 1,
                       border.thickness = 0.4)+
    theme(legend.position = c(x.text, y.text),
          legend.text = element_text(size = legend.size),
          legend.title = element_text(size = legend.size+2),
          legend.background = element_blank(),
          legend.key.size = unit(3, "mm"))+
    coord_cartesian(ylim = y.limits)+
    scale_x_discrete(expand = c(0.05, 0),
                     name = "dose [ng/mL]")+
    scale_y_continuous(expand = c(0.05, 0),
                       breaks = y.breaks,
                       name = "response [a.u.]")+
    guides(shape = guide_legend(override.aes = list(size=2)))
  # plots[[paste(stim.type.chosen, "_percentile", sep = "")]]
  i = i + 1
}

i=1
for(stim.type.chosen in chosen.stimulants.trajectory){
  
  
  if(stim.type.chosen == "IFNG"){
    my.pal <- IFN.pal
    time.chosen <- 15
    paste("dose", stim.type.chosen, " [ng/mL]")
    xtitle <- bquote("IFN-"*gamma ~ "dose [ng/mL]")
  } else if(stim.type.chosen == "OSM") { 
    my.pal <- OSM.pal
    time.chosen <- 30
    xtitle <- bquote("OSM dose [ng/mL]")
  } else if(stim.type.chosen == "TNF") { 
    my.pal <- TNF.pal
    time.chosen <- 30
    xtitle <- bquote("TNF-"*alpha ~ "dose [ng/mL]")
  } else if(stim.type.chosen == "TNF_h") { 
    stim.type.chosen <- "TNF"
    my.pal <- TNF.pal
    time.chosen <- 30
    xtitle <- bquote("TNF-"*alpha ~ "dose [ng/mL]")
  }
  sd.coefficient <- sd.coefficients[[stim.type.chosen]]
  if(i == 4) {stim.type.chosen = "TNF_h" }
  
  CM.capacity <- confusion.matrix.gamma(data = bio.traj,
                                   stim_type_chosen = stim.type.chosen,
                                   time.chosen = time.chosen,
                                   sd.coefficient = sd.coefficient,
                                   percentiles = seq(0.05, 0.95, by = 0.01),
                                   n.sampled = 100)
  bio.CM <- CM.capacity$CM
  gamma.capacity <- CM.capacity$capacity
  gamma.capacity.continuous <- CM.capacity$capacity.continuous
  
  if(0 %in% bio.CM$prob){
    bio.CM[bio.CM$prob == 0, ]$prob <- NA 
  }
   
  plots[[paste(stim.type.chosen, "CM")]] <- ggplot(bio.CM,
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
                       name = "dose for which response is typical [ng/mL]")+
    scale_x_continuous(breaks = NULL,
                       name = xtitle)+
    geom_vline(xintercept = 1.45,
               colour = "black") +
    scale_fill_manual(values = my.pal,
                      guide = "none")+
    labs(caption = bquote(atop(C^"*"* " = " * .(round(gamma.capacity, 2)),
                               C[cont.]^"*"* " = " * .(round(gamma.capacity.continuous, 2)))))+
    theme_trajectories(axis.name.size = 8)+
    theme(panel.border = element_blank(),
          panel.spacing.x=unit(0, "mm"),
          panel.spacing.y=unit(0, "mm"),
          plot.caption = element_text(size = 20,
                                      hjust = 0.5))
  plots[[paste(stim.type.chosen, "CM")]]
  i = i + 1
}

arranged.plots <- arrangeGrob(grobs = plots, layout_matrix = rbind(c(1, 2),
                                                                   c(3, 4),
                                                                   c(5, 6)))

width <- 5
height <- 9
filename.suffix <- "_sd_capacity2.pdf"

## for OSM, IFNG, TNFm and TNFh ##
# arranged.plots <- arrangeGrob(grobs = plots, layout_matrix = rbind(c(1, 2, 3, 4), 
#                                                                    c(5, 6, 7, 8), 
#                                                                    c(9, 10, 11, 12)))
# width <- 7.5*1.5
# height <- 8
# filename.suffix <- "_sd_TNF.pdf"

# ## for TNFm only ##
# arranged.plots <- arrangeGrob(grobs = plots, layout_matrix = rbind(c(1, 2, 3)))
# width <- 7.3
# height <- 2.5
# filename.suffix <- "_sd_TNFm_only.pdf"

path.to.save <- paste(path, "/", project, "/output/", normaliz.chosen, sep = "")
push.dir(path.to.save)

cairo_pdf(paste(path.to.save, "/", project, 
                filename.suffix, sep = ""),
          width = width, height = height)
grid.arrange(arranged.plots)
dev.off()

# jpeg(paste(path.to.save, "/", project, "_", merge.chosen,
#            "_sd_TNF2.jpeg", sep = ""),
#      width = 12*300, height = 8*300,
#      quality = 100,
#      res = 300,
#      type = "cairo")
# grid.arrange(arranged.plots)
# dev.off()

# 5; 15