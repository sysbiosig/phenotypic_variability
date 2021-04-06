# A script for reproducing the Figure 3 of a paper:

# "Phenotypic variability and not noise accounts for most of the cell-to-cell
# heterogeneity in cytokine signaling", Topolewski et al. 2021

#### set-up of directories ####
base.path <- "D:/Piotrek/publications/syncytia_noise/phenotypic_variability"
input.path <- paste(base.path, "/data", sep = "")

#### loading libraries and functions ####
source(paste(base.path, "/Topolewski_auxillary_functions.R", sep = ""))
# all code noted as "aux code" is first introduced in the above auxillary file

#### data loading ####
bio.traj <- read.csv(paste(input.path,
                           "/trajectories_paper.csv", sep = ""),
                     header = TRUE, sep = ",")

bio.traj <- bio.traj %>% 
  dplyr::filter(nuclearity == "single-nuclear")

merge.chosen <- "no"

bio.AB <- read.csv(paste(input.path,
                         "/fused_paper.csv", sep = ""), 
                   header = TRUE, sep = ",")

#### Plotting the Figure 3 ####
sd.coefficients <- list()

chosen.stimulants.AB <- c("IFNG", "OSM")
plots <- list()
font.chosen <- "sans" #  
set.seed(12)

for(stim.type.chosen in chosen.stimulants.AB){
  
  if(stim.type.chosen == "IFNG"){
    x.limits <- c(0, 60)
    y.limits <- c(0, 60)*0.5
    x.text <- x.limits[1] + (x.limits[2] - x.limits[1])/10
    y.text <- y.limits[1] + (y.limits[2] - y.limits[1])*4/6
    x.breaks <- seq(x.limits[1], x.limits[2]-1, 20)
    y.breaks <- seq(y.limits[1], y.limits[2]-1, 20)
    my.pal <- IFN.pal # aux code
    time.chosen <- 15
    ggtitle <- bquote("IFN-"*gamma ~ .(time.chosen)~"min")
  } else if(stim.type.chosen == "OSM") { 
    x.limits <- c(0, 40)
    y.limits <- c(0, 40)*0.5
    x.text <- x.limits[1] + (x.limits[2] - x.limits[1])/10
    y.text <- y.limits[1] + (y.limits[2] - y.limits[1])*4/6
    
    x.breaks <- seq(x.limits[1], x.limits[2]-1, 15)
    y.breaks <- seq(y.limits[1], y.limits[2]-1, 15)
    my.pal <- OSM.pal # aux code
    time.chosen <- 30
    ggtitle <- bquote("OSM" ~ .(time.chosen)~"ng/mL")
  } 
  
  bio.sd <- bio.AB %>%
    dplyr::filter(merged_nuclei == merge.chosen,
           stim_type == stim.type.chosen) %>%
    mutate(sd = sd.AB(A, B), # aux code
           mean = (A + B)/2)
  
  bio.lm <- bio.sd %>%
    summarise(a = summary(lm(sd ~ 0 + mean))$coefficients[[1]],
              sd_err = summary(lm(sd ~ 0 + mean))$coefficients[[2]])
  sd.coefficients[[stim.type.chosen]] <- bio.lm$a
  bio.plot <- bio.sd
  
  alpha <- 1
  size <- 1.5
  shape <- 16

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
                                round(a, 3),  " ~~ italic(bar(y)[i]))", 
                                sep = "")),
              hjust = 0,
              size = 4,
              family=font.chosen)+
    geom_text(data = bio.lm,
              parse = TRUE,
              x = x.text + x.limits[2]*0.25,
              y = y.text - y.limits[2]*0.04,
              aes(label = paste("phantom(0)%+-%", round(sd_err, 4), sep = "")),
              hjust = 0,
              size = 2.5,
              family=font.chosen)+
    geom_abline(data = bio.lm, aes(slope = a, intercept = 0), 
                linetype = 1, 
                color = "grey70",
                size = 0.75)+
    theme_trajectories(aspect.ratio = 1)+
    theme(axis.title.y = element_text(family = font.chosen, 
                                      margin = margin(0,-5,0,0), size = 12),
          axis.title.x = element_text(family = font.chosen))+
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

}

#### percentile responses ####
i = 1
chosen.stimulants.trajectory <- c("IFNG", "OSM")
for(stim.type.chosen in chosen.stimulants.trajectory){
  if(stim.type.chosen == "IFNG"){
    y.limits <- c(0, 17)
    my.pal <- IFN.pal # aux code
    time.chosen <- 15
    
    x.text <- 0.4
    y.text <- 0.85
    
    y.breaks <- seq(x.limits[1], x.limits[2], 6)
    x.breaks <- seq(y.limits[1], y.limits[2], 0.5)
    
    ggtitle <- bquote("IFN-"*gamma ~ .(time.chosen)~"min")
  } else if(stim.type.chosen == "OSM") { 
    y.limits <- c(0, 25)
    my.pal <- OSM.pal # aux code
    time.chosen <- 30
    x.text <- 0.4
    y.text <- 0.85
    
    y.breaks <- seq(x.limits[1], x.limits[2], 10)
    x.breaks <- seq(y.limits[1], y.limits[2], 0.6)
    ggtitle <- bquote("OSM" ~ .(time.chosen)~"min")
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
                       border.thickness = 0.4)+ # aux code
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
                       name = "relative fluorescence intensity")+
    guides(shape = guide_legend(override.aes = list(size=2)))
  i = i + 1
}

i=1
for(stim.type.chosen in chosen.stimulants.trajectory){
  
  
  if(stim.type.chosen == "IFNG"){
    my.pal <- IFN.pal # aux code
    time.chosen <- 15
    paste("dose", stim.type.chosen, " [ng/mL]")
    xtitle <- bquote("IFN-"*gamma ~ "dose [ng/mL]")
  } else if(stim.type.chosen == "OSM") { 
    my.pal <- OSM.pal # aux code
    time.chosen <- 30
    xtitle <- bquote("OSM dose [ng/mL]")
  } 
  sd.coefficient <- sd.coefficients[[stim.type.chosen]]
  
  bio.CM <- confusion.matrix.gamma(data = bio.traj,
                                   stim_type_chosen = stim.type.chosen,
                                   time.chosen = time.chosen,
                                   sd.coefficient = sd.coefficient,
                                   percentiles = seq(0.05, 0.95, by = 0.01)) 
  # aux code
  
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
    theme_trajectories(axis.name.size = 8)+
    theme(panel.border = element_blank(),
          panel.spacing.x=unit(0, "mm"),
          panel.spacing.y=unit(0, "mm"),
          plot.caption = element_text(size = 20,
                                      hjust = 0.5))
  i = i + 1
}

arranged.plots <- arrangeGrob(grobs = plots, layout_matrix = rbind(c(1, 2),
                                                                   c(3, 4),
                                                                   c(5, 6)))

width <- 5.25
height <- 8
filename.suffix <- "Figure_3.pdf"

path.to.save <- paste(base.path, "/figures", sep = "")
push.dir(path.to.save)

cairo_pdf(paste(path.to.save, "/",
                filename.suffix, sep = ""),
          width = width, height = height)
grid.arrange(arranged.plots)
dev.off()
