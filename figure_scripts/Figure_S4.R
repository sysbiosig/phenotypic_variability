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

#### plotting the Figure S4 ####
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
                     aspect.ratio = 1)+ # aux. code
  theme(axis.title=element_text(size=8, face = "plain"))+
  scale_color_manual(values = c("lightgoldenrod2",
                                "mediumorchid4"),
                     name = "merged nuclei")+
  geom_abline(slope = 1, intercept = -differ.factor,
              linetype = 2, color = "grey40")+
  geom_abline(slope = 1, intercept = differ.factor,
              linetype = 2, color = "grey40")+
  scale_x_continuous(breaks = c(differ.factor, 1, seq(0, 4, 2)), 
                     limits = ratio.limits, expand = c(0, 0),
                     name = bquote(frac("nucleus#1", "nucleus#2") ~ 
                                     "inter-nuclear ratio of Dye A"))+
  scale_y_continuous(breaks = c(differ.factor, 1, seq(0, 4, 2)), 
                     limits = ratio.limits, expand = c(0, 0),
                     name = bquote(frac("nucleus#1", "nucleus#2") ~ 
                                     "inter-nuclear ratio of Dye B"))+
  guides(color = guide_legend(override.aes = list(size=8)))
plots[["mixing"]]

path.to.save <- paste(base.path, "/figures/", sep = "")
push.dir(path.to.save)
pdf(file = paste(path.to.save, "/Figure_S4.pdf", sep = ""),
    width = 4, height = 3)
plots
dev.off()