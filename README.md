# phenotypic_variability
The data sets, R scripts for channel capacity calculations and R scripts for figure reproducing of a paper:

"Phenotypic variability, not noise, accounts for most of the cell-to-cell heterogeneity in IFN-γ and oncostatin M signaling responses" 

Topolewski et al. 2022, Science Signaling

All required libraries are either at the beginning of scripts or in "Topolewski_auxiliary_functions.R" script. 
Loading of "Topolewski_auxiliary_functions.R" script is at the beginning of the plotting script if needed.

Each script is designed to be independent, that is no need for pre-running of any other scripts manually, 
all needed packages/functions are loaded together with the script execution. 

## Requirements - Software
The main software requirement is the installation of the R environment (version: >= 3.2), which can be downloaded from [R project website](https://www.r-project.org) and is distributed for all common operating systems.
The use of a dedicated Integrated development environment (IDE), e.g. [RStudio](https://www.rstudio.com) is recommended. 

Apart from a base installation of R, the sucessful runnig of all scripts requires the following R packages:
+ "SLEMI", which can be installed according to https://github.com/sysbiosig/SLEMI

CRAN packages should be installed automatically during running "Topolewski_auxiliary_functions.R" script.
If problem occurs, the following CRAN packages can be installed via: 
> "install.packages("name_of_a_package")":

+ "ggplot2", 
+ "gridExtra",
+ "data.table",
+ "reshape2",
+ "dplyr",
+ "scales",
+ "RColorBrewer",
+ "rlist",

## Description of folders:
- data: containing the data required to reconstruct the figures
- figure_scripts: contatining all R scripts for preparing the paper quantitative figures 
- figures: contatining all paper quantitative figures plotted in R. 
- capacity_scripts: containing scripts for channel capacity calculations

## Base directory setting up
It is recommended to run all scripts in "phenotypic_variability" R project, placed in the main directory, which assures the correct path structure. Otherwise, `base.path` has to be provided in each script

## Licence
The repository is released under the GNU 3.0 licence and is freely available.
