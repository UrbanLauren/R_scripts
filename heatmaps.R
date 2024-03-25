## RRBS manuscript
# ------ HOMER plots ------------

# Install and load libraries 
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
install.packages("readxl")
install.packages("writexl")
install.packages("tidyverse")
install.packages("viridis")
install.packages("paletteer")
install.packages("superheat")

# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
# usage
packages <- c("colorspace", "ggsci", "lattice", "superheat", "devtools", "heatmaply",
              "ggrepel", "gplots","ggplot2", "ggpubr", "dplyr", "plyr", "paletteer", 
              "tidyverse", "readxl","readr", "tibble", "readr", "reshape2",
              "RColorBrewer", "scales","sm", "viridis", "writexl")
ipak(packages)
#colors <- rev(brewer.pal(9, "Greys")) 

# Import data
homer_hypo <- data.frame(homer_hypo)
homer_hyper <- data.frame(homer_hyper)
# Replace missing values with "NA"
homer_hyper <- data.frame(homer_hyper,  na.strings = "NA") # Make sure to import strings as "NA" for the Unique DMCs (Hyper)**
View(homer_hypo)

# reorder by column index
homer_hypo <- homer_hypo[, c(3,2,1)] # leave the row index blank to keep all rows
# set rows names and remove that column 
rownames(homer_hypo) <- homer_hypo$terms
df_homer_hypo<-subset(homer_hypo, select = -c(terms))

# convert into matrix
mat_df_homer_hypo_Top30 <- as.matrix(df_homer_hypo)
View(mat_df_homer_hypo_Top30)

# create the heat map (Hypo DMCs)
superheat(mat_df_homer_hypo_Top30, left.label.text.size = 2, bottom.label.text.size = 4,
          heat.lim= c(40, 510),
          #heat.pal = viridis::inferno(50, direction =-1),
          heat.col.scheme = "grey",
          #heat.pal = "colors", 
          scale = FALSE,
          bottom.label.text.angle = 90, grid.hline.col = "white",
          grid.vline.col = "white", 
          heat.na.col = "black")

# Repeat for Unique DMCs (Hyper)
View(homer_hyper)
#reorder by column index
homer_hyper <- homer_hyper[, c(3,2,1,4)] # leave the row index blank to keep all rows

# set rows names and remove that column 
rownames(homer_hyper) <- homer_hyper$terms
df_homer_hyper<-subset(homer_hyper, select = -c(terms, na.strings))
# convert into matrix
mat_df_homer_hyper_Top30 <- as.matrix(df_homer_hyper)
View(mat_df_homer_hyper_Top30)

# create the heat map
superheat(mat_df_homer_hyper_Top30, left.label.text.size = 2, bottom.label.text.size = 4,
          heat.lim= c(5, 150),
          #heat.pal = viridis::inferno(50, direction =-1),
          heat.col.scheme = "grey",
          #heat.pal = "colors", 
          scale = FALSE,
          bottom.label.text.angle = 90, grid.hline.col = "white",
          grid.vline.col = "white", 
          heat.na.col = "black")


