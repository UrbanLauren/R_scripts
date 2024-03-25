# RRBS plots
# Install and load libraries 
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
install.packages("readxl")
install.packages("writexl")
install_github("dgrtwo/broom")
# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
# usage
packages <- c("ggsci","devtools", "ggplot2", "ggpubr", "dplyr", "plyr",
              "tidyverse", "readxl","readr", "tibble", "readr", "reshape2",
              "RColorBrewer", "scales","sm", "writexl")
ipak(packages)

# Make a correlation scatterplot
# create a list of genes you are interested in
top20labels = (list of genes)

# Label just the top genes in results table
df_scatter_d1_dep <- mutate(df_scatter_d1_dep, labels = ifelse(GeneList %in% top20labels, GeneList, ""))

# make a correlation scatterplot
library(ggpubr)
sp <- ggscatter(df_scatter_d1_dep, x="methdiff_day1", y = "value",
                add = "reg.line",               # Add regression line
                conf.int = TRUE,  # Add confidence interval
                color = "DEP")+  # Color by groups "cyl"
  # Change point shape by groups "cyl"
  stat_cor(aes(color = DEP), label.x = 3)       # Add correlation coefficient
sp +   facet_wrap(~loc)+
  geom_text_repel(aes(label = labels)) + geom_vline(xintercept = c(-20, 20), linetype="dotted", 
                                                    color = "black", size=0.5)+ 
  geom_hline(yintercept = c(-1, 1), linetype="dotted", color = "black", size=0.5)+ 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


## data.matrix
## columns are individual samples (i.e. cells)
## rows are measurements taken for all the samples (i.e. genes)
data.matrix <- matrix(nrow=100, ncol=10)
colnames(data.matrix) <- c(
  paste("wt", 1:5, sep=""),
  paste("ko", 1:5, sep=""))
rownames(data.matrix) <- paste("gene", 1:100, sep="")
for (i in 1:100) {
  wt.values <- rpois(5, lambda=sample(x=10:1000, size=1))
  ko.values <- rpois(5, lambda=sample(x=10:1000, size=1))
  
  data.matrix[i,] <- c(wt.values, ko.values)
}
head(data.matrix)
dim(data.matrix)

pca <- prcomp(t(data.matrix), scale=TRUE) 

## plot pc1 and pc2
plot(pca$x[,1], pca$x[,2])

## make a scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

## now make a fancy looking plot that shows the PCs and the variation:
library(ggplot2)

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("My PCA Graph")

## get the name of the top 10 measurements (genes) that contribute
## most to pc1.
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes ## show the names of the top 10 genes

pca$rotation[top_10_genes,1] ## show the scores (and +/- sign)
