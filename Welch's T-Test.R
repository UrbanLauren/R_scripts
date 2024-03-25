# Welch's t test - doesn't assume homogeneity of variances 

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
packages <- c("broom","devtools", "ggplot2", "ggpubr", "grid", "plyr",
              "tidyverse", "readxl","rstatix", "tibble", "readr", "reshape2",
              "RColorBrewer", "scales", "writexl")
ipak(packages)

# Import data
# Set path to location of excel file
Fig2b <- read_excel(path = xl_data, sheet = "Fig2b", col_names = TRUE, col_types = "numeric")

# Remove excel columns and rows 
Fig2b <-as.data.frame(Fig2b[1:4,3:5]) #Remove condition label & Avg/SEM rows

#Step1 - transform the data
NewFig2b <- t(Fig2b)
dimnames(NewFig2b)

# Add conditions as column name 
colnames(NewSfig2) <- c("PBGT", "Unstimulated", "PMA Stimulated", "PBGT+PMA")
colnames(NewFig2b) <- c("MinA", "NaOCl", "NaOCl + BME 1uM", "NaOCl + BME 10uM", "NaOCl + BME 100uM",
                        "NaOCl + Methionine 1uM","NaOCl + Methionine 10uM","NaOCl + Methionine 100uM",
                        "NaOCl + Cysteine 0.1uM", "NaOCl + Cysteine 1uM", "NaOCl +Cysteine 10uM")
# Convert to dataframe
NewFig2b<-as.data.frame(NewFig2b)

#Normality test - if the p-value > 0.05 then you can assume normality
shapiro.test(NewFig2B$`MinA+NAOCl + FeCl2`) #p-value = 0.7567
shapiro.test(NewFig2B$`MinA+NAOCl + FeCl3`) #p-value = 0.6842

# Use F Test to Compare Two Variances
# the more this ratio deviates from 1, the stronger the evidence for unequal variances.
var.test(NewFig2B$`MinA 1:1000`, NewFig2B$NaCl) #F = 0.042759
var.test(NewFig2B$`MinA 1:1000`, NewFig2B$HCl) #F = 0.5518
var.test(NewSfig2$`MinA 1:1000`, NewSfig2$HNO3) #F = 1.8843
var.test(NewSfig2$`MinA 1:1000`, NewSfig2$H2SO4) #F = 0.0020528

# t-test for New Supplemental Figs
t1<-t.test(NewSfig2$NaOCl, NewSfig2$`NaOCl + BME 1uM`) #p-value = 0.167
t2<-t.test(NewSfig2$NaOCl, NewSfig2$`NaOCl + Methionine 1uM`) #p-value = 0.167
t3<-t.test(NewSfig2$NaOCl, NewSfig2$`NaOCl + Cysteine 1uM`) #p-value = 0.167

# Export welch t-test output 
capture.output(t1, file = "200831_NewSupFig_WelchTest.csv", append = TRUE)
cat("\n\n", file = "200831_NewSupFig_WelchTest.csv", append = TRUE)
capture.output(t2, file = "200831_NewSupFig_WelchTest.csv", append = TRUE)
cat("\n\n", file = "200831_NewSupFig_WelchTest.csv", append = TRUE)
capture.output(t3, file = "200831_NewSupFig_WelchTest.csv", append = TRUE)
cat("\n\n", file = "200831_NewSupFig_WelchTest.csv", append = TRUE)