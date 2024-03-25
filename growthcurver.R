#LU - growthcurver 
#5/28/2020

#Copy excel data for Mac OSX system
#Step 1) Select and copy the data (Cmd + c)
#Step 2) Use the function pipe(pbpaste) to import the data you’ve copied (with Cmd + c):
my_data <- read.table(pipe("pbpaste"), sep="\t", header = TRUE, stringsAsFactors = FALSE)

#time needs to be in hours (0-20hrs), if you need to convert minutes -> hours see link below: 
#https://cran.r-project.org/web/packages/growthcurver/vignettes/Growthcurver-vignette.html#background-correction
#install program 
install.packages("growthcurver")
library(growthcurver)

#rename column named Time to time 
#d <-my_data %>% rename(time = Time)

#Normalize to the lowest value for each observation 
gc_out <- SummarizeGrowthByPlate(my_data, plot_fit = TRUE, plot_file = "synergy.pdf")
#look at the first few rows
head(gc_out)
#set the working directory
setwd("~/")
#save the entire data table to a tab-separated file that can be imported into Excel.
output_file_name <- "2020-09-09_RN4220_USA300_TSB_H2A.txt"
write.table(gc_out, file = "2020-10-30_USA300_TSB_Ramopl.txt", quote = FALSE, sep = "\t", row.names = FALSE)




