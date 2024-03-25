# GO Enrichment analysis for RRBS paper
# Install and load libraries 
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
install.packages("readxl")
install.packages("writexl")
install.packages("devtools")
devtools::install_github("slowkow/ggrepel")
library(kableExtra)
options(knitr.table.format = "html")
install.packages("gprofiler2")

# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
# usage
packages <- c("ggsci","devtools", "hrbrthemes", "ggplot2", "ggpubr", "dplyr", "ggrepel",
              "tidyverse", "readxl","readr", "tibble", "readr", "reshape2", "gprofiler2",
              "RColorBrewer", "scales","sm", "viridis", "writexl", "GGally")
ipak(packages)

# Import data
w1_hypo_genesmp <- data.frame(fdep_w1_hypo_genesmp)
w1_hyper_genesmp <- data.frame(fdep_w1_hyper_genesmp)
d_hypo_genesmp <- data.frame(fdep_d_hypo_genesmp)
d_hyper_genesmp <- data.frame(fdep_d_hyper_genesmp)
dw_hypo_genesmp <- data.frame(fdep_dw_hypo_genesmp)
dw_hyper_genesmp <- data.frame(fdep_dw_hyper_genesmp)

# Add jitter points 
par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(5, 1, 0))

d_hypo_genesmp<-d_hypo_genesmp %>% arrange(d_hypo_genesmp, 'methdiff_day1')
# Violin Plot
p<-ggviolin(d_hypo_genesmp, x='GeneList', y='methdiff_day1',
         add = "jitter",orientation = "horiz", color = 'GeneList', 
  font.label = list(size = 4, color = "black")) +
  theme(legend.position = 'none')    #remove the legend

p<-ggviolin(w1_hyper_genesmp, x='GeneList', y='methdiff_wk1',
            add = "jitter",orientation = "horiz", color = 'GeneList', 
            font.label = list(size = 4, color = "black")) +
  theme(legend.position = 'none')    #remove the legend

# Order data
tmp <- dw_hyper_genesmp %>%
  arrange(desc('methdiff_d1wk1')) %>%
  mutate('GeneList'=factor('GeneList', 'GeneList'))

# Geom_plot
ggplot(tmp, aes(x=as.factor('GeneList'), y='methdiff_d1wk1')) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", fill=alpha("#69b3a2", 0.8)) +
  ylim(-7000,13000) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar(start = 0) 
x = dw_hyper_genesmp$GeneList
y = dw_hyper_genesmp$methdiff_d1wk1

p <- ggplot(dw_hyper_genesmp, aes(x=x, y=y)) +
  geom_bar()+ylim(-100,120)

ggplot(dw_hypo_genesmp, aes(x='GeneList', y='methdiff_d1wk1')) +
         geom_segment( aes(x='GeneList', xend='GeneList', y=-20, yend='methdiff_d1wk1'), color = 'GeneList') +
         geom_point( color = 'GeneList', size=4, alpha=0.6) +
  coord_flip() +
  theme_light() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

# GO Enrichment analysis
gostres <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)

pp <- publish_gostplot(p, highlight_terms = c("GO:0048013", "REAC:R-HSA-3928663"), 
                       width = NA, height = NA, filename = NULL )

# Import MethylDiff data------------------- 
df_UniqueDMCs_Wk1_base_Hypo <- as.data.frame(UniqueDMCs_Wk1_base_Hypo)
df_genelist <- as.data.frame(genelist_Wk1_base_Hypo)
# remove duplicate rows 
df_UniqueDMCs_Wk1_base_Hypo <- unique(df_UniqueDMCs_Wk1_base_Hypo)
# Combine the dataframes
UniqueDMCs_Wk1_base_Hypo <- left_join(df_UniqueDMCs_Wk1_base_Hypo,
                                               df_genelist,
                                               by = 'Gene')
# drop na's for now
UniqueDMCs_Wk1_base_Hypo <- na.omit(UniqueDMCs_Wk1_base_Hypo)
UniqueDMCs_Wk1_base_Hypo <- unique(UniqueDMCs_Wk1_base_Hypo)
UniqueDMCs_Wk1_base_Hypo$Pathway <- as.factor(UniqueDMCs_Wk1_base_Hypo$Pathway)

# Subset of data
df_UniqueDMCs_Wk1_base_Hypo <-subset(UniqueDMCs_Wk1_base_Hypo,
                    Pathway=="EGF receptor signaling pathway",
                    select = c("nsc_methylation", "d1_methylation", "wk1_methylation", "Gene"))
# Make into tidy data
tidydf_merged_EGF <-merged_EGF %>%
  gather(group, methylation, 'nsc_methylation':'wk1_methylation', convert = TRUE)
tidydf_merged_EGF$group <- as.factor(tidydf_merged_EGF$group)

plot1 <- ggplot(data=tidydf_merged_EGF, aes(x=group, y=methylation, colour=Gene))+ 
  geom_jitter()+
  #geom_dotplot(binwidth = 1.5, stackratio = .7)+
  theme(axis.text.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")+
  ylab("Methylation level")+
  ggtitle("EGF receptor signaling pathway")

# Plot coordinates
ggparcoord(merged_EGF,
           columns = 1:3, groupColumn = 4,
           showPoints = TRUE, 
           scale="globalminmax",
           title = "Parallel Coordinate Plot",
           alphaLines = 0.3) + 
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum()+
  theme(
    plot.title = element_text(size=10))



