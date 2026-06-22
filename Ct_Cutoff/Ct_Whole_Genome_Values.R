#######################################################
# Conclusion Chapter TD Thesis
# Ct Cutoff Value - Whole Genome (yes/no)
# TD 02.26.26
#######################################################

# load required libraries
library(tidyverse)
library(readxl)
library(ggplot2)

#######################################################

# read in data in spreadsheet
Ct_data <- read_excel("Cts_WGs.xlsx")

# set x axis order
Ct_data$WG <- factor(Ct_data$WG, levels = c("0", "<10", "10"))

# make plot
ggplot() +
  geom_jitter(width = 0.05, height = NULL, data = Ct_data, aes(x = WG, y = PAN_Ct, color = WG, shape = sample_type), size = 4) +
  scale_color_manual(values = c("10" = "gold", "<10" = "darkblue", "0" = "maroon")) + 
  labs(
    title = "BTV Genomes Generated & Corresponding Ct values",
    x = "Number of Segment Sequences",
    y = "BTV Segment 10 Ct",
    shape="Sample Type", 
    colour="Number of Segment\nSequences") + 
  theme_minimal()

# save plot to PDF
ggsave("Ct_Whole_Genome_Values_blue_red_gold_031926.pdf", height=6, width=8, units="in")
