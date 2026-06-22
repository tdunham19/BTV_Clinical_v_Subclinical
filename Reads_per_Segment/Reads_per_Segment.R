#######################################################
# Aim 3 Results TD Thesis
# Reads per Segment - Whole Genome (yes/no)
# TD 03.13.26
#######################################################

# load required libraries
library(tidyverse)
library(readxl)
library(ggplot2)

#######################################################

# read in data in spreadsheet
# metadata <- read_excel("A3_reads_per_segment.xlsx", sheet = "Reads Per Segment Clean No C1_3")
metadata <- read_excel("A3_reads_per_segment.xlsx", sheet = "Reads-Segment Clean NoC1-3 YC19")

# turn it into long format
metadata_long <- metadata %>% pivot_longer(cols = starts_with("segment_"),
                                           names_to = "segment",
                                           values_to = "reads_aligned" )
# set x axis order
metadata_long$segment <- factor(metadata_long$segment, levels = c("Segment_1", "Segment_2", "Segment_3",
                                                                  "Segment_4", "Segment_5", "Segment_6",
                                                                  "Segment_7", "Segment_8", "Segment_9",
                                                                  "Segment_10"))

# change labels to be more clear
## yes = obtained all 10 segment sequences
## partial = obtained < 10 segment sequences (3-5/10)
## no = could not obtain a consensus for any segment
metadata_long$WG <- factor(metadata_long$WG,
                 levels = c("Yes", "Partial", "No"),
                 labels = c("10", "<10", "0"))

# make plot
ggplot() +
  geom_jitter(width = 0.08, height = NULL, 
              data = metadata_long, aes(x = segment, y = reads_aligned, color = WG), size = 4) +
  scale_y_log10() +
  scale_color_manual(values = c("10" = "gold", "<10" = "darkblue", "0" = "maroon")) + 
  labs(title = "Number of Reads to Generate Whole Genome", x = "Segment", y = "Reads Aligned (log10)", 
       colour="Number of Segment\nSequences") +
  theme_minimal(base_size = 15) + 
  scale_x_discrete(
    labels = c(
      "Segment_1" = "1",
      "Segment_2" = "2",
      "Segment_3" = "3",
      "Segment_4" = "4",
      "Segment_5" = "5",
      "Segment_6" = "6",
      "Segment_7" = "7",
      "Segment_8" = "8",
      "Segment_9" = "9",
      "Segment_10" = "10"))

# save plot to PDF
# ggsave("Reads_Per_Segment_blue_red_gold_NO_C1_C3_032126.pdf", height=10, width=15, units="in")
ggsave("Reads_Per_Segment_blue_red_gold_NO_C1_C3_Yes_C19_032126.pdf", height=10, width=15, units="in")
