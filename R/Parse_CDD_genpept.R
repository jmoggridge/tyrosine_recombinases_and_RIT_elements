## Part 2

library(tidyverse)
# Conserved domains specific proteins...
gp <- './data/CDD/sample.gp'
gp <- read_file(gp)

gp_split <- str_split(gp, '\\n//\\n')

str_count(gp_split, '\\')

length(gp_split)
gp_split
