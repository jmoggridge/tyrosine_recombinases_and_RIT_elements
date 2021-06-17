library(tidyverse)
library(tidymodels)


## downsample to size of smallest class
set.seed(1)

# smallest class size
smallest <- full_data |> group_by(subfamily) |> count() |> pull(n) |> min()
# downsample data
df <- full_data |> 
  group_by(subfamily) |> 
  sample_n(size =  smallest, replace = F) |> 
  ungroup()

glimpse(df)
rm(smallest, integrases, non_integrases, full_data)

