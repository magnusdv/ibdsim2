### Decode map 2019: Halldorsson et al ### 
# Data sets S1 and S2 downloaded from 
# https://science.sciencemag.org/content/363/6425/eaau1043/tab-figures-data

library(tidyverse)
library(ibdsim2)

mraw = read_tsv("data-raw/downloads/aau1043_DataS1.gz", comment = "#")
fraw = read_tsv("data-raw/downloads/aau1043_DataS2.gz", comment = "#")

# Merge and fix chromosome sorting
all = mraw %>% 
  full_join(fraw, by=c("Chr", "Begin", "End")) %>% 
  select(-matches("cMper")) %>%
  arrange(suppressWarnings(as.numeric(sub("chr", "", Chr))), Begin) %>%
  mutate(Chr = fct_inorder(Chr)) %>%
  print

# Thin: Preserve all jumps bigger than 0.05 cM in either males or females
thinned = all %>% 
  group_by(Chr) %>%
  mutate(chromEnd = row_number() %in% c(1, n())) %>%
  ungroup() %>%
  mutate(
    incM = cM.x - lag(cM.x),
    incF = cM.y - lag(cM.y),
    upM = incM >= 0.05, 
    upF = incF >= 0.05,
    jumpM = upM | lead(upM),
    jumpF = upF | lead(upF),
  ) %>%
  filter(chromEnd | jumpM | jumpF) %>%
  print

# Chrom-wise:
chroms = thinned %>%
  group_split(Chr) %>%
  map(~tibble(Mb = c(.$Begin[1], .$End)/1e6,
              male = c(0, .$cM.x),
              female = c(0, .$cM.y)))

# Autosomal chromosome maps
maps = chroms[1:22] %>%
  imap(~ibdsim2:::chromMap(male = select(.x, Mb, male), 
                           female = select(.x, Mb, female),
                           chrom = .y))
# X
maps[[23]] = ibdsim2:::chromMap(male = NULL,
                                female = select(chroms[[23]], Mb, female),
                                chrom = "X")


legacy_decode19 = ibdsim2:::genomeMap(maps)

usethis::use_data(legacy_decode19, overwrite = TRUE)


### QC ###
# Check max error
test1 = mraw |> 
  filter(Chr =="chr1") |> 
  mutate(approx = ibdsim2::convertPos(Mb = End/1e6, map = maps[[1]]$male),
         err = cM - approx)
test1 |> pull(err) |> summary()

# Check chrom order
barplot(sapply(legacy_decode19, mapLen), beside = T)
