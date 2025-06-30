### Decode map 2019: Halldorsson et al ### 
# Data sets S1 and S2 downloaded from 
# https://science.sciencemag.org/content/363/6425/eaau1043/tab-figures-data

# Revised 2025: Fix MB endpoints


# Physical endpoints ------------------------------------------------------

# Download physical chromosome lengths from Ensembl
L = GenomeInfoDb::getChromInfoFromEnsembl("GRCh38", release = "114") |> 
  subset(name %in% c(1:22, "X"))
ord = order(suppressWarnings(as.numeric(L$name)))
BPLEN = setNames(L$length, L$name)[ord]
BPLEN



# Decode map (raw) ----------------------------------------------------------

library(tidyverse)
options(tibble.print_min = 10)

mraw = read_tsv("data-raw/downloads/aau1043_DataS1.gz", comment = "#")
fraw = read_tsv("data-raw/downloads/aau1043_DataS2.gz", comment = "#")

# Merge and fix chromosome sorting
all = mraw |> 
  full_join(fraw, by = c("Chr", "Begin", "End")) |>
  rename(male = `cM.x`, female = `cM.y`) |> 
  select(-matches("cMper")) |>
  mutate(Chr = sub("chr", "", Chr)) |>
  arrange(suppressWarnings(as.numeric(Chr)), Begin) |>
  print()


# Thinning ----------------------------------------------------------------

# Ramer–Douglas–Peucker algorithm (extended to two y-values):
# Thins (x, y1, y2) to a minimal subset so both y1 and y2 can be linearly 
# interpolated from x within epsilon at all positions.
# Returns indices to keep.
rdp2 = function(x, y1, y2, epsilon) {
  
  # Avoid error for X
  if(is.na(y1[3]))
    y1 = rep(0, length(x))
  
  # Initialise with endpoints
  keep = c(1, length(x))
  
  repeat {
    # Linear interpolation from current subset for both y1 and y2
    interp1 = approx(x[keep], y1[keep], x)$y
    interp2 = approx(x[keep], y2[keep], x)$y
    
    # Max deviation at each point
    dev = pmax(abs(y1 - interp1), abs(y2 - interp2))
    mx = which.max(dev)
    
    # Stop if within tolerance
    if(dev[mx] <= epsilon) 
      break 
    
    # Otherwise, add point with largest deviation
    keep = sort.default(c(keep, mx))     
  }  
  message(sprintf("%d -> %d", length(x), length(keep)))
  keep
}

# Process each chrom separately
maps = list()

for(i in c(1:22, "X")) {
  chri = filter(all, Chr == i)
  ni = nrow(chri)
  
  mp0 = tibble(
    Mb = c(chri$Begin[1], chri$End)/1e6,
    male = c(0, chri$male),
    female = c(0, chri$female)
  ) 
  
  # Thin
  mp = mp0 |> slice(rdp2(Mb, male, female, epsilon = 0.2))
  
  # Add physical endpoints
  mpExt = mp |> 
    add_row(Mb = 0, male = 0, female = 0, .before = 1) |> 
    add_row(Mb = BPLEN[[i]]/1e6, 
            male = chri$male[ni], 
            female = chri$female[ni])

  if(i == "X")
    mpExt$male = NA_real_
  maps[[i]] = mpExt
}

# Create genomeMap as list of chromMap's
decode19 = lapply(1:22, function(i) {
  mp = maps[[i]]
  ibdsim2:::chromMap(
    male   = mp[c("Mb", "male")],
    female = mp[c("Mb", "female")],
    chrom  =  i)
})

# X
decode19[[23]] = ibdsim2:::chromMap(chrom = "X", male = NULL, 
                          female = maps[[23]][c("Mb", "female")])

decode19 = ibdsim2:::genomeMap(decode19)

usethis::use_data(decode19, internal = TRUE, overwrite = TRUE)



# QC tests: RMS ---------------------------------------------------------------


rms = function(orig, thinned) {
  fitM = approx(thinned$Mb, thinned$male, orig$Mb)$y
  fitF = approx(thinned$Mb, thinned$female, orig$Mb)$y
  c(male = sqrt(mean((orig$male - fitM)^2)),
    female = sqrt(mean((orig$female - fitF)^2)))
}

# RMS previous map

RMS1 = sapply(1:22, function(i) {print(i)
  orig = all |> filter(Chr == i) |> 
    mutate(Mb = End/1e6) |> select(Mb, male, female)
  mp = ibdsim2::loadMap("decode19_legacy", chrom = i)[[1]]
  thinned = tibble(Mb = mp$male[,1], male = mp$male[,2], female = mp$female[,2])
  round(c(chrom = i, rms(orig, thinned)),3)
})

# RMS new map
RMS2 = sapply(1:22, function(i) {print(i)
  orig = all |> filter(Chr == i) |> 
    mutate(Mb = End/1e6) |> select(Mb, male, female)
  mp = ibdsim2::loadMap("decode19", chrom = i)[[1]]
  thinned = orig |> slice(rdp2(Mb, male, female, epsilon = 0.2))
  round(c(chrom = i, rms(orig, thinned)),3)
})
