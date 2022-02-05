
## Simulations of length of shared segments from inbreeding

library(AlphaSimR)
library(dplyr)
library(purrr)
library(readr)
library(tibble)

source("R/simulation_functions.R")

## Set up simulation

founders <- runMacs(nInd = 10,
                    nChr = 25)


simparam <- SimParam$new(founders)

simparam$setTrackRec(TRUE)


founderpop <- newPop(founders,
                     simParam = simparam)


## Pedigrees

ped_fullsib <- read_csv("pedigrees/inbreeding_fullsib.txt")
ped_halfsib <- read_csv("pedigrees/inbreeding_halfsib.txt")
ped_cousin <- read_csv("pedigrees/inbreeding_cousin.txt")


## Target individuals for each pedigree

parents <- list(c(3, 4), c(4, 5), c(7, 8))
inbred <- c(5, 6, 9)


## Simulations

sim_fullsib <- replicate(100,
                         simulate_pedigree(ped = ped_fullsib,
                                           target_individuals = parents[[1]],
                                           focal_individual = inbred[1],
                                           founderpop = founderpop,
                                           simparam = simparam),
                         simplify = FALSE)
             
sim_halfsib <- replicate(100,
                         simulate_pedigree(ped = ped_halfsib,
                                           target_individuals = parents[[2]],
                                           focal_individual = inbred[2],
                                           founderpop = founderpop,
                                           simparam = simparam),
                         simplify = FALSE)

sim_cousin <- replicate(100,
                        simulate_pedigree(ped = ped_cousin,
                                          target_individuals = parents[[3]],
                                          focal_individual = inbred[3],
                                          founderpop = founderpop,
                                          simparam = simparam),
                        simplify = FALSE)



## Save

dir.create("simulations")

saveRDS(sim_fullsib,
        file = "simulations/sim_fullsib.Rds")

saveRDS(sim_halfsib,
        file = "simulations/sim_halfsib.Rds")

saveRDS(sim_cousin,
        file = "simulations/sim_cousin.Rds")
