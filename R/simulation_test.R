
## Plot the pedigrees and run some checks

library(AlphaSimR)
library(dplyr)
library(GeneticsPed)
library(kinship2)
library(purrr)
library(readr)
library(tibble)

source("R/simulation_functions.R")


ped_fullsib <- read_csv("pedigrees/inbreeding_fullsib.txt")
ped_halfsib <- read_csv("pedigrees/inbreeding_halfsib.txt")
ped_cousin <- read_csv("pedigrees/inbreeding_cousin.txt")


## Plot pedigrees

plot_ped <- function(ped) {
  ped_kinship2 <- pedigree(ped$id,
                           ped$father,
                           ped$mother,
                           ped$sex)
  
  plot(ped_kinship2)
}


dir.create("figures")

pdf("figures/ped_fullsib.pdf",
    height = 4,
    width = 4)
plot_ped(ped_fullsib)
dev.off()

pdf("figures/ped_halfsib.pdf",
    height = 4,
    width = 4)
plot_ped(ped_halfsib)
dev.off()

pdf("figures/ped_cousin.pdf",
    height = 4,
    width = 4)
plot_ped(ped_cousin)
dev.off()


## GenticsPed::inbreeding

inbreeding_ped <- function(ped) {
  
  inbreeding(Pedigree(ped))
  
}

print(map(list(ped_fullsib, ped_halfsib, ped_cousin), inbreeding_ped))



## Set up simulation

founders <- runMacs(nInd = 10,
                    nChr = 25)


simparam <- SimParam$new(founders)

simparam$setTrackRec(TRUE)


founderpop <- newPop(founders,
                     simParam = simparam)



pops <- map(list(ped_fullsib, ped_halfsib, ped_cousin),
            function(ped) pedigreeCross(founderPop = founderpop,
                                        id = ped$id,
                                        mother = ped$mother,
                                        father = ped$father,
                                        simParam = simparam))


parent_comp <- compare_individuals(pops[[1]],
                                   c(3, 4),
                                   simparam)

self_comp <- compare_self(pops[[1]],
                          5,
                          simparam)



ibd <- map(pops, pullIbdHaplo, simParam = simparam, chr = 1)

map <- simparam$genMap[[1]]



ind1 <- ibd[[1]][c("3_1", "3_2"),]
ind2 <- ibd[[1]][c("4_1", "4_2"),]

h1 <- as.vector(ind1[1,])
h2 <- as.vector(ind2[1,]) 


compare_individuals_chr(ind1, ind2, map)



