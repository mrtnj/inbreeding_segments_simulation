## Find shared segments between two haplotypes expressed as vectors
## map is a vector of marker positions

compare_haplotypes <- function(h1, h2, map) {
  sharing <- h1 == h2
  
  runs <- rle(sharing)
  end <- cumsum(runs$lengths)
  start <- c(1, end[-length(end)] + 1)
  
  segments <- tibble(start = start,
                     end = end,
                     start_pos = map[start],
                     end_pos = map[end],
                     segment_length = end_pos - start_pos,
                     value = runs$values)
  
  segments[segments$value,]
}



## Find shared segments between two individuals (expressed as
## matrices of haplotypes) for one chromosome

compare_individuals_chr <- function(ind1, ind2, map) {
  
  h1_1 <- as.vector(ind1[1,])
  h1_2 <- as.vector(ind1[2,])
  
  h2_1 <- as.vector(ind2[1,])
  h2_2 <- as.vector(ind2[2,])
  
  sharing1 <- compare_haplotypes(h1_1, h2_1, map)
  sharing2 <- compare_haplotypes(h1_1, h2_2, map)
  sharing3 <- compare_haplotypes(h1_2, h2_1, map)
  sharing4 <- compare_haplotypes(h1_2, h2_2, map)
  
  bind_rows(sharing1, sharing2, sharing3, sharing4)  
}


## Find shared segments between two target individuals in a
## population

compare_individuals <- function(pop,
                                target_individuals,
                                simparam) {
  
  n_chr <- simparam$nChr
  
  ind1_ix <- paste(target_individuals[1], c("_1", "_2"), sep = "")
  ind2_ix <- paste(target_individuals[2], c("_1", "_2"), sep = "")
  
  ibd <- pullIbdHaplo(pop,
                      simParam = simparam)
  
  map <- simparam$genMap
  loci_per_chr <- map_dbl(map, length)
  
  chr_ends <- cumsum(loci_per_chr)
  chr_starts <- c(1, chr_ends[-n_chr] + 1)
  
  results <- vector(mode = "list",
                    length = n_chr)
  
  for (chr_ix in 1:n_chr) {
    ##print(chr_ix)
    
    ind1 <- ibd[ind1_ix, chr_starts[chr_ix]:chr_ends[chr_ix]]
    ind2 <- ibd[ind2_ix, chr_starts[chr_ix]:chr_ends[chr_ix]]
    
    results[[chr_ix]] <- compare_individuals_chr(ind1, ind2, map[[chr_ix]])
    results[[chr_ix]]$chr <- chr_ix
  } 
  
  bind_rows(results)
}


## Find shared segments between the two haplotypes carried by
## a target individual in a population

compare_self <- function(pop,
                         individual,
                         simparam) {
  
  n_chr <- simparam$nChr
  
  ind_ix <- paste(individual, c("_1", "_2"), sep = "")
  
  ibd <- pullIbdHaplo(pop,
                      simParam = simparam)
  
  map <- simparam$genMap
  loci_per_chr <- map_dbl(map, length)
  
  chr_ends <- cumsum(loci_per_chr)
  chr_starts <- c(1, chr_ends[-n_chr] + 1)
  
  results <- vector(mode = "list",
                    length = n_chr)
  
  for (chr_ix in 1:n_chr) {
    ##print(chr_ix)
    
    ind <- ibd[ind_ix, chr_starts[chr_ix]:chr_ends[chr_ix]]
    
    results[[chr_ix]] <- compare_haplotypes(ind[1, ], ind[2, ], map[[chr_ix]])
    results[[chr_ix]]$chr <- chr_ix
  } 
  
  bind_rows(results)
}


## Run the simulation for a pedigree one replicate

simulate_pedigree <- function(ped,
                              target_individuals,
                              focal_individual,
                              founderpop,
                              simparam) {
  pop <- pedigreeCross(founderPop = founderpop,
                       id = ped$id,
                       mother = ped$mother,
                       father = ped$father,
                       simParam = simparam)
  shared_parents <- compare_individuals(pop,
                                        target_individuals,
                                        simparam)
  shared_inbred <- compare_self(pop,
                                focal_individual,
                                simparam)
  list(population = pop,
       shared_segments_parents = shared_parents,
       shared_segments_self_inbred = shared_inbred)
}