
# A simulation of shared genome segments in simple pedigrees


Code to accompany the blog post "Simulating shared segments between relatives".

Uses AlphaSimR with recombination tracking to follow IBD segments in small
simulated populations.


## Contents

* `R/simulate_inbreeding_segments.R`: main simulation script, runs three
pedigrees, compares segments between two target individuals and genomic
inbreeding in their offspring.

* `R/simulation_functions.R`: helper functions to simulate a pedigree and 
compare individuals.

* `R/summarise_results.R`: take output from simulation and make graphs.

* `pedigrees/`: three small pedigrees (full sibs, half sibs, cousins)
