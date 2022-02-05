
library(dplyr)
library(ggplot2)
library(purrr)


## Read simulation results

sim_fullsib <- readRDS("simulations/sim_fullsib.Rds")
sim_halfsib <- readRDS("simulations/sim_halfsib.Rds")
sim_cousin <- readRDS("simulations/sim_cousin.Rds")


## Get proportions of the genome shared between selected parents, and the proportion
## shared between the two genome copies in the inbred individual

get_stats <- function(sim) {
  
  parents_segment_length <-map_dbl(sim,
                                   function(x) sum(x$shared_segments_parents$segment_length))
  
  inbred_segment_length <- map_dbl(sim,
                                   function(x) sum(x$shared_segments_self_inbred$segment_length))
  
  tibble(mean_inbred_sharing = mean(inbred_segment_length/25),
         sd_inbred_sharing = sd(inbred_segment_length/25),
         mean_parents_sharing = mean(parents_segment_length/25/2),
         sd_parents_sharing = sd(parents_segment_length/25/2))

}

stats <- rbind(transform(get_stats(sim_fullsib),
                         case = "full-sib"),
               transform(get_stats(sim_halfsib),
                         case = "half-sib"),
               transform(get_stats(sim_cousin),
                          case = "cousin"))



## Write out table of proportion shared

dir.create("tables")


stats$inbred_self_sharing <- paste(signif(stats$mean_inbred_sharing, 2),
                                   " (",
                                   signif(stats$sd_inbred_sharing, 2),
                                   ")",
                                   sep = "")

stats$parent_sharing <- paste(signif(stats$mean_parents_sharing, 2),
                              " (",
                              signif(stats$sd_parents_sharing, 2),
                              ")",
                              sep = "")

write.csv(stats[, c("case", "inbred_self_sharing", "parent_sharing")],
          file = "tables/proportion_sharing.csv",
          row.names = FALSE,
          quote = FALSE)




## Get counts of segments with particular minimum length

get_segment_counts <- function(sim) {

  length_grid <- seq(from = 0, to = 1, by = 0.01) 
  segment_counts <- map_dfr(sim, function(s) {
    map_dfr(length_grid,
            function(l) {
              tibble(l = l,
                     n_shared_segments = sum(s$shared_segments_parents$segment_length > l))
              
            })
  }, .id = "rep")
  
  segment_counts
}

counts <- rbind(transform(get_segment_counts(sim_fullsib),
                          case = "full-sib"),
                transform(get_segment_counts(sim_halfsib),
                          case = "half-sib"),
                transform(get_segment_counts(sim_cousin),
                          case = "cousin"))



## Minimum numbers of segments of at least 20 cM shared between cousins

print(min(filter(counts, case == "cousin" & l == 0.2)$n_shared_segments))


## Plot of sharing

plot_segment_sharing <- qplot(x = l * 100, y = n_shared_segments, group = rep, geom = "line", data = counts) +
  facet_wrap(~ factor(case, levels = c("full-sib", "half-sib", "cousin"))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  xlab("Minimum segment length (cM)") +
  ylab("Shared segments") +
  ggtitle("Number of shared segments between close relatives")



## Make a stats table of the number of shared segments

counts_stats <- summarise(group_by(filter(counts, l %in% c(0.01, 0.1, 0.2, 0.3, 0.4)), case, l),
                          mean_shared = mean(n_shared_segments),
                          sd_shared = sd(n_shared_segments))

counts_stats$pretty_stats <- paste(signif(counts_stats$mean_shared, 2),
                                   " (",
                                   signif(counts_stats$sd_shared, 2),
                                   ")",
                                   sep = "")
counts_stats$l_cM <- paste(counts_stats$l * 100, "cM")


pretty_table <- pivot_wider(counts_stats[, c("case", "l_cM", "pretty_stats")],
                            values_from = c("pretty_stats"),
                            names_from = c("l_cM"))

pretty_table <- pretty_table[match(c("full-sib", "half-sib", "cousin"), pretty_table$case),]


write.csv(pretty_table,
          file = "tables/mean_segment_count.csv",
          row.names = FALSE,
          quote = FALSE)
