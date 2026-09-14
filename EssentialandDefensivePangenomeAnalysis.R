# Pangenome guide-match analysis
#
# Inputs:
#   BLAST result files with the columns requested by the
#   companion batch-search script.
# Outputs:
#   A data frame and plot showing the number and fraction of genomes with
#   at least one perfect guide-sequence match. The plot can be exported with
#   ggsave() or inspected interactively.


library(tidyverse)
library(janitor)
library(reshape2)


num_genomes <- 16601 # Number of SA genomes in database

### Import Data
ST_output1 <- read_tsv("guide_blast_results_1.tsv")
ST_output2 <- read_tsv("guide_blast_results_2.tsv")
ST_output3 <- read_tsv("guide_blast_results_3.tsv")

full_df <- bind_rows(ST_output1, ST_output2, ST_output3) # Combine data

# Add column names
col_names <- c("query_id", "subject_id", "pct_identity",
               "align_length", "mismatches", "gap_opens",
               "q_start", "q_end", "s_start", "s_end",
               "value", "bit_score")
colnames(full_df) <- col_names

ST_seed_coverage <- full_df %>%
  select(. %in% col_names) %>%
  filter(align_length == 16 & mismatches == 0) %>% # Filter to perfect matches
  group_by(query_id) %>% # Group by ST guide
  # Summarize the number of unique strains with perfect matches
  summarize(strains_hit = n_distinct(subject_id),
            fraction_hit = strains_hit / num_genomes)

ST_seed_coverage$query_id <- factor(ST_seed_coverage$query_id,
                                    levels = ST_seed_coverage$query_id[order(ST_seed_coverage$strains_hit, decreasing = TRUE)]) # Arrange guides in decreasing order

### Plot everything
perfect_matches <- ggplot(ST_seed_coverage, aes(x = query_id, y = strains_hit)) +
  # The main line
  geom_hline(yintercept = num_genomes, linetype = "dashed", color = "gray50") +
  geom_point(size = 3) + 
  # The dashed "Total N Strains" reference line
  annotate("text", x = 30, y = num_genomes, label = "Total Strains = 16,601",
           vjust = -0.5, hjust = 1, color = "gray50", size = 5) +
  # Formatting to match the sketch
  # scale_x_continuous(breaks = thresholds, name = "N Mismatches") +
  scale_y_continuous(name = "Number Strains", limits = c(0, num_genomes * 1.1)) +
  labs(x = "Seed Sequence", title = expression("Strains with" >= "1 Perfect Seed Match")) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(family = "Arial", size = 13,
                              color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(angle = 0),
    axis.title = element_text(family = "Arial", size = 14),
    plot.title = element_text(face = "bold", family = "Arial", size = 16)
  )
perfect_matches # View graph

# Data in ST_seed_coverage was manually exported
# to Graphpad Prism for final figure.
