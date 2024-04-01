#' Differentially Dispersing Taxa Analysis
#' 
#' 
#' This script performs analysis to identify differentially dispersing taxa pairs. 
#' It calculates the relative dispersal ratio (RDR) for each pair of taxa and conducts 
#' statistical tests to determine significant differences in dispersal behavior.
#' 
#' @title Differentially Dispersing Taxa Analysis
#' @param data_file Path to the data file containing experimental outcomes.
#'  #' The `fmt_outcome` file should contain four columns: 
#' - `taxon`: Name or identifier of the taxon.
#' - `recipient_post_fmt`: Value representing the recipient community post-treatment.
#' - `donor_pre_fmt`: Value representing the donor community pre-treatment.
#' - `taxon_id`: Identifier for the taxon (can be numeric or alphanumeric).
#' - `experiment_id`: Identifier for the experiment.
#' @param metadata_file Path to the metadata file containing additional subject information.
#' @param min_subjects Minimum number of subjects required for analysis.
#' @param alpha Significance level for determining statistical significance.
#' @param statistical_test Method for statistical testing ("wilcoxon" or "lmer").
#' - 'wilcoxon': Wilcoxon rank-sum test. Can be used in case we dont need to correct for confounders (etc, multiple samples from the same donor/subject).
#' - 'lmer': Linear mixed-effects model. Can be used when we need to correct for confounders (e.g., multiple samples from the same donor/subject).
#' @param stat_signif_pairs_output_file Path to the output file for storing differently dispersing pairs results.
#' @param rdr_output_file Path to the output file for storing RDR data.
#' @return A list of differentially dispersing pairs and RDR data.
#' 
# Set up 

##Settings

# Data file paths
data_file <- "fmt_outcome.csv"
metadata_file <- "donor_metadata.csv"

# Analysis parameters
min_subjects <- 10 
alpha <- 0.1  

# Statistical test method
statistical_test <- "wilcoxon"  

# Output file paths
stat_signif_pairs_output_file <- "stat_signif_pairs.csv"
rdr_output_file <- "rdr_data.csv"

## Check and install required packages
required_packages <- c("stringr", "tidyverse", "forcats", "nlme")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(stringr) # stringr_1.4.0
library(tidyverse) # tidyverse_2.0.0    
library(forcats) # forcats_1.0.0      
library(nlme) # nlme_3.1-157
library(readr) # readr_1.4.0

## Define function to find statistically significant pairs
#' Find statistically significant pairs
#' 
#' This function calculates statistical significance for each pair of taxa.
#' 
#' @param .x Data frame containing relevant variables.
#' @param method Method for significance testing ("lmer" or "wilcoxon").
#' @return A tibble with statistical results.
#' @export

find_stat_signif_pairs <- function(.x, method = "wilcoxon"){
  
  if(method == "lmer"){  
    rme_results <- lmerTest::lmer(rdr ~ 1 + (1 | subject_id) + (1 | donor_id), data = .x,
                                  control=lmerControl(check.nobs.vs.nlev = "ignore",
                                                      check.nobs.vs.rankZ = "ignore",
                                                      check.nobs.vs.nRE="ignore"))
    
    results <- summary(rme_results)$coefficients%>%
      as_tibble(.)%>%
      rename(p_value = "Pr(>|t|)", t_value = "t value")%>%
      mutate(pair_id = .x$pair_id[[1]])%>%
      rename(mean = Estimate)
    
  }else if (method == "wilcoxon"){
    results <- wilcox.test(.x$rdr, data = .x, mu = 0)%>%
      tidy(.)%>%
      mutate(pair_id = .x$pair_id[[1]])%>%
      rename(p_value = p.value)%>%
      mutate(mean  = mean(.x$rdr))%>%
      select(-method, -alternative)
    
  }else{
    results <- "method not implemented"
  }
  
  return(results)
}

## Import data 
fmt_outcome <- read_rds(data_file)

# Main
## Calculate the relative dispersal ratio (RDR) for each pair of taxa

taxa_two <- fmt_outcome%>%
  rename(taxa_two = taxon,
         recipient_post_fmt_two =recipient_post_fmt ,
         donor_pre_fmt_two = donor_pre_fmt,
         taxon_id_two = taxon_id)

rdr_data <- fmt_outcome%>%
    rename(taxa_one = taxon,
           recipient_post_fmt_one = recipient_post_fmt ,
           donor_pre_fmt_one = donor_pre_fmt,
           taxon_id_one = taxon_id)%>%
    left_join(taxa_two, by = "experiment_id")%>%
    filter(taxon_id_one>taxon_id_two)%>%
    rowwise()%>%
    mutate(rdr  = log10(recipient_post_fmt_one/recipient_post_fmt_two)/(donor_pre_fmt_one /donor_pre_fmt_two),
           pair_id = str_c(taxon_id_one, "_", taxon_id_two))%>%
    ungroup()%>%
    separate(experiment_id, into = c("subject_id", "time_point_post"), sep = "_d")


##Filter pairs with low number of samples
pairs_subject_count <- rdr_data%>%
  select(pair_id, subject_id)%>%
  unique()%>%
  count(pair_id)%>%
  rename(n_subjects = n)

pairs_for_analysis <- pairs_subject_count%>%
  filter(n_subjects>9)%>%
  pull(pair_id)

## Find differentially dispersing pairs 

test_stat_signif_pairs<- rdr_data%>%
  filter(pair_id %in% pairs_for_analysis)%>%
  left_join(.,donor_metadata,  by = "subject_id")%>%
  mutate(time_point_post = as.numeric(time_point_post),
         subject_id = as.factor(subject_id))%>%
  #run linear mixed effect model for every pair
  group_by(pair_id)%>%
  group_map(.keep = TRUE, ~{find_stat_signif_pairs(.x, method = "wilcoxon")})%>%
  bind_rows()%>%
  # Define the dominant pair 
  separate(pair_id, into = c("taxa_1", "taxa_2"), sep = "_")%>%
  mutate(taxa_1_tmp = if_else(mean<0, true = taxa_2, false = taxa_1),
         taxa_2 = if_else(mean<0, true = taxa_1, false = taxa_2),
         taxa_1 = taxa_1_tmp,
         mean = if_else(mean<0, true = -mean, false = mean),
         pair_id = str_c(taxa_1, "_", taxa_2))%>%
  select(-taxa_1_tmp)

stat_signif_pairs <- mutate(test_stat_signif_pairs, fdr_p_value = p.adjust(p_value, method = "fdr"))%>%
  filter(fdr_p_value<0.1)

## Save results
write_csv(test_stat_signif_pairs, stat_signif_pairs_output_file)
write_csv(rdr_data, rdr_output_file)
