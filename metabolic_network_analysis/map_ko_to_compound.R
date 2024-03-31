###############################
### Map KEGG KO to compound ###
###############################

# --- (0) Set up ---
# Please set local path before running the code
ko_reaction.list_path <- "set_path"
reaction_mapformula.list <- "set_path"
output_path <- "set_path"

# --- (1) KO to reaction ---

ko_reaction <- read_delim(ko_reaction.list_path, 
           delim = "\t", escape_double = FALSE, 
           col_names = FALSE, trim_ws = TRUE)%>%
  rename(KO = X1, RXN = X2)
ko_reaction$RXN <- gsub("rn:","",ko_reaction$RXN)
ko_reaction$KO <- str_remove(string = ko_reaction$KO, pattern = "ko:")

# --- (2) Reaction to compounds ---

reaction_mapformula <- read_delim(reaction_mapformula.list, 
                                  ":", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
names(reaction_mapformula) <- c("RXN", "PWY", "Formula")
reaction_mapformula <- reaction_mapformula %>%
  select(-PWY) %>%
  distinct()

## Parse formulas - separately for reversible and directional reactions
reaction_mapformula_Producing <- reaction_mapformula %>% 
  filter(grepl("\\s=>", Formula, perl = TRUE)) %>%
  tidyr::separate(col = Formula, into = c("Consuming", "Producing"), sep = "\\s=>") %>%
  mutate(Consuming = trimws(Consuming), Producing = trimws(Producing))%>%
  tidyr::separate_rows(Consuming, sep = " \\+ ")%>%
  tidyr::separate_rows(Producing, sep = " \\+ ")

reaction_mapformula_Consuming <- reaction_mapformula %>% 
  filter(grepl("<=\\s", Formula)) %>%
  tidyr::separate(col = Formula, into = c("Producing", "Consuming"), sep = "<=\\s") %>%
  mutate(Consuming = trimws(Consuming), Producing = trimws(Producing))%>%
  tidyr::separate_rows(Consuming, sep = " \\+ ")%>%
  tidyr::separate_rows(Producing, sep = " \\+ ")

reaction_mapformula_rev <- reaction_mapformula %>% 
  filter(grepl("<=>", Formula)) %>%
  tidyr::separate(col = Formula, into = c("Consuming", "Producing"), sep = "<=>") %>%
  mutate(Consuming = trimws(Consuming), Producing = trimws(Producing)) %>%
  tidyr::separate_rows(Consuming, sep = " \\+ ")%>%
  tidyr::separate_rows(Producing, sep = " \\+ ")

reaction_mapformula_rev_opposite_direction <- reaction_mapformula_rev%>%
  rename(Consuming= Producing, Producing = Consuming) # we change swipe between the columns because the reactions are reversible

rxn_compound <- bind_rows(reaction_mapformula_Producing, reaction_mapformula_Consuming, reaction_mapformula_rev, reaction_mapformula_rev_opposite_direction)%>%
  unique()

# --- (3) Merge all tables ---

ko_compound <- left_join(ko_reaction,rxn_compound , by = "RXN")%>%
  drop_na()%>%
  unique()

# --- (4) save and remove unnecessary data ---
saveRDS(ko_compound, str_c(output_path, "ko_compound.rds"))
remove(ko_reaction, reaction_mapformula, reaction_mapformula_Producing, reaction_mapformula_Consuming, reaction_mapformula_rev, reaction_mapformula_rev_opposite_direction, rxn_compound)
  