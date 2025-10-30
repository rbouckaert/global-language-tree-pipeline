# Table S2 generation script

# 0. Load libraries
library(dplyr)
library(here)
library(magrittr)

# 1. Load and summarise ED data -----------------------------------------------
all_ed <- read.csv(here("outputs", "EDGEscores", "ALL_ED_wTREES.csv")) %>%
  mutate(glottocode = if_else(glottocode == "osse1243", "iron1242", glottocode)) %>%
  group_by(glottocode) %>%
  summarise(ED_mean = mean(ED), .groups = "drop")

# 2. Load Glottolog 4.0 families ----------------------------------------------
lg <- read.csv(here("outputs", "Table_S2", "languoid40.csv"))

lg = lg %>%
  transmute(
    glottocode   = if_else(id == "osse1243", "iron1242", id),
    family_id_4  = family_id,
    parent_id)  %>%  
  left_join(lg %>%  #attach family names
              select(id, name) %>% 
              rename(family_id_4 = id, family_name = name),by = "family_id_4")

# 3. Combine ED with Glottolog and flag isolates ------------------------------
lang_ed <- all_ed %>%
  left_join(lg, by = "glottocode") %>%
  add_count(family_id_4, name = "fam_size") %>%
  mutate(
    fam_size         = if_else(family_id_4 == "", 0L, fam_size),
    True_Isolate     = as.integer(family_id_4 == ""),
    Effective_Isolate = as.integer(fam_size == 1)
  )

# 4. Load threat data ----------------------------------------------------------
threat <- read.csv(here("input_data", "processed_threat_data_frame.csv")) %>%
  mutate(glottocode = if_else(glottocode == "osse1243", "iron1242", glottocode))

# 5. Load top EDGE languages----------------------------------------------------
edge_top <- read.csv(here("outputs", "EDGEscores", "NEW_EDGE_top_100.csv")) %>%
  rename(glottocode = Group.1) %>%
  mutate(glottocode = if_else(glottocode == "osse1243", "iron1242", glottocode))

# 6. Merge all sources and compute final flags --------------------------------
table_s2 <- edge_top %>%
  left_join(lang_ed, by = "glottocode") %>%
  left_join(
    threat %>% select(glottocode, status, glottolog_macroarea),
    by = "glottocode"
  ) %>%
  mutate(
    glottolog_macroarea = recode(glottolog_macroarea,
                                 "Australia" = "Oceania",
                                 "Papunesia" = "Oceania"),
    Language_Isolate    = as.integer(True_Isolate == 1 | Effective_Isolate == 1)
  ) %>%
  select(
    name, glottocode, ED, EDGE, EDGErank,
    ED_mean, family_name, family_id_4,fam_size, 
    True_Isolate, Effective_Isolate, Language_Isolate, status, glottolog_macroarea
  )

# 7. Write output ------------------------------------------------------------
write.csv(table_s2,
          here("outputs", "Table_S2", "S2_New.csv"),
          row.names = FALSE)

tables2_final = table_s2 %>% select(-family_name, -family_id_4, -fam_size, -True_Isolate, -Effective_Isolate)
tables2_final %<>% select(name, glottocode,status, ED, EDGE, EDGErank, glottolog_macroarea, Language_Isolate)

write.csv(tables2_final,
          here("outputs", "Table_S2", "TableS2_New_final.csv"),
          row.names = FALSE)

# 8. Quick summary checks ----------------------------------------------------
table(table_s2$Language_Isolate)
table(table_s2$True_Isolate)
table(table_s2$Effective_Isolate)
