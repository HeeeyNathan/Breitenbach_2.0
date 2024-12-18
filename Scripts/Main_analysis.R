# ========== META DATA ==========
# title:        Breitenbach 2.0
# author:       Viktor Baranov & Nathan Jay Baker
# date created: 2024/12/18
# DOI:          XXXX

# ========== LOAD PACKAGES ==========
library(tidyverse)
library(patchwork)
library(janitor)
library(SpadeR)
library(vegan)
library(codyn)
library(mobr)
library(iNEXT)
library(readxl)
library(writexl)

# ========== LOAD DATA ==========
diptera <- read_xlsx("Data/Breitenbach_community_data_17.12.2024.xlsx") |> # read in excel spreadsheet
  as_tibble() |>                                                           # Make tibble table
  clean_names() |>                                                         # Homogenize column names
  rename_with(~ gsub("^x", "", .x), .cols = matches("^x[0-9]")) |>         # Change column names
  filter(!str_detect(original_name, "Summe")) |>                           # delete rows containing "summe"
  filter(order == "Diptera") |>                                            # Keeps only Dipteran taxa
  select(where(~ !is.numeric(.x) || sum(.x) != 0)) |>                      # Removes empty (0) columns
  mutate(trap_new = case_when(                                             # Homogenize trap names
    trap %in% c('Haus 0') ~ 'O',
    trap %in% c('Haus A', 'A/I') ~ 'A',
    trap %in% c('Haus I') ~ 'I',
    trap %in% c('Haus B', 'B / II', 'B-II-X', 'B/II', 'Haus B/II') ~ 'B',
    trap %in% c('Haus C-IV', 'C / IV', 'C-IV', 'C/IV', 'Haus C', 'Haus C/IV') ~ 'C',
    trap %in% c('Haus III') ~ 'III',
    trap %in% c('Haus D') ~ 'D',
    trap %in% c('Haus E - V', 'E / V', 'E-V', 'E/V', 'Haus E/V') ~ 'E',
    trap %in% c('Haus F - VI', 'F / VI', 'F-VI', 'F-VII', 'F/VI', 'Haus F/VI') ~ 'F',
    trap %in% c('Haus G', 'G-VII', 'G/ VII', 'G/VII', 'Haus G/VII') ~ 'G',
    trap %in% c('Quelle') ~ 'Source',
    TRUE ~ trap)) |>
  select(-c(taxa_id:trap)) |>                                              # Remove unnecessary columns
  mutate(sp_name  = validated_name,                                        # Change column names
         trap     = trap_new,                                              # Change column names
         trap     = factor(trap)) |>                                       # make trap a factor
  arrange(trap, family, sp_name) |>                                        # Sort data
  select(-c(validated_name, trap_new, order)) |>                           # Remove unnecessary columns
  select(trap, sp_name, family, everything()) |>                           # rearrange columns
  print()
diptera |>
  select(where(is.numeric)) |>
  unlist() |>
  sum(na.rm = TRUE)      # check combined abundance of all species

ept <- read_xlsx("Data/Breitenbach_community_data_17.12.2024.xlsx") |>     # read in excel spreadsheet
  as_tibble() |>                                                           # Make tibble table
  clean_names() |>                                                         # Homogenize column names
  rename_with(~ gsub("^x", "", .x), .cols = matches("^x[0-9]")) |>         # Change column names
  filter(!str_detect(original_name, "Summe")) |>                           # delete rows containing "summe"
  filter(order %in% c("Ephemeroptera", "Plecoptera", "Trichoptera")) |>    # Keeps only EPT taxa
  select(where(~ !is.numeric(.x) || sum(.x) != 0)) |>                      # Removes empty (0) columns
  mutate(trap_new = case_when(                                             # Homogenize trap names
    trap %in% c('Haus 0') ~ 'O',
    trap %in% c('Haus A', 'A/I') ~ 'A',
    trap %in% c('Haus I') ~ 'I',
    trap %in% c('Haus B', 'B / II', 'B-II-X', 'B/II', 'Haus B/II') ~ 'B',
    trap %in% c('Haus C-IV', 'C / IV', 'C-IV', 'C/IV', 'Haus C', 'Haus C/IV') ~ 'C',
    trap %in% c('Haus III') ~ 'III',
    trap %in% c('Haus D') ~ 'D',
    trap %in% c('Haus E - V', 'E / V', 'E-V', 'E/V', 'Haus E/V') ~ 'E',
    trap %in% c('Haus F - VI', 'F / VI', 'F-VI', 'F-VII', 'F/VI', 'Haus F/VI') ~ 'F',
    trap %in% c('Haus G', 'G-VII', 'G/ VII', 'G/VII', 'Haus G/VII') ~ 'G',
    trap %in% c('Quelle') ~ 'Source',
    TRUE ~ trap)) |>
  filter(!str_detect(trap_new, "D")) |>                                    # delete rows for trap D (only sampled 1983)
  select(-c(taxa_id:trap)) |>                                              # Remove unnecessary columns
  mutate(sp_name  = validated_name,                                        # Change column names
         trap     = trap_new,                                              # Change column names
         trap     = factor(trap)) |>                                       # make trap a factor
  arrange(trap, family, sp_name) |>                                        # Sort data
  select(-c(validated_name, trap_new, order)) |>                           # Remove unnecessary columns
  select(trap, sp_name, family, everything()) |>                           # rearrange columns
  print()
ept |>
  select(where(is.numeric)) |>
  unlist() |>
  sum(na.rm = TRUE)                                                        # check combined abundance of all species

# ========== REMOVE DUPLICATES ==========
diptera_clean <- diptera |>
  group_by(trap,
           sp_name,
           family) |>                                                      # Grouping by trap, sp_name, and family
  summarise(across(where(is.numeric),
                   ~sum(.x, na.rm = TRUE)),
                   .groups = "drop") |>
  print()
diptera_clean |>
  select(where(is.numeric)) |>
  unlist() |>
  sum(na.rm = TRUE)                                                        # check combined abundance of all species

ept_clean <- ept |>
  group_by(trap,
           sp_name,
           family) |>                                                      # Grouping by trap, sp_name, and family
  summarise(across(where(is.numeric),
                   ~sum(.x, na.rm = TRUE)),
                   .groups = "drop") |>
  print()
ept_clean |>
  select(where(is.numeric)) |>
  unlist() |>
  sum(na.rm = TRUE)

# ========== TRANSFORM DATA ==========
diptera_long <- diptera_clean |>
  pivot_longer(
    cols = starts_with("19") | starts_with("20"),                          # Specify the year column to pivot
    names_to = "year",                                                     # Name column that will hold year values
    values_to = "abundance") |>                                            # Name column that will hold the counts
  mutate(trap_code = paste(trap, year, sep = "_")) |>                      # create new column with trap_code (trap name + year)
  mutate(year = as.numeric(year)) |>                                       # make year variable numeric
  filter(abundance != 0) |>                                                # Remove rows where count is 0
  arrange(trap, year)                                                      # Sort data

ept_long <- ept_clean |>
  pivot_longer(
    cols = starts_with("19") | starts_with("20"),                          # Specify the year column to pivot
    names_to = "year",                                                     # Name column that will hold year values
    values_to = "abundance") |>                                            # Name column that will hold the counts
  mutate(trap_code = paste(trap, year, sep = "_")) |>                      # create new column with trap_code (trap name + year)
  mutate(year = as.numeric(year)) |>                                       # make year variable numeric
  filter(abundance != 0) |>                                                # Remove rows where count is 0
  arrange(trap, year)                                                      # Sort data

# ========== AGGREGATE DATA BY TRAP ==========
diptera_agg <- diptera_long |>
  group_by(sp_name,
           trap) |>                                                        # Group by sp_name and trap
  summarise(abundance = sum(abundance),
            .groups = "drop") |>                                           # Summarise with the sum of "abundance"
  arrange(trap) |>                                                         # Sort data
  print()
sum(diptera_agg$abundance)                                                 # check combined abundance of all species

ept_agg <- ept_long |>
  group_by(sp_name,
           trap) |>                                                        # Group by sp_name and trap
  summarise(abundance = sum(abundance),
            .groups = "drop") |>                                           # Summarise with the sum of "abundance"
  arrange(trap) |>                                                         # Sort data
  print()
sum(ept_agg$abundance)                                                     # check combined abundance of all species

# ========== PIVOT DATA ==========
diptera_wide <- diptera_agg |>
  pivot_wider(names_from  = trap,
              values_from = abundance,
              values_fill = 0,
              values_fn   = sum) |>
  arrange(sp_name) |>                                                      # Sort data
  print()

ept_wide <- ept_agg |>
  pivot_wider(names_from  = trap,
              values_from = abundance,
              values_fill = 0,
              values_fn   = sum) |>
  arrange(sp_name) |>                                                      # Sort data
  print()

# ========== START OF ANALYSIS ==========
## ----------- Estimation of species richness in a community (Chao) -----------
chaoA_dipt      <- ChaoSpecies(diptera_wide$A ,"abundance", k = 2, conf = 0.95)      # Diptera trap A
chaoB_dipt      <- ChaoSpecies(diptera_wide$B ,"abundance", k = 2, conf = 0.95)      # Diptera trap B
chaoC_dipt      <- ChaoSpecies(diptera_wide$C ,"abundance", k = 2, conf = 0.95)      # Diptera trap C
chaoD_dipt      <- ChaoSpecies(diptera_wide$D ,"abundance", k = 2, conf = 0.95)      # Diptera trap D
chaoE_dipt      <- ChaoSpecies(diptera_wide$E ,"abundance", k = 2, conf = 0.95)      # Diptera trap E
chaoF_dipt      <- ChaoSpecies(diptera_wide$F ,"abundance", k = 2, conf = 0.95)      # Diptera trap F
chaoG_dipt      <- ChaoSpecies(diptera_wide$G ,"abundance", k = 2, conf = 0.95)      # Diptera trap G
chaoI_dipt      <- ChaoSpecies(diptera_wide$I ,"abundance", k = 2, conf = 0.95)      # Diptera trap I
chaoIII_dipt    <- ChaoSpecies(diptera_wide$III ,"abundance", k = 2, conf = 0.95)    # Diptera trap III

chaoA_ept      <- ChaoSpecies(ept_wide$A ,"abundance", k = 2, conf = 0.95)           # ept trap A
chaoB_ept      <- ChaoSpecies(ept_wide$B ,"abundance", k = 2, conf = 0.95)           # ept trap B
chaoC_ept      <- ChaoSpecies(ept_wide$C ,"abundance", k = 2, conf = 0.95)           # ept trap C
chaoE_ept      <- ChaoSpecies(ept_wide$E ,"abundance", k = 2, conf = 0.95)           # ept trap E
chaoF_ept      <- ChaoSpecies(ept_wide$F ,"abundance", k = 2, conf = 0.95)           # ept trap F
chaoG_ept      <- ChaoSpecies(ept_wide$G ,"abundance", k = 2, conf = 0.95)           # ept trap G
chaoI_ept      <- ChaoSpecies(ept_wide$I ,"abundance", k = 2, conf = 0.95)           # ept trap I
chaoIII_ept    <- ChaoSpecies(ept_wide$III ,"abundance", k = 2, conf = 0.95)         # ept trap III
chaoO_ept      <- ChaoSpecies(ept_wide$O ,"abundance", k = 2, conf = 0.95)           # ept trap O

## ----------- Univariate statistics -----------
diptera_long$ro.ab <- round(diptera_long$abundance, digits = 0)                      # Round abundance values to whole numbers for rarefaction analysis
diptera_TD <- NULL                                                                   # Initialize empty data frame to store results

for (i in unique(diptera_long$trap)) {                                               # Loop through each unique trap in the dataset
  sub <- diptera_long[diptera_long$trap == i, ]                                      # Create subset for current trap and reshape data from long to wide format
  sub_m <- sub |>
    select(trap_code, sp_name, abundance) |>
    pivot_wider(names_from = sp_name, values_from = abundance, values_fill = 0)
  sub_ta <- sub_m[, -1]                                                              # Remove trap_code column for calculations

  # Calculate diversity indices
  SppRich <- specnumber(sub_ta)                                                      # Species richness (total number of species)
  Simp <- diversity(sub_ta, index = "simpson")                                       # Simpson's diversity (probability two random individuals are different species)
  Shan <- diversity(sub_ta, index = "shannon")                                       # Shannon's diversity (accounts for both abundance and evenness)
  EvenJ <- Shan / log(SppRich)                                                       # Pielou's evenness (how close in numbers each species is)
  E10 <- Shan / SppRich                                                              # Shannon's evenness (alternative evenness measure)
  Abund <- rowSums(sub_ta)                                                           # Total abundance (sum of all individuals)
  S_PIE <- calc_PIE(sub_ta, ENS = TRUE)                                              # Effective number of common species

  # Calculate species turnover metrics between years
  DATA1_Turnover <- codyn::turnover(sub,
    time.var = "year",
    species.var = "sp_name",
    abundance.var = "abundance",
    replicate.var = NA,
    metric = "total")                                                                # Total turnover (appearances + disappearances)
  Turnover <- c("NA", DATA1_Turnover$total)                                          # Add NA for first year

  DATA1_Turnover_app <- codyn::turnover(sub,                                         # Species appearances only
    time.var = "year",
    species.var = "sp_name",
    abundance.var = "abundance",
    replicate.var = NA,
    metric = "appearance")
  Turnover_app <- c("NA", DATA1_Turnover_app$appearance)

  DATA1_Turnover_disapp <- codyn::turnover(sub,                                      # Species disappearances only
    time.var = "year",
    species.var = "sp_name",
    abundance.var = "abundance",
    replicate.var = NA,
    metric = "disappearance")
  Turnover_disapp <- c("NA", DATA1_Turnover_disapp$disappearance)

  # Prepare data for rarefaction analysis
  sub_m_r <- sub |>
    select(trap_code, sp_name, ro.ab) |>
    pivot_wider(names_from = sp_name,
      values_from = ro.ab,
      values_fill = 0)                                                               # Create matrix with rounded abundances
  sub_ta_r <- sub_m_r[, -1]                                                          # Remove trap_code column

  # Calculate rarefied species richness
  rare.SppRich <- if (min(rowSums(sub_ta_r)) > 10) {
    rarefy(sub_ta_r, sample = min(rowSums(sub_ta_r)))
  } else {
    rarefy(sub_ta_r, sample = 10)
  }                                                                                  # If minimum sample size > 10, use that; otherwise use 10 as minimum

  # Combine all calculated metrics into a data frame
  TD.i <- data.frame(sub_m$trap_code,
    SppRich,
    Simp,
    Shan,
    EvenJ,
    E10,
    Abund,
    S_PIE,
    Turnover,
    Turnover_app,
    Turnover_disapp,
    rare.SppRich)

  # Append results to main data frame
  diptera_TD <- rbind(diptera_TD, TD.i)

  # Clean up temporary variables to avoid conflicts in next iteration
  rm(TD.i, sub_m, sub_ta, sub, SppRich, Simp, Shan, EvenJ, E10, Abund, S_PIE,
    DATA1_Turnover, Turnover, DATA1_Turnover_app, Turnover_app,
    DATA1_Turnover_disapp, Turnover_disapp, sub_m_r, sub_ta_r, rare.SppRich)
}

# Clean up dataframe
diptera_TD <- diptera_TD |>
  as_tibble() |>                                                                     # Convert to tibble format for better handling
  rename(trap_code = sub_m.trap_code) |>                                             # Fix column name from earlier processing
  mutate(trap = gsub("_.*", "", trap_code),                                          # Remove everything after _ to get trap ID
         year = gsub(".*_", "", trap_code)) |>                                       # Remove everything before _ to get year
  mutate(across(-c(trap_code, trap), as.numeric)) |>                                 # Convert columns to numeric format
  select(trap_code, trap, year, everything()) |>                                     # Reorder columns: trap_code, trap, year first, then rest
  print()

ept_long$ro.ab <- round(ept_long$abundance, digits = 0)                              # Round abundance values to whole numbers for rarefaction analysis
ept_TD <- NULL                                                                       # Initialize empty data frame to store results

for (i in unique(ept_long$trap)) {                                                   # Loop through each unique trap in the dataset
  sub <- ept_long[ept_long$trap == i, ]                                              # Create subset for current trap and reshape data from long to wide format
  sub_m <- sub |>
    select(trap_code, sp_name, abundance) |>
    pivot_wider(names_from = sp_name, values_from = abundance, values_fill = 0)
  sub_ta <- sub_m[, -1]                                                              # Remove trap_code column for calculations

  # Calculate diversity indices
  SppRich <- specnumber(sub_ta)                                                      # Species richness (total number of species)
  Simp <- diversity(sub_ta, index = "simpson")                                       # Simpson's diversity (probability two random individuals are different species)
  Shan <- diversity(sub_ta, index = "shannon")                                       # Shannon's diversity (accounts for both abundance and evenness)
  EvenJ <- Shan / log(SppRich)                                                       # Pielou's evenness (how close in numbers each species is)
  E10 <- Shan / SppRich                                                              # Shannon's evenness (alternative evenness measure)
  Abund <- rowSums(sub_ta)                                                           # Total abundance (sum of all individuals)
  S_PIE <- calc_PIE(sub_ta, ENS = TRUE)                                              # Effective number of common species

  # Calculate species turnover metrics between years
  DATA1_Turnover <- codyn::turnover(sub,
    time.var = "year",
    species.var = "sp_name",
    abundance.var = "abundance",
    replicate.var = NA,
    metric = "total")                                                                # Total turnover (appearances + disappearances)
  Turnover <- c("NA", DATA1_Turnover$total)                                          # Add NA for first year

  DATA1_Turnover_app <- codyn::turnover(sub,                                         # Species appearances only
    time.var = "year",
    species.var = "sp_name",
    abundance.var = "abundance",
    replicate.var = NA,
    metric = "appearance")
  Turnover_app <- c("NA", DATA1_Turnover_app$appearance)

  DATA1_Turnover_disapp <- codyn::turnover(sub,                                      # Species disappearances only
    time.var = "year",
    species.var = "sp_name",
    abundance.var = "abundance",
    replicate.var = NA,
    metric = "disappearance")
  Turnover_disapp <- c("NA", DATA1_Turnover_disapp$disappearance)

  # Prepare data for rarefaction analysis
  sub_m_r <- sub |>
    select(trap_code, sp_name, ro.ab) |>
    pivot_wider(names_from = sp_name,
      values_from = ro.ab,
      values_fill = 0)                                                               # Create matrix with rounded abundances
  sub_ta_r <- sub_m_r[, -1]                                                          # Remove trap_code column

  # Calculate rarefied species richness
  rare.SppRich <- if (min(rowSums(sub_ta_r)) > 10) {
    rarefy(sub_ta_r, sample = min(rowSums(sub_ta_r)))
  } else {
    rarefy(sub_ta_r, sample = 10)
  }                                                                                  # If minimum sample size > 10, use that; otherwise use 10 as minimum

  # Combine all calculated metrics into a data frame
  TD.i <- data.frame(sub_m$trap_code,
    SppRich,
    Simp,
    Shan,
    EvenJ,
    E10,
    Abund,
    S_PIE,
    Turnover,
    Turnover_app,
    Turnover_disapp,
    rare.SppRich)

  # Append results to main data frame
  ept_TD <- rbind(ept_TD, TD.i)

  # Clean up temporary variables to avoid conflicts in next iteration
  rm(TD.i, sub_m, sub_ta, sub, SppRich, Simp, Shan, EvenJ, E10, Abund, S_PIE,
    DATA1_Turnover, Turnover, DATA1_Turnover_app, Turnover_app,
    DATA1_Turnover_disapp, Turnover_disapp, sub_m_r, sub_ta_r, rare.SppRich)
}

# Clean up dataframe
ept_TD <- ept_TD |>
  as_tibble() |>                                                                     # Convert to tibble format for better handling
  rename(trap_code = sub_m.trap_code) |>                                             # Fix column name from earlier processing
  mutate(trap = gsub("_.*", "", trap_code),                                          # Remove everything after _ to get trap ID
         year = gsub(".*_", "", trap_code)) |>                                       # Remove everything before _ to get year
  mutate(across(-c(trap_code, trap), as.numeric)) |>                                 # Convert columns to numeric format
  select(trap_code, trap, year, everything()) |>                                     # Reorder columns: trap_code, trap, year first, then rest
  print()





### ------------ Visualizations ------------
My_theme <- theme(panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, linewidth = 1.25),
                  strip.background = element_rect(fill = "white",
                                                  color = "white", linewidth = 1.25),
                  legend.position = "bottom",
                  text = element_text(size = 16,
                                      family = "gillsans"))

# Turnover by trap by year
# Diptera
ggplot(diptera_TD |>
         group_by(trap) |>
         filter(n() > 3),
       aes(x = year, y = Turnover)) +
  geom_line() +
  geom_smooth(method = "lm") +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.5,
             linetype = "dotted",
             color = "grey50",
             size = 0.5) +
  facet_wrap(~trap) +
  scale_x_continuous(breaks = seq(min(diptera_TD$year), max(diptera_TD$year), by = 3)) +
  scale_y_continuous(limits = c(0.25, max(diptera_TD$Turnover, na.rm = TRUE))) +
  labs(title = "Dipteran species turnover",
       x = "Year",
       y = "Species turnover",
       caption = "Diptera species turnover over time by trap") +
  My_theme +
  My_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.25, "cm"))

ggsave(
  "Plots/turnover_diptera.png",
  plot = last_plot(),
  scale = 1,
  width = 10,
  height = 10,
  units = "in",
  dpi = 450,
  device = "png",
  bg = "white"
)

# EPT
ggplot(ept_TD |>
         group_by(trap) |>
         filter(n() > 3),
       aes(x = year, y = Turnover)) +
  geom_line() +
  geom_smooth(method = "lm") +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.5,
             linetype = "dotted",
             color = "grey50",
             size = 0.5) +
  facet_wrap(~trap) +
  scale_x_continuous(breaks = seq(min(ept_TD$year), max(ept_TD$year), by = 3)) +
  scale_y_continuous(limits = c(0.25, max(ept_TD$Turnover, na.rm = TRUE))) +
  labs(title = "EPT species turnover",
       x = "Year",
       y = "Species turnover",
       caption = "EPT species turnover over time by trap") +
  My_theme +
  My_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.25, "cm"))

ggsave(
  "Plots/turnover_ept.png",
  plot = last_plot(),
  scale = 1,
  width = 10,
  height = 10,
  units = "in",
  dpi = 450,
  device = "png",
  bg = "white"
)

# Abundance by trap by year
# Diptera
ggplot(diptera_TD |>
         group_by(trap) |>
         filter(n() > 3),
       aes(x = year, y = Abund)) +
  geom_line() +
  geom_smooth(method = "lm") +
  # geom_smooth(method = "gam", col = "red") +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.5,
             linetype = "dotted",
             color = "grey50",
             size = 0.5) +
  facet_wrap(~trap, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(diptera_TD$year), max(diptera_TD$year), by = 3)) +
  labs(title = "Diptera abundance",
       x = "Year",
       y = "Abundance",
       caption = "Diptera abundance through time by trap") +
  My_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.25, "cm"))

ggsave(
  "Plots/abundance_diptera.png",
  plot = last_plot(),
  scale = 1,
  width = 10,
  height = 10,
  units = "in",
  dpi = 450,
  device = "png",
  bg = "white"
)

# EPT
ggplot(ept_TD |>
         group_by(trap) |>
         filter(n() > 3),
       aes(x = year, y = Abund)) +
  geom_line() +
  geom_smooth(method = "lm") +
  # geom_smooth(method = "gam", col = "red") +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.5,
             linetype = "dotted",
             color = "grey50",
             size = 0.5) +
  facet_wrap(~trap, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(ept_TD$year), max(ept_TD$year), by = 3)) +
  labs(title = "EPT abundance",
       x = "Year",
       y = "Abundance",
       caption = "EPT abundance through time by trap") +
  My_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.25, "cm"))

ggsave(
  "Plots/abundance_ept.png",
  plot = last_plot(),
  scale = 1,
  width = 10,
  height = 10,
  units = "in",
  dpi = 450,
  device = "png",
  bg = "white"
)

# Species richness by trap by year
# Diptera
ggplot(diptera_TD %>%
         group_by(trap) %>%
         filter(n() > 3),
       aes(x = year, y = SppRich)) +
  geom_line() +
  geom_smooth(method = "lm") +
  geom_point(size = 2) +
  # geom_hline(yintercept = 0.5,
  #            linetype = "dotted",
  #            color = "grey50",
  #            size = 0.5) +
  facet_wrap(~trap, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(diptera_TD$year), max(diptera_TD$year), by = 3)) +
  labs(title = "Diptera species richness",
       x = "Year",
       y = "Species richness",
       caption = "Diptera species richness through time by trap") +
  My_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.25, "cm"))

ggsave(
  "Plots/species_richness_diptera.png",
  plot = last_plot(),
  scale = 1,
  width = 10,
  height = 10,
  units = "in",
  dpi = 450,
  device = "png",
  bg = "white"
)

# EPT
ggplot(ept_TD %>%
         group_by(trap) %>%
         filter(n() > 3),
       aes(x = year, y = SppRich)) +
  geom_line() +
  geom_smooth(method = "lm") +
  geom_point(size = 2) +
  # geom_hline(yintercept = 0.5,
  #            linetype = "dotted",
  #            color = "grey50",
  #            size = 0.5) +
  facet_wrap(~trap, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(ept_TD$year), max(ept_TD$year), by = 3)) +
  labs(title = "EPT species richness",
       x = "Year",
       y = "Species richness",
       caption = "EPT species richness through time by trap") +
  My_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.25, "cm"))

ggsave(
  "Plots/species_richness_ept.png",
  plot = last_plot(),
  scale = 1,
  width = 10,
  height = 10,
  units = "in",
  dpi = 450,
  device = "png",
  bg = "white"
)

# (Individual-based) rarefied species richness by trap by year
# Diptera
ggplot(diptera_TD %>%
         group_by(trap) %>%
         filter(n() > 3),
       aes(x = year, y = rare.SppRich)) +
  geom_line() +
  geom_smooth(method = "lm") +
  geom_point(size = 2) +
  # geom_hline(yintercept = 0.5,
  #            linetype = "dotted",
  #            color = "grey50",
  #            size = 0.5) +
  facet_wrap(~trap, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(diptera_TD$year), max(diptera_TD$year), by = 3)) +
  labs(title = "(Individual-based) rarefied Diptera species richness",
       x = "Year",
       y = "Rarefied species richness",
       caption = "(Individual-based) rarefied Diptera species richness through time by trap") +
  My_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.25, "cm"))

ggsave(
  "Plots/rarefied_species_richness_diptera.png",
  plot = last_plot(),
  scale = 1,
  width = 10,
  height = 10,
  units = "in",
  dpi = 450,
  device = "png",
  bg = "white"
)

# EPT
ggplot(ept_TD %>%
         group_by(trap) %>%
         filter(n() > 3),
       aes(x = year, y = rare.SppRich)) +
  geom_line() +
  geom_smooth(method = "lm") +
  geom_point(size = 2) +
  facet_wrap(~trap, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(ept_TD$year), max(ept_TD$year), by = 3)) +
  labs(title = "(Individual-based) rarefied EPT species richness",
       x = "Year",
       y = "Rarefied species richness",
       caption = "(Individual-based) rarefied EPT species richness through time by trap") +
  My_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.25, "cm"))

ggsave(
  "Plots/rarefied_species_richness_ept.png",
  plot = last_plot(),
  scale = 1,
  width = 10,
  height = 10,
  units = "in",
  dpi = 450,
  device = "png",
  bg = "white"
)

## ----------- Rarefaction accumulation curves -----------
# diptera_wide_red <- diptera_wide |>
#   select(-D, -F, -Source) |>
#   column_to_rownames(var = "sp_name") |>
#   data.frame() |>
#   print()                                                                             # remove traps with lower quality data: Trap D, F, & Quelle
# diptera_cum.rare = iNEXT(diptera_wide_red, q = 0, datatype = "abundance")             # calculate accumulation curves
# write_rds(diptera_cum.rare, "Outputs/rarefaction_diptera.rds")                        # save output for quicker loading next time
diptera_cum.rare <- readRDS("Outputs/rarefaction_diptera.rds")                           # load saved output
diptera_cum.rare$DataInfo                                                              # check basic outputs
diptera_cum.rare$AsyEst                                                                # Asymptotic diversity estimates

# ept_wide_red <- ept_wide |>
#   select(-O, -F) |>
#   column_to_rownames(var = "sp_name") |>
#   data.frame() |>
#   print()                                                                              # remove traps with lower quality data: Trap D, F, & Quelle
# ept_cum.rare = iNEXT(ept_wide_red, q = 0, datatype = "abundance")                      # calculate accumulation curves
# write_rds(ept_cum.rare, "Outputs/rarefaction_ept.rds")                                 # save output for quicker loading next time
ept_cum.rare <- readRDS("Outputs/rarefaction_ept.rds")                                    # load saved output
ept_cum.rare$DataInfo                                                                   # check basic outputs
ept_cum.rare$AsyEst                                                                     # Asymptotic diversity estimates

### ------------ Visualizations ------------
diptera_p <- ggiNEXT(diptera_cum.rare, type = 1, color.var = "Assemblage")
diptera_p <- diptera_p + ggtitle("Diptera")
diptera_p <- diptera_p + theme(legend.position = "bottom")
diptera_p

diptera_p1 <- ggiNEXT(diptera_cum.rare, facet.var = "Assemblage", grey = T)
diptera_p1 <- diptera_p1 + ggtitle("Diptera")
diptera_p1 <- diptera_p1 + theme(axis.title.x = element_blank(), legend.position = "bottom")
diptera_p1

ept_p <- ggiNEXT(ept_cum.rare, type = 1, color.var = "Assemblage")
ept_p <- ept_p + ggtitle("EPT")
ept_p <- ept_p + theme(axis.title.y = element_blank(), legend.position = "bottom")
ept_p

ept_p1 <- ggiNEXT(ept_cum.rare, facet.var = "Assemblage", grey = T)
ept_p1 <- ept_p1 + ggtitle("EPT")
ept_p1 <- ept_p1 + theme(legend.position = "bottom")
ept_p1

combined_plot <- (diptera_p + ept_p) +
  plot_annotation(
    title = "Comparison of Rarefaction and Extrapolation Curves",
    caption = "iNEXT(dataset, q = 0, datatype = abundance)",
    theme = theme(legend.position='bottom')) +
  plot_layout(ncol = 2, guides = "collect")
combined_plot

ggsave(
  "Plots/Combined_rarefaction_curves.png",
  plot = last_plot(),
  scale = 1,
  width = 15,
  height = 10,
  units = "in",
  dpi = 450,
  device = "png",
  bg = "white"
)

combined_plot1 <- (diptera_p1 + ept_p1) +
  plot_annotation(
    title = "Comparison of Rarefaction and Extrapolation Curves",
    caption = "iNEXT(dataset, q = 0, datatype = abundance)",
    theme = theme(legend.position='bottom')) +
  plot_layout(nrow = 2, guides = "collect")
combined_plot1

ggsave(
  "Plots/Combined_rarefaction_curves_faceted.png",
  plot = last_plot(),
  scale = 1,
  width = 15,
  height = 10,
  units = "in",
  dpi = 450,
  device = "png",
  bg = "white"
)

# ========== CLEAN UP ==========
library(pacman)
# Free unsed memory
gc()
# Clear data
rm(list = ls())  # Removes all objects from environment
# Clear packages
p_unload(all)  # Remove all contributed packages
# Clear plots
graphics.off()  # Clears plots, closes all graphics devices
# Clear console
cat("\014")  # Mimics ctrl+L
# Clear mind :)
