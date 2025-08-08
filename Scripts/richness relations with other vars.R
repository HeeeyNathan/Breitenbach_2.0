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
require(gridExtra)
library(data.table)
library(entropart)
library(future.apply)
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)

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

#======== Separate Trap and Year=====




# Clean up: separate Trap and Year from "Trap_Year" if needed
diptera_long <- diptera_long %>%
  separate(trap_code, into = c("TrapID", "YearID"), sep = "_") %>%
  mutate(YearID = as.integer(YearID))  # Ensure YearID is numeric

#==================================================================================================
#Diptera
# ----
# Step 1: Create species × site × year matrix
# ----

# Aggregate counts (in case there are repeated records)
agg_dipt <- diptera_long %>%
  group_by(TrapID, YearID, sp_name) %>%
  summarise(Count = sum(abundance), .groups = "drop")
keep_levels <- c("A", "C", "G", "III", "B", "E", "I")

# Filter the dataframe
agg_dipt <- agg_dipt %>%
  filter(TrapID %in% keep_levels) %>%
  mutate(TrapID = factor(TrapID, levels = keep_levels))

# Pivot to species abundance matrix per trap and year
comm_matrix <- agg_dipt %>%
  pivot_wider(names_from = sp_name, values_from = Count, values_fill = 0)

# Store trap-year info separately
trap_year_meta <- comm_matrix %>%
  select(TrapID, YearID)

# Drop metadata columns to retain numeric matrix for diversity analysis
species_matrix <- comm_matrix %>%
  select(-TrapID, -YearID)

# ----
# Step 2: Calculate metrics per trap and year
# ----

# Compute N (total individuals), S (richness), PIE
metrics <- species_matrix %>%
  mutate(N = rowSums(.),
         S = specnumber(.),
         PIE = diversity(., index = "simpson")) %>%
  bind_cols(trap_year_meta, .)

# ----


# Step 3 (updated): Compute deltas for each consecutive year per trap

# Function to calculate deltas between consecutive years
compute_all_deltas <- function(df_sub) {
  df_sub <- df_sub %>% arrange(YearID)

  # If fewer than 2 years, skip
  if (nrow(df_sub) < 2) return(NULL)

  # Store results
  result_list <- list()

  for (i in 2:nrow(df_sub)) {
    year1 <- df_sub[i - 1, ]
    year2 <- df_sub[i, ]

    delta_N <- year2$N - year1$N
    delta_S <- year2$S - year1$S
    delta_PIE <- year2$PIE - year1$PIE

    # Extract species abundance data for rarefaction
    species_1 <- year1 %>% select(-TrapID, -YearID, -N, -S, -PIE)
    species_2 <- year2 %>% select(-TrapID, -YearID, -N, -S, -PIE)

    rarefied <- rarefy(rbind(species_1, species_2), sample = min(year1$N, year2$N))

    result_list[[i - 1]] <- tibble(
      TrapID = year1$TrapID,
      Year1 = year1$YearID,
      Year2 = year2$YearID,
      Delta_N = delta_N,
      Delta_S = delta_S,
      S_rare_1 = rarefied[1],
      S_rare_2 = rarefied[2],
      Delta_PIE = delta_PIE
    )
  }

  bind_rows(result_list)
}

# Apply the delta function across traps
delta_results_all_years <- metrics %>%
  group_by(TrapID) %>%
  group_split() %>%
  lapply(compute_all_deltas) %>%
  bind_rows()

# Final output
print(delta_results_all_years)

#Correlations


# Create Delta_S_rare
delta_results_all_years <- delta_results_all_years %>%
  mutate(Delta_S_rare = S_rare_2 - S_rare_1)


# Correlation: Delta_N vs Delta_S
cor_N_S <- cor.test(
  delta_results_all_years$Delta_N,
  delta_results_all_years$Delta_S,
  method = "pearson"
)

# Correlation: Delta_S_rare vs Delta_S
cor_Srare_S <- cor.test(
  delta_results_all_years$Delta_S_rare,
  delta_results_all_years$Delta_S,
  method = "pearson"
)

# Print results
print(cor_N_S)
print(cor_Srare_S)


# Correlations per TrapID
correlations_by_trap <- delta_results_all_years %>%
  group_by(TrapID) %>%
  summarise(
    cor_N_S = cor(Delta_N, Delta_S, method = "pearson", use = "complete.obs"),
    cor_Srare_S = cor(Delta_S_rare, Delta_S, method = "pearson", use = "complete.obs"),
    .groups = "drop"
  )

# View table
print(correlations_by_trap)


#Part 4 plots
# Load required libraries

# Calculate Delta_S_rare (rarefied richness change)
delta_results_all_years <- delta_results_all_years %>%
  mutate(Delta_S_rare = S_rare_2 - S_rare_1)

# Plot 1: ΔS vs ΔN with 95% CI ellipse
p1 <- ggplot(delta_results_all_years, aes(x = Delta_N, y = Delta_S)) +
  geom_point(alpha = 0.6, color = "#E377C2") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  stat_ellipse(type = "norm", level = 0.95, color = "gray40") +
  theme_minimal(base_size = 13) +
  labs(title = "ΔS vs ΔN", x = "ΔN (Change in Abundance)", y = "ΔS (Change in Richness)")
p1 <- p1 +
  annotate("text", x = 0, y = 150, label = "r = 0.54",
           size = 5, fontface = "italic", color = "black")

# Plot 2: ΔS vs ΔS_rare with 95% CI ellipse
p2 <- ggplot(delta_results_all_years, aes(x = Delta_S_rare, y = Delta_S)) +
  geom_point(alpha = 0.6, color = "#FFCE54") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  stat_ellipse(type = "norm", level = 0.95, color = "gray40") +
  theme_minimal(base_size = 13) +
  labs(title = "ΔS vs ΔS_rare", x = "ΔS_rare (Rarefied Richness Change)", y = "ΔS (Observed Richness Change)")

p2<- p2 +
  annotate("text", x = 0, y = 150, label = "r = 0.82",
           size = 5, fontface = "italic", color = "black")
# === DIPTERA ===

# p1: ΔS vs ΔN with facet
p1a <- ggplot(delta_results_all_years, aes(x = Delta_N, y = Delta_S)) +
  geom_point(alpha = 0.6, color = "#E377C2") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  stat_ellipse(level = 0.95, color = "gray40") +
  facet_wrap(~ TrapID) +
  theme_minimal(base_size = 14) +
  labs(title = "ΔS vs ΔN", x = "ΔN (Change in Abundance)", y = "ΔS (Change in Richness)")

# p2: ΔS vs ΔS_rare with facet
p2a <- ggplot(delta_results_all_years, aes(x = Delta_S_rare, y = Delta_S)) +
  geom_point(alpha = 0.6, color = "#FFCE54") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  stat_ellipse(level = 0.95, color = "gray40") +
  facet_wrap(~ TrapID) +
  theme_minimal(base_size = 14) +
  labs(title = "ΔS vs ΔS_rare", x = "ΔS_rare (Rarefied Richness Change)", y = "ΔS (Observed Richness Change)")

#==============================================
#EPT


# ----
# Step 1: Create species × site × year matrix
# ----

# Clean up: separate Trap and Year from "Trap_Year" if needed
ept_long <- ept_long %>%
  separate(trap_code, into = c("TrapID", "YearID"), sep = "_") %>%
  mutate(YearID = as.integer(YearID))  # Ensure YearID is numeric



# Aggregate counts (in case there are repeated records)
agg_ept <- ept_long %>%
  group_by(TrapID, YearID, sp_name) %>%
  summarise(Count = sum(abundance), .groups = "drop")

# Pivot to species abundance matrix per trap and year
comm_matrix <- agg_ept %>%
  pivot_wider(names_from = sp_name, values_from = Count, values_fill = 0)

# Store trap-year info separately
trap_year_meta <- comm_matrix %>%
  select(TrapID, YearID)

# Drop metadata columns to retain numeric matrix for diversity analysis
species_matrix <- comm_matrix %>%
  select(-TrapID, -YearID)

# ----
# Step 2: Calculate metrics per trap and year
# ----

# Compute N (total individuals), S (richness), PIE
metrics <- species_matrix %>%
  mutate(N = rowSums(.),
         S = specnumber(.),
         PIE = diversity(., index = "simpson")) %>%
  bind_cols(trap_year_meta, .)

# ----


# Step 3 (updated): Compute deltas for each consecutive year per trap

# Function to calculate deltas between consecutive years
compute_all_deltas <- function(df_sub) {
  df_sub <- df_sub %>% arrange(YearID)

  # If fewer than 2 years, skip
  if (nrow(df_sub) < 2) return(NULL)

  # Store results
  result_list <- list()

  for (i in 2:nrow(df_sub)) {
    year1 <- df_sub[i - 1, ]
    year2 <- df_sub[i, ]

    delta_N <- year2$N - year1$N
    delta_S <- year2$S - year1$S
    delta_PIE <- year2$PIE - year1$PIE

    # Extract species abundance data for rarefaction
    species_1 <- year1 %>% select(-TrapID, -YearID, -N, -S, -PIE)
    species_2 <- year2 %>% select(-TrapID, -YearID, -N, -S, -PIE)

    rarefied <- rarefy(rbind(species_1, species_2), sample = min(year1$N, year2$N))

    result_list[[i - 1]] <- tibble(
      TrapID = year1$TrapID,
      Year1 = year1$YearID,
      Year2 = year2$YearID,
      Delta_N = delta_N,
      Delta_S = delta_S,
      S_rare_1 = rarefied[1],
      S_rare_2 = rarefied[2],
      Delta_PIE = delta_PIE
    )
  }

  bind_rows(result_list)
}

# Apply the delta function across traps
delta_results_all_years <- metrics %>%
  group_by(TrapID) %>%
  group_split() %>%
  lapply(compute_all_deltas) %>%
  bind_rows()


keep_levels <- c("A", "C", "G", "III", "B", "E", "I")

# Filter the dataframe
delta_results_all_years<- delta_results_all_years %>%
  filter(TrapID %in% keep_levels) %>%
  mutate(TrapID = factor(TrapID, levels = keep_levels))

# Final output
print(delta_results_all_years)



# Create Delta_S_rare
delta_results_all_years <- delta_results_all_years %>%
  mutate(Delta_S_rare = S_rare_2 - S_rare_1)

# Correlation: Delta_N vs Delta_S
cor_N_S <- cor.test(
  delta_results_all_years$Delta_N,
  delta_results_all_years$Delta_S,
  method = "pearson"
)

# Correlation: Delta_S_rare vs Delta_S
cor_Srare_S <- cor.test(
  delta_results_all_years$Delta_S_rare,
  delta_results_all_years$Delta_S,
  method = "pearson"
)

# Print results
print(cor_N_S)
print(cor_Srare_S)


# Correlations per TrapID
correlations_by_trap <- delta_results_all_years %>%
  group_by(TrapID) %>%
  summarise(
    cor_N_S = cor(Delta_N, Delta_S, method = "pearson", use = "complete.obs"),
    cor_Srare_S = cor(Delta_S_rare, Delta_S, method = "pearson", use = "complete.obs"),
    .groups = "drop"
  )

# View table
print(correlations_by_trap)



#4 Plot results

# Calculate Delta_S_rare (rarefied richness change)
delta_results_all_years <- delta_results_all_years %>%
  mutate(Delta_S_rare = S_rare_2 - S_rare_1)

# Plot 1: ΔS vs ΔN with 95% CI ellipse
p3 <- ggplot(delta_results_all_years, aes(x = Delta_N, y = Delta_S)) +
  geom_point(alpha = 0.6, color = "#E377C2") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  stat_ellipse(type = "norm", level = 0.95, color = "gray40") +
  theme_minimal(base_size = 13) +
  labs(title = "ΔS vs ΔN", x = "ΔN (Change in Abundance)", y = "ΔS (Change in Richness)")
p3 <- p3 +
  annotate("text", x = 0, y = 15, label = "r = 0.25",
           size = 5, fontface = "italic", color = "black")

# Plot 2: ΔS vs ΔS_rare with 95% CI ellipse
p4 <- ggplot(delta_results_all_years, aes(x = Delta_S_rare, y = Delta_S)) +
  geom_point(alpha = 0.6, color = "#FFCE54") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  stat_ellipse(type = "norm", level = 0.95, color = "gray40") +
  theme_minimal(base_size = 13) +
  labs(title = "ΔS vs ΔS_rare", x = "ΔS_rare (Rarefied Richness Change)", y = "ΔS (Observed Richness Change)")
p4 <- p4 +
  annotate("text", x = 0, y = 15, label = "r = 0.84",
           size = 5, fontface = "italic", color = "black")


# p1: ΔS vs ΔN with facet
p3a <- ggplot(delta_results_all_years, aes(x = Delta_N, y = Delta_S)) +
  geom_point(alpha = 0.6, color = "#E377C2") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  stat_ellipse(level = 0.95, color = "gray40") +
  facet_wrap(~ TrapID) +
  theme_minimal(base_size = 14) +
  labs(title = "ΔS vs ΔN", x = "ΔN (Change in Abundance)", y = "ΔS (Change in Richness)")

# p2: ΔS vs ΔS_rare with facet
p4a <- ggplot(delta_results_all_years, aes(x = Delta_S_rare, y = Delta_S)) +
  geom_point(alpha = 0.6, color = "#FFCE54") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  stat_ellipse(level = 0.95, color = "gray40") +
  facet_wrap(~ TrapID) +
  theme_minimal(base_size = 14) +
  labs(title = "ΔS vs ΔS_rare", x = "ΔS_rare (Rarefied Richness Change)", y = "ΔS (Observed Richness Change)")

# Combine both plots in a grid

# Calculate Delta_S_rare
delta_results_all_years <- delta_results_all_years %>%
  mutate(Delta_S_rare = S_rare_2 - S_rare_1)

# Labels for A-D
label_A <- textGrob("A", x = unit(0, "npc"), y = unit(1, "npc"),
                    just = c("left", "top"),
                    gp = gpar(fontsize = 16, fontface = "bold"))
label_B <- textGrob("B", x = unit(0, "npc"), y = unit(1, "npc"),
                    just = c("left", "top"),
                    gp = gpar(fontsize = 16, fontface = "bold"))
label_C <- textGrob("C", x = unit(0, "npc"), y = unit(1, "npc"),
                    just = c("left", "top"),
                    gp = gpar(fontsize = 16, fontface = "bold"))
label_D <- textGrob("D", x = unit(0, "npc"), y = unit(1, "npc"),
                    just = c("left", "top"),
                    gp = gpar(fontsize = 16, fontface = "bold"))

# Create plots with labels

# Section titles
top_label <- textGrob("Diptera", gp = gpar(fontsize = 18, fontface = "bold"))
bottom_label <- textGrob("EPT", gp = gpar(fontsize = 18, fontface = "bold"))

# Arrange top row and bottom row


# Arrange top row and bottom row
top_row <- arrangeGrob(p1, p2, ncol = 2)
bottom_row <- arrangeGrob(p3, p4, ncol = 2)

# Final grid with section titles
grid.arrange(
  top_label,
  top_row,
  bottom_label,
  bottom_row,
  ncol = 1,
  heights = c(0.5, 10, 0.5, 10)
)
