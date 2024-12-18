# Calculate winners and losers for all groups (families?) with mblm and
# then use that metric to compare rates of biodiv accumulation

##########################################
# Package  installation
library(tidyverse)
library(janitor)

##########################################
## Subsetting and preparing Breitenbach datasets for analysis
# Load ept
ept <- as_tibble(read.csv("ept breitenbach_all traps.csv",sep=",")) %>%
  clean_names() %>%
  rename_with(~ gsub("^x", "", .x), .cols = matches("^x[0-9]")) %>%
  mutate(across(-c(ept, breitenbach, x), ~ as.numeric(.x))) %>%
  mutate(across(-c(ept, breitenbach, x), ~ replace_na(.x, 0)))

# Load diptera
diptera <- as_tibble(read.csv("breit_main_diptera1_updated_names.csv",sep=","))  %>%
# Homogenize column names
  clean_names() %>%
# Change column names
  rename_with(~ gsub("^x", "", .x), .cols = matches("^x[0-9]")) %>%
# convert NAs to 0
  mutate(across(-c(diptera, trap, family), ~ as.numeric(.x))) %>%
  mutate(across(-c(diptera, trap, family), ~ replace_na(.x, 0))) %>%
# Homogenize trap names
  mutate(trap = case_when(
    trap %in% c('A/I') ~ 'A',
    trap %in% c('Haus I') ~ 'I', # is this correct?
    trap %in% c('B / II', 'B-II-X', 'B/II', 'Haus B/II') ~ 'B',
    trap %in% c('C / IV', 'C-IV', 'C/IV', 'Haus C', 'Haus C/IV') ~ 'C',
    trap %in% c('Haus III') ~ 'III', # is this correct?
    trap %in% c('E / V', 'E-V', 'E/V', 'Haus E/V') ~ 'E',
    trap %in% c('F / VI', 'F-VI', 'F-VII', 'F/VI', 'Haus F/VI') ~ 'F',
    trap %in% c('G-VII', 'G/ VII', 'G/VII', 'Haus G/VII') ~ 'G',
    trap %in% c('Quelle') ~ 'Quelle', # is this correct?
    TRUE ~ trap)) %>%
# Change column name
  mutate(sp_name = diptera) %>%
# # Create new column with taxa names without sex marker
#   mutate(sp_name = gsub(' k', '', sp_name),
#          sp_name = gsub(' l', '', sp_name),
#          sp_name = gsub(' m', '', sp_name)) %>%
# Reorder columns
  select(diptera, sp_name, family, trap, everything()) %>%
# Sort data
  arrange(sp_name) %>%
# Remove diptera column
  select(-diptera)

# Aggregate data according to species names, traps, and families  - THIS NEEDS WORK DUE TO NAMING ERRORS IN ORIGINAL DATA!!!!
diptera_clean <- aggregate(diptera[, 4:40],
                           list(sp_name = diptera$sp_name,
                                trap = diptera$trap,
                                family = diptera$family),
                           sum) %>%
  # Sort data
  arrange(trap, sp_name)

# Convert diptera dataset from wide to long format
diptera_long <- diptera_clean %>%
  pivot_longer(
    cols = starts_with("19") | starts_with("20"),  # Specify the year column to pivot
    names_to = "year",        # Name column that will hold year values
    values_to = "abundance"   # Name column that will hold the counts
  ) %>%
# create new column with trap_code (trap name + year)
  mutate(trap_code = paste(trap, year, sep = "_")) %>%
# make year variable numeric
  mutate(year = as.numeric(year)) %>%
# Remove rows where count is 0
  filter(abundance != 0) %>%
# Sort data
  arrange(sp_name)

# Reaggregate data according to species names and traps
diptera_agg <- aggregate(diptera_long[, 5], list(sp_name = diptera_long$sp_name,
                                               trap = diptera_long$trap), sum)

# Check if the columns of both the initial and aggregated dataset are the same
identical(sum(diptera_agg$abundance, na.rm = F), sum(diptera_long$abundance, na.rm = F))

# Create species by site (trap) matrix in wide form
library(reshape2)
diptera_wide <- dcast(diptera_agg, sp_name ~ trap, value.var = "abundance", fun = sum)

# Check if the columns of both the initial and aggregated dataset are the same
identical(sum(diptera_agg$abundance, na.rm = F), sum(diptera_long$abundance, na.rm = F), sum(diptera_wide[, c(2:11)], na.rm = F))

# Isolate trap abundances
rownames(diptera_wide) <- diptera_wide$sp_name # change row names to sp_names
diptera_wide <- as.data.frame(diptera_wide[,2:11])

# calculate indices
library(SpadeR)
chaoA_dipt      <- ChaoSpecies(diptera_wide$A ,"abundance", k = 2, conf = 0.95)      # Diptera trap A
chaoB_dipt      <- ChaoSpecies(diptera_wide$B ,"abundance", k = 2, conf = 0.95)      # Diptera trap B
chaoC_dipt      <- ChaoSpecies(diptera_wide$C ,"abundance", k = 2, conf = 0.95)      # Diptera trap C
chaoE_dipt      <- ChaoSpecies(diptera_wide$E ,"abundance", k = 2, conf = 0.95)      # Diptera trap E
chaoG_dipt      <- ChaoSpecies(diptera_wide$G ,"abundance", k = 2, conf = 0.95)      # Diptera trap G
chaoI_dipt      <- ChaoSpecies(diptera_wide$I ,"abundance", k = 2, conf = 0.95)      # Diptera trap I
chaoIII_dipt    <- ChaoSpecies(diptera_wide$III ,"abundance", k = 2, conf = 0.95)    # Diptera trap III

# Univariate diversity metrics
library(vegan)
library(codyn)
library(mobr)

# Recreate long dataset
diptera_long <- diptera_clean %>%
  pivot_longer(
    cols = starts_with("19") | starts_with("20"),  # Specify the year column to pivot
    names_to = "year",        # Name column that will hold year values
    values_to = "abundance"   # Name column that will hold the counts
  ) %>%
  # create new column with trap_code (trap name + year)
  mutate(trap_code = paste(trap, year, sep = "_")) %>%
  # make year variable numeric
  mutate(year = as.numeric(year)) %>%
  # Remove rows where count is 0
  filter(abundance != 0) %>%
  # Sort data
  arrange(sp_name)

# Individual-based rarefication of species richness
diptera_long$ro.ab <- round(diptera_long$abundance, digits = 0) #rarefication only works on integers

TD <- NULL
for(i in unique(diptera_long$trap)){
  sub <- diptera_long[diptera_long$trap == i, ]
  sub.m <- dcast(sub, trap_code ~ sp_name, sum, value.var = "abundance") 	# matrix form
  sub.ta <- subset(sub.m[,c(2:length(sub.m))])          # subset matrix to remove row names
  SppRich <- specnumber(sub.ta)  							          # taxonomic richness
  Simp <- diversity(sub.ta, index = "simpson") 					# Simpson's taxonomic diversity
  Shan <- diversity(sub.ta, index = "shannon")					# Shannon's taxonomic diversity
  EvenJ <- Shan/log(SppRich) 							    	        # Pielou's evenness (J)
  E10 <- Shan/SppRich 								  	              # Shannon's evenness (E10)
  Abund <- rowSums (sub.ta) 								            # Total abundance
  S_PIE <- calc_PIE(sub.ta, ENS = TRUE)						      # effective number of common species
  DATA1_Turnover <- turnover(sub, time.var = "year", species.var = "sp_name", abundance.var = "abundance", replicate.var = NA, metric = "total")
  Turnover <- c("NA", DATA1_Turnover$total) 					  # Turnover per yr And first yr is "NA"
  DATA1_Turnover_app <- turnover(sub, time.var = "year", species.var = "sp_name", abundance.var = "abundance", replicate.var = NA, metric = "appearance")
  Turnover_app <- c("NA", DATA1_Turnover_app$appearance) 					  # Turnover per yr And first yr is "NA"
  DATA1_Turnover_disapp <- turnover(sub, time.var = "year", species.var = "sp_name", abundance.var = "abundance", replicate.var = NA, metric = "disappearance")
  Turnover_disapp <- c("NA", DATA1_Turnover_disapp$disappearance) 					  # Turnover per yr And first yr is "NA"
  sub.m_r <- dcast(sub, trap_code ~ sp_name, sum, value.var = "ro.ab")    # matrix form for rarefaction with rounded richness
  sub.ta_r <- subset(sub.m_r[,c(2:length(sub.m_r))])                  	# subset matrix to remove row names for rounded sppRich
  rare.SppRich <- if (min(rowSums(sub.ta_r)) > 10) {
    rarefy(sub.ta_r, sample = min(rowSums(sub.ta_r)))			# rarefy based on min abundance
  } else {rarefy(sub.ta_r, sample = 10)} 					        # rarefy based on abund = 10 if min is less
  TD.i <- data.frame(sub.m$trap_code, SppRich, Simp, Shan, EvenJ, E10, Abund, S_PIE, Turnover, Turnover_app, Turnover_disapp, rare.SppRich)
  TD <- rbind(TD, TD.i) ; rm(TD.i, sub.m, sub.ta, sub, SppRich, Simp, Shan,
                             EvenJ, E10, Abund, S_PIE, DATA1_Turnover, Turnover, DATA1_Turnover_app, Turnover_app, DATA1_Turnover_disapp, Turnover_disapp, sub.m_r, sub.ta_r, rare.SppRich)
} ; rm(i)

colnames(TD)[1] <- "trap_code"

# Add "trap" and "year" columns by splitting the "trap_code"
TD <- TD %>%
  mutate(trap = gsub("_.*", "", trap_code),
         year = gsub(".*_", "", trap_code)) %>%
  # mutate(EvenJ = ifelse(is.nan(EvenJ), 0, EvenJ)) %>%
  mutate(across(-c(trap_code, trap), as.numeric)) %>%
  # mutate(across(starts_with("Turnover"), ~ replace_na(., 0))) %>%
  select(trap_code, trap, year, everything())

# Lets do some plotting
# Define 'My_theme' for plotting
My_theme <- theme(panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, linewidth = 1.25),
                  strip.background = element_rect(fill = "white",
                                                  color = "white", linewidth = 1.25),
                  legend.position = "bottom",
                  text = element_text(size = 16,
                                      family = "gillsans"))

# Turnover by year
ggplot(TD %>%
         group_by(trap) %>%
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
  scale_x_continuous(breaks = seq(min(TD$year), max(TD$year), by = 3)) +
  scale_y_continuous(limits = c(0.25, max(TD$Turnover, na.rm = TRUE))) +
  labs(title = "Turnover",
       x = "Year",
       y = "Turnover",
       caption = "Species turnover over time by trap") +
  My_theme +
  My_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.25, "cm"))

# Turnover appearences by year
ggplot(TD %>%
         group_by(trap) %>%
         filter(n() > 3),
       aes(x = year, y = Turnover_app)) +
  geom_line() +
  geom_smooth(method = "lm") +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.5,
             linetype = "dotted",
             color = "grey50",
             size = 0.5) +
  facet_wrap(~trap) +
  scale_x_continuous(breaks = seq(min(TD$year), max(TD$year), by = 3)) +
  scale_y_continuous(limits = c(0, max(TD$Turnover, na.rm = TRUE))) +
  labs(title = "Turnover: appearences",
       x = "Year",
       y = "Turnover appearences",
       caption = "Species appearences over time by trap") +
  My_theme +
  My_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.25, "cm"))

# Turnover disappearences by year
ggplot(TD %>%
         group_by(trap) %>%
         filter(n() > 3),
       aes(x = year, y = Turnover_disapp)) +
  geom_line() +
  geom_smooth(method = "lm") +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.5,
             linetype = "dotted",
             color = "grey50",
             size = 0.5) +
  facet_wrap(~trap) +
  scale_x_continuous(breaks = seq(min(TD$year), max(TD$year), by = 3),
                     minor_breaks = seq(min(TD$year), max(TD$year), by = 1)) +
  scale_y_continuous(limits = c(0, max(TD$Turnover, na.rm = TRUE))) +
  labs(title = "Turnover: disappearences",
       x = "Year",
       y = "Turnover Disappearences",
       caption = "Species dissappearences through time by traps") +
  My_theme +
  My_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.25, "cm"))

# Abundance by trap by year
ggplot(TD %>%
         group_by(trap) %>%
         filter(n() > 3),
       aes(x = year, y = Abund)) +
  geom_line() +
  geom_smooth(method = "lm") +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.5,
             linetype = "dotted",
             color = "grey50",
             size = 0.5) +
  facet_wrap(~trap, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(TD$year), max(TD$year), by = 3)) +
  labs(title = "Abundance",
       x = "Year",
       y = "Abundance",
       caption = "Abundance through time by trap") +
  My_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.25, "cm"))

# Log abundance by trap by year
ggplot(TD %>%
         group_by(trap) %>%
         filter(n() > 3),
       aes(x = year, y = log(Abund + 1))) +
  geom_line() +
  geom_smooth(method = "lm") +
  geom_point(size = 2) +
  # geom_hline(yintercept = 0.5,
  #            linetype = "dotted",
  #            color = "grey50",
  #            size = 0.5) +
  facet_wrap(~trap, scales = "free_y") +
  scale_y_continuous(breaks = seq(min(0), max(10), by = 2)) +
  scale_x_continuous(breaks = seq(min(TD$year), max(TD$year), by = 3)) +
  labs(title = "Log abundance",
       x = "Year",
       y = "log(abundance + 1)",
       caption = "Log abundance through time by trap") +
  My_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.25, "cm"))

# Species richness by trap by year
ggplot(TD %>%
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
  scale_x_continuous(breaks = seq(min(TD$year), max(TD$year), by = 3)) +
  labs(title = "Species richness",
       x = "Year",
       y = "Species richness",
       caption = "Species richness through time by trap") +
  My_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.25, "cm"))

# Rarefied species richness by trap by year
ggplot(TD %>%
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
  scale_x_continuous(breaks = seq(min(TD$year), max(TD$year), by = 3)) +
  labs(title = "Rarefied species richness",
       x = "Year",
       y = "Rarefied species richness",
       caption = "Rarefied species richness through time by trap") +
  My_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.25, "cm"))

# Calculate rarefaction as cumulative curve (i.e. Cumulative rarefaction)
library(iNEXT)
# remove traps with lower quality data: Trap D, F, & Quelle
diptera_wide_red <- diptera_wide %>%
  select(-D, -F, -Quelle)
# calculate curves
diptera_cum.rare = iNEXT(diptera_wide_red, q = 0, datatype = "abundance")
# check outputs
diptera_cum.rare$DataInfo  # show basic data information.
diptera_cum.rare$AsyEst    # show asymptotic diversity estimates.
# plotting accumulation curves
# set plotting theme
My_theme <- theme(panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, linewidth = 1.25),
                  strip.background = element_rect(fill = "white",
                                                  color = "white", linewidth = 1.25),
                  legend.position = "bottom",
                  text = element_text(size = 16,
                                      family = "gillsans"),
                  axis.text.x = element_text(angle = 45, hjust = 1))
# plot all accumulation curves in the same graph coloured by assemblages
p <- ggiNEXT(diptera_cum.rare, type = 1, color.var = "Assemblage")
p <- p + My_theme
p
# plot all accumulation curves in seperate graphs (facet_wrap)
p1 <- ggiNEXT(diptera_cum.rare, facet.var = "Assemblage", grey = T)
p1 <- p1 + My_theme
p1


###############################################################
# CLEAN UP
library(pacman)
# Clear data
rm(list = ls())  # Removes all objects from environment
# Clear packages
p_unload(all)  # Remove all contributed packages
# Clear plots
graphics.off()  # Clears plots, closes all graphics devices
# Clear console
cat("\014")  # Mimics ctrl+L
# Clear mind :)
