# ========== LOAD PACKAGES ==========
library(tidyverse)
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


#===============================================================================================

# R scripts for computing diversity (Hill numbers) profile using individual-based abundance data or sampling-unit-based incidence data.
# In all functions, param x is a vector of species sample frequencies (for abundance data) or incidence-based sample frequencies (for incidence data).
# For incidence data, the first entry of x must be the number of sampling units.
# In all functions, param q is the diversity order; the suggested range for q is [0, 3].
# If you use the scripts for publishing papers, please cite Chao and Jost 2015 MEE paper (Appendix S8).

#-----------------------------------------------
# Diversity profile estimator (abundance data)
#-----------------------------------------------
#' Chao_Hill_abu(x, q) is a function of obtaining estimators of Hill numbers of order q based on abundance data.
#' @param x a vector of species sample frequencies.
#' @param q a numeric or a vector of diversity order; The suggested range for q is [0, 3].
#' @return a numerical vector of diversity.

Chao_Hill_abu = function(x,q){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  r <- 0:(n-1)
  Sub <- function(q){
    if(q==0){
      sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }
    else if(q==1){
      A <- sum(x/n*(digamma(n)-digamma(x)))
      B <- ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(-log(p1)-sum(sapply(1:(n-1), function(r)(1-p1)^r/r))))
      exp(A+B)
    }else if(abs(q-round(q))==0){
      A <- sum(exp(lchoose(x,q)-lchoose(n,q)))
      ifelse(A==0,NA,A^(1/(1-q)))
    }else {
      sort.data = sort(unique(x))
      tab = table(x)
      term = sapply(sort.data,function(z){
        k=0:(n-z)
        sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
      })
      A = sum(tab*term)
      B = ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(p1^(q-1)-sum(choose(q-1,r)*(p1-1)^r)))
      (A+B)^(1/(1-q))
    }
  }
  sapply(q, Sub)
}

#-----------------------------------------------
# Diversity profile estimator (incidence data)
#-----------------------------------------------
#' Chao_Hill_inc(x, q) is a function of obtaining estimators of Hill numbers of order q based on incidence data.
#' @param x a vector of species incidence-based sample frequencies. The first entry of x must be the number of sampling units.
#' @param q a numeric or a vector of diversity order.
#' @return a numerical vector.

Chao_Hill_inc = function(x,q){
  n = x[1]
  x = x[-1];x = x[x>0]
  U = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  r <- 0:(n-1)
  Sub <- function(q){
    if(q==0){
      sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }
    else if(q==1){
      A <- sum(x/U*(digamma(n)-digamma(x)))
      B <- ifelse(f1==0|p1==1,0,f1/U*(1-p1)^(1-n)*(-log(p1)-sum(sapply(1:(n-1), function(r)(1-p1)^r/r))))
      exp(A+B)*U/n
    }else if(abs(q-round(q))==0){
      A <- sum(exp(lchoose(x,q)-lchoose(n,q)))
      ifelse(A==0,NA,((n/U)^q*A)^(1/(1-q)))
    }else {
      sort.data = sort(unique(x))
      tab = table(x)
      term = sapply(sort.data,function(z){
        k=0:(n-z)
        sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
      })
      A = sum(tab*term)
      B = ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(p1^(q-1)-sum(choose(q-1,r)*(p1-1)^r)))
      ((n/U)^q*(A+B))^(1/(1-q))
    }
  }
  sapply(q, Sub)
}

#' Chao_Hill(x, q,datatype) combines Chao_Hill_abu and Chao_Hill_inc given a specified datatype (either abundance data or incidence data).
#' @param x a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param q a numeric or a vector of diversity order.
#' @param datatype a character of data type,"abundance" or "incidence".
#' @return a numerical vector.

Chao_Hill = function(x,q,datatype = c("abundance","incidence")){
  datatype = match.arg(datatype,c("abundance","incidence"))
  if(datatype == "abundance"){
    est = Chao_Hill_abu(x,q)
  }else{
    est = Chao_Hill_inc(x,q)
  }
  return(est)
}

#-----------------------
# The empirical profile
#-----------------------
#' Hill(x, q, datatype) is a function of obtaining the empirical Hill numbers of order q based on abundance data or incidence data.
#' @param x a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param q a numeric or a vector of diversity order.
#' @param datatype a character of data type,"abundance" or "incidence".
#' @return a numerical vector.

Hill <- function(x,q,datatype = c("abundance","incidence")){
  if(datatype=="incidence"){x = x[-1]}
  p <- x[x>0]/sum(x)
  Sub <- function(q){
    if(q==0) sum(p>0)
    else if(q==1) exp(-sum(p*log(p)))
    else exp(1/(1-q)*log(sum(p^q)))
  }
  sapply(q, Sub)
}

#-----------------------------------------
# The bootstrap method for obtaining s.e.
#-----------------------------------------
#' Bt_prob_abu(x) is a function of estimating the species probabilities in the bootstrap assemblage based on abundance data.
#' @param x a vector of species sample frequencies.
#' @return a numeric vector.

Bt_prob_abu = function(x){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  C = 1 - f1/n*ifelse(f2>0,(n-1)*f1/((n-1)*f1+2*f2),ifelse(f1>0,(n-1)*(f1-1)/((n-1)*(f1-1)+2),0))
  W = (1-C)/sum(x/n*(1-x/n)^n)

  p.new = x/n*(1-W*(1-x/n)^n)
  f0 = ceiling(ifelse(f2>0,(n-1)/n*f1^2/(2*f2),(n-1)/n*f1*(f1-1)/2))
  p0 = (1-C)/f0
  p.new=c(p.new,rep(p0,f0))
  return(p.new)
}

#' Bt_prob_inc(x) is a function of estimating the species incidence probabilities in the bootstrap assemblage based on incidence data.
#' @param x a vector of incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @return a numeric vector.

Bt_prob_inc = function(x){
  n = x[1]
  x = x[-1]
  U = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  A = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  C=1-f1/U*(1-A)
  W=U/n*(1-C)/sum(x/n*(1-x/n)^n)

  p.new=x/n*(1-W*(1-x/n)^n)
  f0 = ceiling(ifelse(f2>0,(n-1)/n*f1^2/(2*f2),(n-1)/n*f1*(f1-1)/2))
  p0=U/n*(1-C)/f0
  p.new=c(p.new,rep(p0,f0))
  return(p.new)
}

#' Bt_prob(x,datatype) combines the two functions Bt_prob_abu and Bt_prob_inc for a specified datatype.
#' @param x a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param datatype a character of data type,"abundance" or "incidence".
#' @return a numeric vector.

Bt_prob = function(x,datatype = c("abundance","incidence")){
  datatype = match.arg(datatype,c("abundance","incidence"))
  if(datatype == "abundance"){
    prob = Bt_prob_abu(x)
  }else{
    prob = Bt_prob_inc(x)
  }
  return(prob)
}

#' Bootstrap.CI(x,q,B,datatype,conf) is a function of calculating the bootsrapping standard error based on abundance data or incidence data.
#' @param x a vector of species sample frequencies (for abundance data) or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param q a numeric or a vector of diversity order.
#' @param B an integer to specify the number of replications in the bootstrap procedure, B = 1000 is suggested for constructing confidence intervals;
#'  To save running time, use a smaller value (e.g. B = 200)..
#' @param datatype a character of data type,"abundance" or "incidence".
#' @param conf a confidence coefficient between 0 and 1.
#' @return a list, consisting of 3 matrices including respectively the difference between the average and lower confidence bound of the B bootstrap estimates,
#'  the difference between the upper confidence bound and the average of the B bootstrap estimates, and the bootstrap standard error of the diversity estimate.
#'  In each matrix, the first row gives the results for the empirical diversity, and the second row gives the results for the proposed diversity estimates.
#'  Columns give the results for different orders of q.

Bootstrap.CI = function(x,q,B = 1000, datatype = c("abundance","incidence"),conf = 0.95){
  datatype = match.arg(datatype,c("abundance","incidence"))
  p.new = Bt_prob(x,datatype)
  n = ifelse(datatype=="abundance",sum(x),x[1])
  # set.seed(456)
  if(datatype=="abundance"){
    data.bt = rmultinom(B,n,p.new)
  }else{
    data.bt = rbinom(length(p.new)*B,n,p.new)
    data.bt = matrix(data.bt,ncol=B)
    data.bt = rbind(rep(n,B),data.bt)
  }

  mle = apply(data.bt,2,function(x)Hill(x,q,datatype))
  pro = apply(data.bt,2,function(x)Chao_Hill(x,q,datatype))

  mle.mean = rowMeans(mle)
  pro.mean = rowMeans(pro)

  LCI.mle =  -apply(mle,1,function(x)quantile(x,probs = (1-conf)/2)) + mle.mean
  UCI.mle = apply(mle,1,function(x)quantile(x,probs = 1-(1-conf)/2)) - mle.mean

  LCI.pro =  -apply(pro,1,function(x)quantile(x,probs = (1-conf)/2)) + pro.mean
  UCI.pro = apply(pro,1,function(x)quantile(x,probs = 1-(1-conf)/2)) - pro.mean

  LCI = rbind(LCI.mle,LCI.pro)
  UCI = rbind(UCI.mle,UCI.pro)

  sd.mle = apply(mle,1,sd)
  sd.pro = apply(pro,1,function(x)sd(x,na.rm = T))
  se = rbind(sd.mle,sd.pro)

  return(list(LCI=LCI,UCI=UCI,se=se))

}

#------------------------
# Main function ChaoHill
#------------------------
#' ChaoHill(dat, datatype, from, to, interval, B, conf) is the function of calculating the empirical and the proposed diversity profile,
#' their bootsrap standard errors and confidance intervals.
#' @param dat a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param datatype a character of data type,"abundance" or "incidence".
#' @param from a numeric number of diversity order q (the start order of profile).
#' @param to a numeric number of diversity order q (the end order of profile).
#' @param interval a numeric number to specify each increment of q from the start to end order.
#' @param B an integer to specify the number of bootstrap replications, B = 1000 is suggested for constructing confidence intervals;
#'  To save running time, use a smaller value (e.g. B = 200).
#' @param conf a confidence coefficient between 0 and 1.
#' @return a list, consisting of 4 matrices including respectively diversity estimates, bootstrap standard errors, lower confidence bounds, and upper confidence bounds.
#'  In each matrix, the first row gives the results for the empirical diversity, and the second row gives the results for the proposed diversity estimates.
#'  Columns give the results for different orders of q.

ChaoHill <- function(dat, datatype=c("abundance", "incidence"), from=0, to=3, interval=0.1, B=1000, conf=0.95){
  datatype = match.arg(datatype,c("abundance","incidence"))
  # for real data estimation

  if (is.matrix(dat) == T || is.data.frame(dat) == T){
    if (ncol(dat) != 1 & nrow(dat) != 1)
      stop("Error: The data format is wrong.")
    if (ncol(dat) == 1){
      dat <- dat[, 1]
    } else {
      dat <- dat[1, ]
    }
  }
  dat <- as.numeric(dat)
  q <- seq(from, to, by=interval)

  #-------------
  #Estimation
  #-------------
  MLE=Hill(dat,q,datatype)

  qD_pro=Chao_Hill(dat,q,datatype)

  CI_bound = Bootstrap.CI(dat,q,B,datatype,conf)
  se = CI_bound$se
  #-------------------
  #Confidence interval
  #-------------------
  tab.est=data.frame(rbind(MLE,qD_pro))

  LCI <- tab.est - CI_bound$LCI
  UCI <- tab.est + CI_bound$UCI

  colnames(tab.est) <- colnames(se) <- colnames(LCI) <- colnames(UCI) <- paste("q = ", q, sep="")
  rownames(tab.est) <- rownames(se) <- rownames(LCI) <- rownames(UCI) <- c("Observed", "Chao_2013")
  return(list(EST = tab.est,
              SD = se,
              LCI = LCI,
              UCI = UCI))

}


#----------------------------
# Plot of confidence interval
#----------------------------
#' conf.reg(x_axis,LCL,UCL,...) is a function to plot the confidence region.
#'
#' \code{conf.reg} uses polygon to draw a confidence band plot
#'
#' @param x_axis a vector of diversity orders.
#' @param LCL a vector of lower confidence bounds.
#' @param UCL a vector of upper confidence bounds.
#' @param ... further arguments to be passed to \code{polygon}
#' @return a polygon plot

conf.reg=function(x_axis,LCL,UCL,...) {
  x.sort <- order(x_axis)
  x <- x_axis[x.sort]
  LCL <- LCL[x.sort]
  UCL <- UCL[x.sort]
  polygon(c(x,rev(x)),c(LCL,rev(UCL)), ...)
}

#=====================================================================
#Apply Chao Hill to all the traps Diptera and EPT

# ========== PARALLEL PROCESSING SETUP ==========
plan(multisession, workers = 10)  # Use 10 of your 20 available threads

# Monitor progress and increase memory limits
options(future.globals.maxSize = 2000 * 1024^2)  # 2GB limit for data transfer

# Define the traps of interest
traps <- c("A", "B", "C", "E", "G", "I", "III")

# Initialize a list to store incidence vectors for each trap level
incidence_vectors <- list()

# Loop over each specified trap level
for (trap_level in traps) {
  # Subset the data for the current trap level
  diptera_subset <- diptera_long %>%
    filter(trap == trap_level)

  # Build the incidence matrix for the current trap level
  incidence_matrix <- diptera_subset %>%
    filter(abundance > 0) %>%
    distinct(sp_name, year) %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = year, values_from = present, values_fill = 0)

  # Count in how many years each species was present
  incidence_freqs <- rowSums(incidence_matrix[,-1])  # exclude sp_name

  # Count number of years (i.e., number of sampling units)
  number_of_years <- ncol(incidence_matrix) - 1  # exclude sp_name

  # Final incidence vector for the current trap level
  incidence_vector <- c(number_of_years, incidence_freqs)

  # Store the incidence vector in the list
  incidence_vectors[[trap_level]] <- incidence_vector
}

# The incidence_vectors list now contains an incidence vector for each specified trap level
# You can access them using the trap level names, e.g., incidence_vectors[["A"]]

# # Initialize a list to store the results of ChaoHill for each trap level
# chaohill_results <- list()
#
# # Loop over each trap level in the incidence_vectors list
# for (trap_level in names(incidence_vectors)) {
#   # Get the incidence vector for the current trap level
#   incidence_vector <- incidence_vectors[[trap_level]]
#
#   # Apply the ChaoHill function to the incidence vector
#   chaohill_result <- ChaoHill(incidence_vector)
#
#   # Store the result in the chaohill_results list
#   chaohill_results[[trap_level]] <- chaohill_result
# }

# # Parallel processing for ChaoHill - Diptera
# chaohill_results <- future_lapply(incidence_vectors, function(incidence_vector) {
#   ChaoHill(incidence_vector)  # Keeping B=1000 (default)
# }, future.seed = TRUE)
#
# # The chaohill_results list now contains the results of ChaoHill for each trap level
# # You can access them using the trap level names, e.g., chaohill_results[["A"]]s

# Clean up parallel processing
plan(sequential)  # Return to sequential processing

# saveRDS(chaohill_results, "Outputs/chaohill_results_diptera.rds")  # Save results for Diptera

chaohill_results <- readRDS(chaohill_results, "Outputs/chaohill_results_diptera.rds")  # read results for Diptera for quicker analyses

######################################################################################
# EPT ChaoHill Per trap

# ========== PARALLEL PROCESSING SETUP ==========
plan(multisession, workers = 10)  # Use 10 of your 20 available threads

# Monitor progress and increase memory limits
options(future.globals.maxSize = 2000 * 1024^2)  # 2GB limit for data transfer

# Define the traps of interest, excluding "D"
traps <- c("A", "B", "C", "E", "G", "I", "III")

# Initialize a list to store incidence vectors for each trap level
incidence_vectors <- list()

# Loop over each specified trap level
for (trap_level in traps) {
  # Subset the data for the current trap level
  ept_subset <- ept_long %>%
    filter(trap == trap_level)

  # Build the incidence matrix for the current trap level
  incidence_matrix <- ept_subset %>%
    filter(abundance > 0) %>%
    distinct(sp_name, year) %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = year, values_from = present, values_fill = 0)

  # Count in how many years each species was present
  incidence_freqs <- rowSums(incidence_matrix[,-1])  # exclude sp_name

  # Count number of years (i.e., number of sampling units)
  number_of_years <- ncol(incidence_matrix) - 1  # exclude sp_name

  # Final incidence vector for the current trap level
  incidence_vector <- c(number_of_years, incidence_freqs)

  # Store the incidence vector in the list
  incidence_vectors[[trap_level]] <- incidence_vector
}

# The incidence_vectors list now contains an incidence vector for each specified trap level
# You can access them using the trap level names, e.g., incidence_vectors[["A"]]

# # Initialize a list to store the results of ChaoHill for each trap level
# ept_chaohill_results <- list()
#
# # Loop over each trap level in the incidence_vectors list
# for (trap_level in names(incidence_vectors)) {
#   # Get the incidence vector for the current trap level
#   incidence_vector <- incidence_vectors[[trap_level]]
#
#   # Apply the ChaoHill function to the incidence vector
#   chaohill_result <- ChaoHill(incidence_vector)
#
# #   # Store the result in the ept_chaohill_results list
# #   ept_chaohill_results[[trap_level]] <- chaohill_result
# # }
#
# # Parallel processing for ChaoHill - EPT
# ept_chaohill_results <- future_lapply(incidence_vectors, function(incidence_vector) {
#   ChaoHill(incidence_vector)  # Keeping B=1000 (default)
# }, future.seed = TRUE)

# The ept_chaohill_results list now contains the results of ChaoHill for each trap level
# You can access them using the trap level names, e.g., ept_chaohill_results[["A"]]

# Clean up parallel processing
plan(sequential)  # Return to sequential processing

# saveRDS(ept_chaohill_results, "Outputs/chaohill_results_ept.rds")  # Save results for EPT

ept_chaohill_results <- readRDS("Outputs/chaohill_results_ept.rds")  # read results for EPT for quicker analyses

####################################################
png("Plots/Hill_numbers_long.png", width = 8, height = 20, unit = "in", res = 300, bg = "white")

q <- seq(0, 3, 0.1)

# Define colors to match your aesthetic
diptera_dark <- "#8B0000"      # Very dark red (proposed)
diptera_light <- "#CD5C5C"     # Lighter red (empirical)
ept_dark <- "#000080"          # Very dark blue (proposed)
ept_light <- "#4169E1"         # Lighter blue (empirical)

# Get trap names (assuming they match in order)
trap_names <- names(chaohill_results)
num_traps <- length(trap_names)

# Set up the plotting layout: 7 rows, 2 columns
par(mfrow = c(num_traps, 2),
    mai = c(0.1, 0.3, 0.1, 0.2),  # Increased bottom margin, reduced horizontal margins
    oma = c(7, 5, 3, 0.5))          # Increased outer margins: bottom, left, top, right

# Create plots for each trap
for (i in 1:num_traps) {
  trap_level <- trap_names[i]

  # DIPTERA PLOT (left column)
  dipt_inc <- chaohill_results[[trap_level]]
  ymin <- min(dipt_inc$LCI)
  ymax <- max(dipt_inc$UCI)

  # Determine if this is the bottom row for x-axis labels
  show_x_labels <- (i == num_traps)
  # Determine if this is the top row for column header
  show_col_header <- (i == 1)

  plot(q, dipt_inc$EST[1, ], type = "n",
       xlab = "", ylab = "",
       lty = 2, col = diptera_light, ylim = c(ymin, ymax),
       main = "",
       xaxt = ifelse(show_x_labels, "s", "n"),
       cex.axis = 1.5,
       las = 1)

  # Add x-axis tick marks without labels for non-bottom rows
  if (!show_x_labels) {
    axis(1, labels = FALSE, cex.axis = 1.5)
  }

  # Add confidence intervals
  polygon(c(q, rev(q)), c(dipt_inc$LCI[1, ], rev(dipt_inc$UCI[1, ])),
          border = NA, col = adjustcolor(diptera_light, 0.2))
  polygon(c(q, rev(q)), c(dipt_inc$LCI[2, ], rev(dipt_inc$UCI[2, ])),
          border = NA, col = adjustcolor(diptera_dark, 0.2))

  # Add main lines
  lines(q, dipt_inc$EST[2, ], type = "l", col = diptera_dark, lwd = 2)
  lines(q, dipt_inc$EST[1, ], type = "l", col = diptera_light, lty = 2, lwd = 2)

  # Add trap letter in top-right corner
  text(max(q) * 0.995, ymax * 0.98, trap_level,
       cex = 2.5, font = 2, adj = c(1, 1))

  # Add column header for first row
  if (show_col_header) {
    mtext("DIPTERA", side = 3, line = 1, cex = 1.5, font = 2)
  }

  # Add box around plot
  box(lwd = 1.2)

  # EPT PLOT (right column)
  ept_inc <- ept_chaohill_results[[trap_level]]
  ymin <- min(ept_inc$LCI)
  ymax <- max(ept_inc$UCI)

  plot(q, ept_inc$EST[1, ], type = "n",
       xlab = "", ylab = "",
       lty = 2, col = ept_light, ylim = c(ymin, ymax),
       main = "",
       xaxt = ifelse(show_x_labels, "s", "n"),
       cex.axis = 1.5,
       las = 1)

  # Add x-axis tick marks without labels for non-bottom rows
  if (!show_x_labels) {
    axis(1, labels = FALSE, cex.axis = 1.5)
  }

  # Add confidence intervals
  polygon(c(q, rev(q)), c(ept_inc$LCI[1, ], rev(ept_inc$UCI[1, ])),
          border = NA, col = adjustcolor(ept_light, 0.2))
  polygon(c(q, rev(q)), c(ept_inc$LCI[2, ], rev(ept_inc$UCI[2, ])),
          border = NA, col = adjustcolor(ept_dark, 0.2))

  # Add main lines
  lines(q, ept_inc$EST[2, ], type = "l", col = ept_dark, lwd = 2)
  lines(q, ept_inc$EST[1, ], type = "l", col = ept_light, lty = 2, lwd = 2)

  # Add trap letter in top-right corner
  text(max(q) * 0.995, ymax * 0.98, trap_level,
       cex = 2.5, font = 2, adj = c(1, 1))

  # Add column header for first row
  if (show_col_header) {
    mtext("EPT", side = 3, line = 1, cex = 1.5, font = 2)
  }

  # Add box around plot
  box(lwd = 1.2)
}

# Add single overall axis labels
mtext("Hill numbers", side = 2, outer = TRUE, line = 2, cex = 1.5, font = 2)
mtext("Order q", side = 1, outer = TRUE, line = 3, cex = 1.5, font = 2)

# Add a single legend below the axis titles
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')

# Create legend positioned lower to be below axis titles
legend("bottom",
       legend = c("Proposed", "Empirical"),
       col = c("black", "black"),
       lty = c(1, 2),
       lwd = 2,
       bty = "n",
       horiz = TRUE,
       cex = 1.5,  # 14pt text
       x.intersp = 2,
       text.width = 0.15,
       inset = c(0, 0.001))  # Increased inset to push legend further down

dev.off()

####################################################
# For plot wider format
png("Plots/Hill_numbers_wide.png", width = 20, height = 7, unit = "in", res = 300, bg = "white")

q <- seq(0, 3, 0.1)

# Define colors to match your aesthetic
diptera_dark <- "#8B0000"      # Very dark red (proposed)
diptera_light <- "#CD5C5C"     # Lighter red (empirical)
ept_dark <- "#000080"          # Very dark blue (proposed)
ept_light <- "#4169E1"         # Lighter blue (empirical)

# Get trap names (assuming they match in order)
trap_names <- names(chaohill_results)
num_traps <- length(trap_names)

# Set up the plotting layout: 2 rows, 7 columns
par(mfrow = c(2, num_traps),
    mai = c(0.1, 0.4, 0.1, 0.1),  # Increased left margin for y-axis space
    oma = c(7.5, 5, 3, 0.5))        # Increased outer margins: bottom, left, top, right

# Create plots - First row: DIPTERA, Second row: EPT
for (i in 1:num_traps) {
  trap_level <- trap_names[i]

  # DIPTERA PLOT (top row)
  dipt_inc <- chaohill_results[[trap_level]]
  ymin <- min(dipt_inc$LCI)
  ymax <- max(dipt_inc$UCI)

  # Show y-axis labels on all plots since they have unique scales
  show_row_header <- (i == 1)

  plot(q, dipt_inc$EST[1, ], type = "n",
       xlab = "", ylab = "",
       lty = 2, col = diptera_light, ylim = c(ymin, ymax),
       main = "",
       xaxt = "n",  # No x-axis labels on top row
       cex.axis = 1.5,
       las = 1)

  # Add x-axis tick marks without labels
  axis(1, labels = FALSE, cex.axis = 1.5)

  # Add confidence intervals
  polygon(c(q, rev(q)), c(dipt_inc$LCI[1, ], rev(dipt_inc$UCI[1, ])),
          border = NA, col = adjustcolor(diptera_light, 0.2))
  polygon(c(q, rev(q)), c(dipt_inc$LCI[2, ], rev(dipt_inc$UCI[2, ])),
          border = NA, col = adjustcolor(diptera_dark, 0.2))

  # Add main lines
  lines(q, dipt_inc$EST[2, ], type = "l", col = diptera_dark, lwd = 2)
  lines(q, dipt_inc$EST[1, ], type = "l", col = diptera_light, lty = 2, lwd = 2)

  # Add column header (trap name)
  mtext(trap_level, side = 3, line = 1, cex = 1.5, font = 2)

  # Add row header for first column
  if (show_row_header) {
    mtext("DIPTERA", side = 2, line = 3.5, cex = 1.5, font = 2)
  }

  # Add box around plot
  box(lwd = 1.2)
}

# Second row: EPT plots
for (i in 1:num_traps) {
  trap_level <- trap_names[i]

  # EPT PLOT (bottom row)
  ept_inc <- ept_chaohill_results[[trap_level]]
  ymin <- min(ept_inc$LCI)
  ymax <- max(ept_inc$UCI)

  # Show y-axis labels on all plots since they have unique scales
  show_row_header <- (i == 1)

  plot(q, ept_inc$EST[1, ], type = "n",
       xlab = "", ylab = "",
       lty = 2, col = ept_light, ylim = c(ymin, ymax),
       main = "",
       xaxt = "s",  # Show x-axis labels on bottom row
       cex.axis = 1.5,
       las = 1)

  # Add confidence intervals
  polygon(c(q, rev(q)), c(ept_inc$LCI[1, ], rev(ept_inc$UCI[1, ])),
          border = NA, col = adjustcolor(ept_light, 0.2))
  polygon(c(q, rev(q)), c(ept_inc$LCI[2, ], rev(ept_inc$UCI[2, ])),
          border = NA, col = adjustcolor(ept_dark, 0.2))

  # Add main lines
  lines(q, ept_inc$EST[2, ], type = "l", col = ept_dark, lwd = 2)
  lines(q, ept_inc$EST[1, ], type = "l", col = ept_light, lty = 2, lwd = 2)

  # Add row header for first column
  if (show_row_header) {
    mtext("EPT", side = 2, line = 3.5, cex = 1.5, font = 2)
  }

  # Add box around plot
  box(lwd = 1.2)
}

# Add single overall axis labels
mtext("Hill numbers", side = 2, outer = TRUE, line = 2, cex = 1.5, font = 2)
mtext("Order q", side = 1, outer = TRUE, line = 3, cex = 1.5, font = 2)

# Add a single legend below the axis titles
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')

# Create legend positioned lower to be below axis titles
legend("bottom",
       legend = c("Proposed", "Empirical"),
       col = c("black", "black"),
       lty = c(1, 2),
       lwd = 2,
       bty = "n",
       horiz = TRUE,
       cex = 1.5,
       x.intersp = 2,
       text.width = 0.15,
       inset = c(0, 0.005))  # Better positioned legend

dev.off()
