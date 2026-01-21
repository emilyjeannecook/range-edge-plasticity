#Author: Emily Cook
#Project: PERSIST plasticity
#Purpose: Supplementary Table S1, the number of full-sibling families 
#with at least 2 reps in each garden, and the sample size of individuals for each population, cohort, trait combo

library(tidyverse)

#load data
#dat <- read.csv("PERSIST_2023_data.csv") from DRYAD: https://doi.org/10.5061/dryad.18931zd7g


#first flower date####
#1.####
#1. find number of full-sib familieswith at least 2 reps in each garden and subset dataframe for each trait separately

#create table that counts # of families in each garden, and makes a list
table <- dat%>%
  filter(!is.na(first_flower_doy)) %>% #only include rows that have an entry
  group_by(garden, Dam_ID, Year)%>%
  summarise(n = n())  

#create table which sums which families have more than 2 in each garden
table2 <- table %>%  group_by(Dam_ID, Year) %>%
  summarise(N_north = length(which(n[garden=="north"]>=2)), #do sire/dams have more than 2 in north
            N_south = length(which(n[garden=="south"]>=2)), #do sire/dams have more than 2 in south
            N_central = length(which(n[garden=="center"]>=2)), #do sire/dams have more than 2 in south
            both = ifelse(N_north>0 & N_south>0 & N_central>0, 1,0)) # do BOTH have more than 2

phen_sum <- table2$Dam_ID[table2$both ==1] #list of families
table(table2$both) #how many total families?

#subset data to include only families within phen_sum (have at least 2 reps in each garden)
subset_phen <- dat %>%
  filter(Dam_ID %in% phen_sum, !is.na(first_flower_doy))

str(subset_phen)

#2.####
#2. get the sample size of number of 
#families and the sample size of number of individuals per pop x cohort
samp_size_ff <- subset_phen %>%
  group_by(Population, Year) %>%
  summarise(
    n_families   = n_distinct(Dam_ID),
    n_individuals = n(),
    .groups = "drop"
  )%>%
  mutate(trait = "First flower date plasticity")

samp_size_ff


#SLA####
#1.####
#1. find number of full-sib families with at least 2 reps in each garden and subset dataframe for each trait separately
#Filter out rows where sla_cm2 is NA
tableLF <- dat %>%
  filter(!is.na(sla_cm2_per_g)) %>%   # keep only rows with a valid sla_cm2
  group_by(garden, Dam_ID, Year) %>%
  summarise(n = n(), .groups = "drop")

#create table which sums which families have more than 2 in each garden
table2LF <- tableLF %>%
  group_by(Dam_ID, Year) %>%
  summarise(
    N_north = sum(garden == "north" & n >= 2),
    N_south = sum(garden == "south" & n >= 2),
    both    = if_else(N_north > 0 & N_south > 0, 1, 0),
    .groups = "drop"
  )
leaf_sum <- table2LF$Dam_ID[table2LF$both ==1]
table(table2LF$both) #768 total

#subset data to include only families within phen_sum (have at least 2 reps in each garden)
subset_leaf <- dat %>%
  filter(Dam_ID %in% leaf_sum, !is.na(sla_cm2_per_g))

#2.####
#2. get the sample size of number of 
#families and the sample size of number of individuals per pop x cohort
samp_size_sla <- subset_leaf %>%
  group_by(Population, Year) %>%
  summarise(
    n_families   = n_distinct(Dam_ID),
    n_individuals = n(),
    .groups = "drop"
  )%>%
  mutate(trait = "SLA plasticity")

samp_size_sla


#LDMC####
#1.####
#1. find number of full-sib families with at least 2 reps in each garden and subset dataframe for each trait separately
#Filter out rows where LDMC is NA
tableldmc <- dat %>%
  filter(!is.na(ldmc)) %>%   # keep only rows with a valid ldmc
  group_by(garden, Dam_ID, Year) %>%
  summarise(n = n(), .groups = "drop")

#create table which sums which families have more than 2 in each garden
table2ldmc <- tableldmc %>%
  group_by(Dam_ID, Year) %>%
  summarise(
    N_north = sum(garden == "north" & n >= 2),
    N_south = sum(garden == "south" & n >= 2),
    both    = if_else(N_north > 0 & N_south > 0, 1, 0),
    .groups = "drop"
  )
ldmc_sum <- table2ldmc$Dam_ID[table2ldmc$both ==1]
table(table2ldmc$both) #774total

subset_ldmc <- dat %>%
  filter(Dam_ID %in% ldmc_sum, !is.na(ldmc))

#add pop_year to subset_leaf
subset_ldmc <- subset_ldmc%>%
  mutate(pop_year=paste(Population,"-",Year, sep='')) #add additional pop-year column

#2.####
#2. get the sample size of number of 
#families and the sample size of number of individuals per pop x cohort
samp_size_ldmc <- subset_ldmc %>%
  group_by(Population, Year) %>%
  summarise(
    n_families   = n_distinct(Dam_ID),
    n_individuals = n(),
    .groups = "drop"
  )%>%
  mutate(trait = "LDMC plasticity")

samp_size_ldmc


#combine tables
TableS1 <- rbind(samp_size_ldmc, samp_size_sla, samp_size_ff)







