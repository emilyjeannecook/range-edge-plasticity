#Author: Emily Cook
#Project: PERSIST plasticity
#Purpose: Calculate phenotypic plasticity for first flower date, 
#specific leaf area (SLA), and leaf dry matter content (LDMC)
#end product: three dataframes for each trait plasticity

#packages####
library(tidyverse)
library(lme4)
library(lmerTest)
library(tictoc)
library(brms)
library(tidybayes)
library(emmeans)

# Assign colors per population
cols <- c("N1" = "#5E4FA2", 
          "N2" = "dodgerblue4", 
          "C1" = "gold", 
          "C2" = "goldenrod1", 
          "S1" = "red2", 
          "S2" = "#9E0142")



#setwd
##setwd("")

#upload data####
#dat <- read.csv("PERSIST_2023_data.csv") from DRYAD: https://doi.org/10.5061/dryad.18931zd7g


######## COLUMN DESCRIPTIONS ####

# garden            : experimental garden (north, central, south)
# block             : experimental randomized block in garden (1 - 10 in north, 1 - 11 in central, 1 - 10 in south)
# garden_block      : variable that combines garden and block
# row               : row (y-coordinate) of experimental garden; with 4 rows per block; rows 101 - 104 are in block 1; rows 1101 - 1104 are in block 11, etc.
# position          : position (x-coordinate) of experimental garden; corresponds to a unique plant ID (Cross_ID_Rep), or has no plant (NA)
# rowPosition       : variable that combines row and position
# Cross_ID_Rep      : variable that combines unique ID for each full-sibling family and replicate of that family within a particular garden 
# Cross_ID          : unique ID for each full-sibling family; each Cross_ID has a unique mom and dad 
# Sire_ID           : unique ID for sire (father); plants with the same Sire_ID and different Dam_IDs are half-sibs
# Dam_ID            : unique ID for sire (mother); dams are nested within sires to yield a nested half-sib/full-sib design
# Population        : Population of scarlet monkeyflower (N1 and N2: northern populations; C1 and C2: central populations; S1 and S2: southern populations)
# Year              : Year that seeds were collected in the field (2010 ancestors and 2017 descendants)
# Year1             : Alternate coding for year corresponding to "ancestor" and "descendant"
# Date_early        : Date of early-season li-600 data collection 
# Time1_early       : Time of early-season li-600 data collection, centered and scaled 
# log_VPDleaf_early : Log-transformed leaf vapor pressure deficit at the time of early-season li-600 data collection
# gsw_early         : Early-season stomatal conductance to water vapor, measured at the leaf level with a li-600 porometer in units of mmol/m²/s
# freshMass_g       : Fresh leaf mass in grams 
# dryMass_g         : Oven-dried leaf mass in grams 
# leafArea_cm2      : Leaf area in square centimeters, derived from leaf scans
# lma_g_per_m2      : Dry leaf mass in grams per area in meters squared 
# sla_cm2_per_g     : Specific leaf area (leaf area in squared centimeters divided by dry leaf mass in grams)
# ldmc              : Leaf dry matter content, measured as dry leaf mass in grams divided by fresh leaf mass in grams
# L1                : Length of the primary or longest stem at first flower in centimeters
# L2                : Length of the second longest stem at first flower in centimeters
# L3                : Length of the third longest stem at first flower in centimeters
# totalStemLen      : Sum of the lengths of the three longest stems at first flower in centimeters
# first_flower_date : Date of first flower
# first_flower_doy  : Day of year of first flower
# last_flower_date  : Date of last flower; note this is not reliable because we did not continue collecting data after a certain point in the growing season
# last_flower_doy   : Day of year of last flower; note this is not reliable because we did not continue collecting data after a certain point in the growing season
# flowering_duration: Duration of flowering expressed as the difference between the date of last flower and the date of first flower; note this is not reliable because we did not continue collecting data after a certain point in the growing season
# Date_late         : Date of early-season li-600 data collection
# Time1_late        : Time of early-season li-600 data collection, centered and scaled (CDM pelase clarify) 
# log_VPDleaf_late  : Log-transformed leaf vapor pressure deficit at the time of late-season li-600 data collection
# gsw_late          : Late-season stomatal conductance to water vapor, measured at the leaf level with a li-600 porometer in units of mmol/m²/s
# maxHeight         : Maximum stem height in centimeters at the end of the growing season
# repBranchN        : Number of major reproductive branches at the end of the growing season
# RScount1          : Number of reproductive structures (flowers, fruits, buds, and pedicels) counted on the stem with the most reproductive structures at the end of the growing season
# RScount2          : Number of reproductive structures (flowers, fruits, buds, and pedicels) counted on a representative major reproductive branch at the end of the growing season
# RScount3          : Number of reproductive structures (flowers, fruits, buds, and pedicels) counted on a representative major reproductive branch at the end of the growing season
# totalRS           : An estimate of the total number of reproductive structures (flowers, fruits, buds, and pedicels) on a plant, calculated as RScount1 + (RScount2 + RScount3)/2 * (repBranchN - 1); calculations were more complex when counts were incomplete



#workflow for first flower date####
#1. find number of full-sib familieswith at least 2 reps in each garden and subset dataframe for each trait separately
#2. get the sample size of number of families/individuals per pop x cohort
#3. "adjust" date of first flower: subtract each DOY by first date of flower in each garden
#4. visualize raw vs. adjusted first flower date
#5. calculate plasticity
#find families, and subset dataframe for each trait separately
#calculate which garden has largest value for each trait, since plasticity subtraction needs average largest value first
#For each population and year:
#take all possible differences and divide by the global mean of that trait
#allow negative values, genotypes which have negative plasticity are those that go in the opposite direction of the average
#i.e. if plants in north on average flower later, subtract north - south, if the value is negative, that means south flowered later than north
#6. take family level averages of those differences for family level plasticity value

#first flower date plasticity####
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
  )

samp_size_ff
#save table


#3.####
#3. "adjust" date of first flower: subtract each DOY by first date of flower in each garden

subset_phen <- subset_phen %>%
  group_by(garden) %>%
  mutate(transplant_DOY = first_flower_doy - min(first_flower_doy)) %>%
  ungroup()

#set factors
subset_phen <- subset_phen %>% mutate(Population=fct_relevel(Population, "N1","N2","C1", "C2", "S1","S2" )) #relevel
subset_phen <- subset_phen %>% mutate(garden=fct_relevel(garden,"north","center","south")) #relevel

#make a pop-year combined column
subset_phen <- subset_phen%>%
  mutate(pop_year=paste(Population,"-",Year, sep='')) #add additional pop-year column



#4. ####
#Visualize raw vs. adjusted first flower date
#adjusted first flower date across gardens
ggplot(subset_phen, aes(x = Population, y = first_flower_doy, fill = Population)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~garden)+
  geom_jitter(width = 0.2, shape = 21, color = "black", size = 1) +
  #geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Source Population", y = "Adjusted First Flower", alpha = "Year") +
  scale_fill_manual(values = cols) +
  theme_bw()

#raw flower date across gardens
ggplot(subset_phen, aes(x = Population, y = transplant_DOY, fill = Population)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~garden)+
  geom_jitter(width = 0.2, shape = 21, color = "black", size = 1) +
  #geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Source Population", y = "Adjusted First Flower", alpha = "Year") +
  scale_fill_manual(values = cols) +
  theme_bw()


#5.####
#5. calculate plasticity

#plasticity is the average of all the possible differences 
#between the three gardens, each divided by the global mean adjusted day of year

#Find which garden mean is largest
#To calculate plasticity, the differences must have the greater value first
DOY_garden_mean <- subset_phen %>%
  group_by(garden)%>%
  summarise(DOY_mean = mean(transplant_DOY)) #using transplant date
DOY_garden_mean

#north = ~30
#south = ~21
#central = ~18
#so on average, you would expect plants in the north garden to flower latest, south intermediate, central earliest
#using the larger value first use this difference pattern:
#N-S
#S-C
#N-C

#First flower date plasticity
#will loop through each pop-year combo:
pop_year <- c("N1-2010", "N1-2017", "N2-2010", "N2-2017", 
              "S1-2010", "S1-2017", "S2-2010", "S2-2017",
              "C1-2010", "C2-2010", "C1-2017", "C2-2017")

#Create an empty dataframe to store results
doy_pp <- data.frame(Dam_ID = character(), 
                     pop_year = character(),
                     garden_1 = character(), #subtract garden_1 - garden_2
                     garden_2 = character(),
                     doy_1 = numeric(),
                     doy_2 = numeric(),
                     difference = numeric(),
                     ratio = numeric(),
                     stringsAsFactors = FALSE)

#Loop through each pop-year combinations
for (pop in pop_year) {
  
  #Filter for the current pop-year
  pop_subset <- subset(subset_phen, pop_year == pop)
  
  #Filter subsets for each garden
  north_garden <- pop_subset[pop_subset$garden == "north", ]
  south_garden <- pop_subset[pop_subset$garden == "south", ]
  central_garden <- pop_subset[pop_subset$garden == "center", ]
  
  #Count replicates in each garden
  north_counts <- north_garden %>% group_by(Dam_ID) %>% summarise(north_reps = n())
  south_counts <- south_garden %>% group_by(Dam_ID) %>% summarise(south_reps = n())
  central_counts <- central_garden %>% group_by(Dam_ID) %>% summarise(central_reps = n())
  
  #Merge the counts
  combined_counts <- reduce(
    list(north_counts, south_counts, central_counts),
    ~ merge(.x, .y, by = "Dam_ID", all = TRUE)
  )
  
  #Filter for valid families (Dam_IDs)
  valid_sires <- combined_counts %>%
    filter(north_reps >= 2 & south_reps >= 2 & central_reps >= 2)
  
  #Loop through valid families ( `Dam_ID` )
  for (sire in valid_sires$Dam_ID) {
    sub_df <- subset(pop_subset, Dam_ID == sire)
    
    #Extract `transplant_DOY` for each garden
    north_doy <- sub_df$transplant_DOY[sub_df$garden == "north"]
    south_doy <- sub_df$transplant_DOY[sub_df$garden == "south"]
    central_doy <- sub_df$transplant_DOY[sub_df$garden == "center"]
    
    if (length(north_doy) > 0 && length(south_doy) > 0 && length(central_doy) > 0) {
      
      #Generate pairwise combinations
      north_south <- expand.grid(north_doy, south_doy) #north minus south
      north_central <- expand.grid(north_doy, central_doy) #north minus central
      south_central <- expand.grid(south_doy, central_doy) #south minus central
      
      #Calculate the global mean
      global_mean <- mean(subset_phen$transplant_DOY, na.rm = TRUE)
      
      #Append pairwise comparisons
      doy_pp <- rbind(
        doy_pp,
        data.frame(
          Dam_ID = sire,
          pop_year = pop,
          garden_1 = "north",
          garden_2 = "south",
          doy_1 = north_south$Var1, #default names from expand.grid()
          doy_2 = north_south$Var2,
          difference = north_south$Var1 - north_south$Var2,
          ratio = (north_south$Var1 - north_south$Var2) / global_mean
        ),
        data.frame(
          Dam_ID = sire,
          pop_year = pop,
          garden_1 = "north",
          garden_2 = "center",
          doy_1 = north_central$Var1,
          doy_2 = north_central$Var2,
          difference = north_central$Var1 - north_central$Var2,
          ratio = (north_central$Var1 - north_central$Var2) / global_mean
        ),
        data.frame(
          Dam_ID = sire,
          pop_year = pop,
          garden_1 = "south",
          garden_2 = "center",
          doy_1 = south_central$Var1,
          doy_2 = south_central$Var2,
          difference = south_central$Var1 - south_central$Var2,
          ratio = (south_central$Var1 - south_central$Var2) / global_mean
        )
      )
    }
  }
}

#6. ####
#Get family level average plasticity values

#group by family (Dam_ID) to take average of ratio column
doy_pp_mn <-  doy_pp %>% 
  group_by(Dam_ID, pop_year)%>% 
  summarize(global_diff = mean(ratio)) 

#split pop-year into population and year
doy_pp_mn <- doy_pp_mn %>%
  separate(pop_year, into = c("Population", "Year"), sep = "-", 
           remove = FALSE)

#calculate average days shifted in days of year DOY by taking the:
#difference (days)=plasticity×global mean DOY
mean_plasticity <- mean(doy_pp_mn$global_diff) #mean plasticity value
mean_DOY <- mean(subset_phen$transplant_DOY, na.rm = TRUE) #mean DOY
mean_plasticity * mean_DOY #8.16 days





#leaf traits####

#workflow for SLA & LDMC####
#1. find number of full-sib familieswith at least 2 reps in each garden and subset dataframe for each trait separately
#2. get the sample size of number of families/individuals per pop x cohort
#3. visualize SLA & LDMC
#4. calculate plasticity
#find families, and subset dataframe for each trait separately
#calculate which garden has largest value for each trait, since plasticity subtraction needs average largest value first
#For each population and year:
#take all possible differences and divide by the global mean of that trait
#allow negative values, genotypes which have negative plasticity are those that go in the opposite direction of the average
#i.e. if plants in north on average flower later, subtract north - south, if the value is negative, that means south flowered later than north
#5. take family level averages of those differences for family level plasticity value

#SLA plasticity####
#1.####
#1. find number of full-sib familieswith at least 2 reps in each garden and subset dataframe for each trait separately
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

#add pop_year to subset_leaf
subset_leaf <- subset_leaf%>%
  mutate(pop_year=paste(Population,"-",Year, sep='')) #add additional pop-year column

#2.####
#2. get the sample size of number of 
#families and the sample size of number of individuals per pop x cohort
samp_size_sla <- subset_leaf %>%
  group_by(Population, Year) %>%
  summarise(
    n_families   = n_distinct(Dam_ID),
    n_individuals = n(),
    .groups = "drop"
  )

samp_size_sla
#save table


#3. ####
#Visualize SLA across gardens
ggplot(subset_leaf, aes(x = Population, y = sla_cm2_per_g, fill = Population)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~garden)+
  geom_jitter(width = 0.2, shape = 21, color = "black", size = 1) +
  #geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Source Population", y = "SLA", alpha = "Year") +
  scale_fill_manual(values = cols) +
  theme_bw()

#4.####
#4. calculate plasticity
#Find which garden mean is largest
#To calculate plasticity, the differences must have the greater value first
SLA_garden_mean <- subset_leaf %>%
  group_by(garden)%>%
  summarise(SLA_mean = mean(sla_cm2_per_g, na.rm=TRUE)) 

#north = ~242
#south = ~320
#central = NA, no data....
#using the larger value first use this difference pattern:
#S-N

#SLA plasticity
#will loop through each pop-year combo:
pop_year <- c("N1-2010", "N1-2017", "N2-2010", "N2-2017", 
              "S1-2010", "S1-2017", "S2-2010", "S2-2017",
              "C1-2010", "C2-2010", "C1-2017", "C2-2017")

#Initialize an empty dataframe to store results
sla_pp <- data.frame(Dam_ID = character(), 
                     pop_year = character(),
                     garden_1 = character(),
                     garden_2 = character(),
                     doy_1 = numeric(),
                     doy_2 = numeric(),
                     difference = numeric(),
                     ratio = numeric(),
                     stringsAsFactors = FALSE)

#Loop through pop-year combinations
for (pop in pop_year) {
  
  #Filter for the current pop-year
  pop_subset <- subset(subset_leaf, pop_year == pop)
  
  #Filter subsets for each garden
  north_garden <- pop_subset[pop_subset$garden == "north", ]
  south_garden <- pop_subset[pop_subset$garden == "south", ]
  #central_garden <- pop_subset[pop_subset$garden == "center", ]
  
  #Count replicates in each garden
  north_counts <- north_garden %>% group_by(Dam_ID) %>% summarise(north_reps = n())
  south_counts <- south_garden %>% group_by(Dam_ID) %>% summarise(south_reps = n())
  #central_counts <- central_garden %>% group_by(Dam_ID) %>% summarise(central_reps = n())
  
  #Merge counts
  combined_counts <- reduce(
    list(north_counts, south_counts),
    ~ merge(.x, .y, by = "Dam_ID", all = TRUE)
  )
  
  #Filter for valid families (Dam_IDs)
  valid_sires <- combined_counts %>%
    filter(north_reps >= 2 & south_reps >= 2)
  
  #Loop through valid `Dam_ID`
  for (sire in valid_sires$Dam_ID) {
    sub_df <- subset(pop_subset, Dam_ID == sire)
    
    #Extract SLA for each garden
    north_doy <- sub_df$sla_cm2_per_g[sub_df$garden == "north"] 
    south_doy <- sub_df$sla_cm2_per_g[sub_df$garden == "south"]
   #central_doy <- sub_df$transplant_DOY[sub_df$garden == "center"]
    
    if (length(north_doy) > 0 && length(south_doy) > 0)  {
      #Generate pairwise combinations
        #this is where i control direction of subtraction
        #in this case, listing south first means its the first variable
        #so it will go first in the subtraction
      north_south <- expand.grid(south_doy, north_doy) #south minus north
        #north_central <- expand.grid(north_doy, central_doy) #north minus central
        #south_central <- expand.grid(south_doy, central_doy) #south minus central
      
      #Calculate global mean
      global_mean <- mean(subset_leaf$sla_cm2_per_g, na.rm = TRUE)
      
      #Append pairwise comparisons
      sla_pp <- rbind(
        sla_pp,
        data.frame(
          Dam_ID = sire,
          pop_year = pop,
          garden_1 = "south",#even if i change the order of these, you still need to change
            #the above expand.grid() to actually change the direction of subtraction!
          garden_2 = "north",
          doy_1 = north_south$Var1, #default names from expand.grid()
          doy_2 = north_south$Var2,
          difference = north_south$Var1 - north_south$Var2,
          ratio = (north_south$Var1 - north_south$Var2) / global_mean
        )
      )
    }
  }
}


#5. ####
#Get family level average plasticity values

#group by family (Dam_ID) to take average of ratio column
sla_pp_mn <-  sla_pp %>% #call data frame
  group_by(Dam_ID, pop_year)%>% 
  summarize(global_diff = mean(ratio, na.rm=TRUE))#sum all ratios/divide number of distances

#split pop-year into population and year
sla_pp_mn <- sla_pp_mn %>%
  separate(pop_year, into = c("Population", "Year"), sep = "-", 
           remove = FALSE)





#LDMC plasticity####

#1.####
#1. find number of full-sib familieswith at least 2 reps in each garden and subset dataframe for each trait separately
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
  )

samp_size_ldmc
#save table

#3. ####
#Visualize ldmc across gardens
ggplot(subset_ldmc, aes(x = Population, y = ldmc, fill = Population)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~garden)+
  geom_jitter(width = 0.2, shape = 21, color = "black", size = 1) +
  #geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Source Population", y = "LDMC", alpha = "Year") +
  scale_fill_manual(values = cols) +
  theme_bw()

#4.####
#4. calculate plasticity
#Find which garden mean is largest
#To calculate plasticity, the differences must have the greater value first
ldmc_garden_mean <- subset_ldmc %>%
  group_by(garden)%>%
  summarise(ldmc_mean = mean(ldmc, na.rm=TRUE)) #using transplant date

#north = ~0.17
#south = ~012
#central = NA, no data....
#using the larger value first use this difference pattern:
#N-S

#SLA plasticity

# Initialize an empty dataframe to store results
ldmc_pp <- data.frame(Dam_ID = character(), 
                     pop_year = character(),
                     garden_1 = character(),
                     garden_2 = character(),
                     doy_1 = numeric(),
                     doy_2 = numeric(),
                     difference = numeric(),
                     ratio = numeric(),
                     stringsAsFactors = FALSE)

#Loop through pop-year combinations
for (pop in pop_year) {
  
  #Filter for the current pop-year
  pop_subset <- subset(subset_ldmc, pop_year == pop)
  
  #Filter subsets for each garden
  north_garden <- pop_subset[pop_subset$garden == "north", ]
  south_garden <- pop_subset[pop_subset$garden == "south", ]
  #central_garden <- pop_subset[pop_subset$garden == "center", ]
  
  #Count replicates in each garden
  north_counts <- north_garden %>% group_by(Dam_ID) %>% summarise(north_reps = n())
  south_counts <- south_garden %>% group_by(Dam_ID) %>% summarise(south_reps = n())
  #central_counts <- central_garden %>% group_by(Dam_ID) %>% summarise(central_reps = n())
  
  #Merge counts
  combined_counts <- reduce(
    list(north_counts, south_counts),
    ~ merge(.x, .y, by = "Dam_ID", all = TRUE)
  )
  
  #Filter for valid families (Dam_IDs)
  valid_sires <- combined_counts %>%
    filter(north_reps >= 2 & south_reps >= 2)
  
  #Loop through valid `Dam_ID`
  for (sire in valid_sires$Dam_ID) {
    sub_df <- subset(pop_subset, Dam_ID == sire)
    
    #Extract `ldmc` for each garden
    north_doy <- sub_df$ldmc[sub_df$garden == "north"]
    south_doy <- sub_df$ldmc[sub_df$garden == "south"]
    # central_doy <- sub_df$transplant_DOY[sub_df$garden == "center"]
    
    if (length(north_doy) > 0 && length(south_doy) > 0)  {
      # Generate pairwise combinations
        #this is where i control direction of subtraction
        #in this case, listing south first means its the first variable
        #so it will go first in the subtraction
      north_south <- expand.grid(north_doy, south_doy) #north minus south
        #north_central <- expand.grid(north_doy, central_doy) #north minus central
        # south_central <- expand.grid(south_doy, central_doy) #south minus central
      
      #Calculate global mean
      global_mean <- mean(subset_ldmc$ldmc, na.rm = TRUE)
      
      #Append pairwise comparisons
      ldmc_pp <- rbind(
        ldmc_pp,
        data.frame(
          Dam_ID = sire,
          pop_year = pop,
          garden_1 = "north",#even if i change the order of these, you still need to change
          #the above expand.grid() to actually change the direction of subtraction!
          garden_2 = "south",
          doy_1 = north_south$Var1, #default names from expand.grid()
          doy_2 = north_south$Var2,
          difference = north_south$Var1 - north_south$Var2,
          ratio = (north_south$Var1 - north_south$Var2) / global_mean
        )
      )
    }
  }
}




#5. ####
#Get family level average plasticity values

#group by family (Dam_ID) to take average of ratio column
ldmc_pp_mn <-  ldmc_pp %>% 
  group_by(Dam_ID, pop_year)%>% 
  summarize(global_diff = mean(ratio, na.rm=TRUE))#sum all ratios/divide number of distances

ldmc_pp_mn <- ldmc_pp_mn %>%
  separate(pop_year, into = c("Population", "Year"), sep = "-", 
           remove = FALSE)




#save to csv####
#date of first flower
doy_pp_mn <- doy_pp_mn %>% 
  rename(
    doy_plasticity = global_diff)
#write.csv(doy_pp_mn,"doy_pp.csv", row.names = FALSE)#write
#doy_pp_mn <- read.csv("doy_pp.csv") #read

#SLA
sla_pp_mn <- sla_pp_mn %>% 
  rename(
   sla_plasticity = global_diff)
#write.csv(sla_pp_mn,"sla_pp.csv", row.names = FALSE) #write
#sla_pp_mn <- read.csv("sla_pp.csv") #read

#LDMC
ldmc_pp_mn <- ldmc_pp_mn %>% 
  rename(
    ldmc_plasticity = global_diff)
#write.csv(ldmc_pp_mn,"ldmc_pp.csv", row.names = FALSE) #write
#ldmc_pp_mn <- read.csv("ldmc_pp.csv") #read



