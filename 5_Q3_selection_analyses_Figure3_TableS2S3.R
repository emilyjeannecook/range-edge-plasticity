#Author: Emily Cook
#Project: PERSIST plasticity
#Purpose: 
#1. Model selection (linear vs quadratic models) with loo compare and across garden selection analyses
#2. Plot Supplemental Table S2 (loo compare output) and S3 (selection coefficients)
#3. Plot figure 3 using raw data and estimated relationship with fitness and plasticity using linear models
#    and plot supplemental figures S7, S8, S9 visualizing estimated effect of plasticity on average (across garden) fitness

#traits: first flower date plasticity, SLA plasticity, and LDMC plasticity


#packages####
library(brms)
library(tidybayes)
library(dplyr)
library(emmeans)
library(tidyverse)
library(lme4)
library(lmerTest)
library(tictoc)
library(GGally)
library(ggplot2)
library(purrr)
library(tibble)
library(posterior)


#theme
theme_set(theme_bw())

#colors
cols <- c("N1" = "#5E4FA2", 
          "N2" = "dodgerblue4", 
          "C1" = "gold", 
          "C2" = "goldenrod1", 
          "S1" = "red2", 
          "S2" = "#9E0142")


#load plasticity data from 2_calculate_plasticity.R
#doy_pp_mn <- read.csv(") #read
#sla_pp_mn <- read.csv("sla_pp.csv") #read
#ldmc_pp_mn <- read.csv("ldmc_pp.csv") #read

#load fitness data
#dat <- read.csv("PERSIST_2023_data.csv") from DRYAD: https://doi.org/10.5061/dryad.18931zd7g


#prep dataframes for selection analyses####
#join plasticity dataframes####
#ff, sla, ldmc
pp_traits <- left_join(doy_pp_mn, sla_pp_mn, by = "Dam_ID") %>%
  dplyr::select(Dam_ID,
         doy_plasticity,
         sla_plasticity,
         Population = Population.x, #both df's have Population and Year so need to rename
         Year = Year.x)

#full df with all plastic traits
pp_traits <- left_join(pp_traits, ldmc_pp_mn, by = "Dam_ID") %>%
  dplyr::select(Dam_ID,
         doy_pp = doy_plasticity,
         sla_pp = sla_plasticity,
         ldmc_pp = ldmc_plasticity,
         Population = Population.x, #both df's have Population and Year so need to rename
         Year = Year.x)


#create quadratic (squared) terms####
pp_traits <- pp_traits %>%
  mutate(
    doy_pp2 = doy_pp^2,
    sla_pp2 = sla_pp^2,
    ldmc_pp2 = ldmc_pp^2,
  )

#avg fit traits across gardens ####
full<- dat %>%
  group_by(Dam_ID, Population, Year)%>%
  summarize(mean_rs = mean(totalRS, na.rm=TRUE), #mean total repro structures by genotype
            #mean_doy = mean(transplant_doy, na.rm = TRUE), #mean adjusted DOY
            #mean_sla = mean(sla_cm2_per_g, na.rm = TRUE), #mean SLA
            #mean_ldmc = mean(ldmc, na.rm = TRUE), #mean LDMC
            .groups = "drop") 

#join average fitness with plasticity
full <- left_join(full, pp_traits, by = "Dam_ID") %>%
  dplyr::select(doy_pp, sla_pp, ldmc_pp,  
                doy_pp2, sla_pp2, ldmc_pp2, 
                mean_rs, 
                #mean_sla, mean_ldmc, mean_doy, 
                Population = Population.x, #both df's have Population and Year so need to rename
                Year = Year.x, 
                Dam_ID)

#relativize fitness, standardize traits ####
full_z <- full %>%
  mutate(across(c(doy_pp, sla_pp, ldmc_pp, #plasticity
                  doy_pp2, sla_pp2, ldmc_pp2, #quadratic terms
                  #mean_rs, fitness
                  #mean_sla, mean_ldmc, mean_doy #trait means
                  ),
                ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE), #standardize traits
                .names = "{.col}_z"),
         fitness_rel = mean_rs / mean(mean_rs, na.rm = TRUE))#Relativize fit by across pop/cohorts family mean (mean_rs), already across gardens...



#check factors
str(full)
full$Population <- as.factor(full$Population)
full$Year <- as.factor(full$Year)
full$Dam_ID <- as.factor(full$Dam_ID)

#check factors for full_z
str(full_z)
full_z$Population <- as.factor(full_z$Population)
full_z$Year <- as.factor(full_z$Year)
full_z$Dam_ID <- as.factor(full_z$Dam_ID)




#set priors####
#for z-scored data
class_priors <- c(
  prior(normal(1, 1), class = Intercept), #was (0,1) for zscored fitness, and tried (1,1) for relativized fitness
  prior(normal(0, 2), class = b),# applies to ALL b coefficients
  prior(exponential(1), class = sigma)
)


#1. Model selection####
#across garden FLWR plasticity####

#modlinff: Plasticity effects vary with pop, but only linearly
#vs
#modquadff: Plasticity effects vary with pop, but with quadratic terms


modlinff <- brm(fitness_rel ~ Population * Year + #fitness ~ pop*year
                   doy_pp_z + #flowering plasticity
                   Population:doy_pp_z, #+  #interaction with population
                   #mean_doy_z, #mean flowering time
                   data = full_z,
                 family = gaussian(),                           # Gaussian distribution for response variable
                 chains = 4,                                    # Number of Markov chains
                 cores = 4,                                     # Number of cores for parallel computation
                 iter = 4000,                                   # Number of iterations per chain
                 prior = class_priors,
                 seed = 123
)
summary(modlinff)

modquadff <- brm(fitness_rel ~ Population * Year + #fitness ~ pop*year
              doy_pp_z + doy_pp2_z + #flowering plasticity and quadratic term
              Population:doy_pp_z + Population:doy_pp2_z, #+ #interaction with population
              #mean_doy_z, #mean flowering time
            data = full_z,
            family = gaussian(),                           # Gaussian distribution for response variable
            chains = 4,                                    # Number of Markov chains
            cores = 4,                                     # Number of cores for parallel computation
            iter = 4000,                                   # Number of iterations per chain
            prior = class_priors,
            seed = 123
)
summary(modquadff)
plot(modquadff)

#saveRDS(modlinff, file = "modlinff.rds")
#modlinff <- readRDS("modlinff.rds")
#saveRDS(modquadff, file = "modquadff.rds")
#modquadff <- readRDS("modquadff.rds")



#across garden SLA plasticity####

modlinsla <- brm(fitness_rel ~ Population * Year + #fitness ~ pop*year
                  sla_pp_z + #LDMC pp
                  Population:sla_pp_z, #+  #interaction with population
                #mean_doy_z, #mean flowering time
                data = full_z,
                family = gaussian(),                           # Gaussian distribution for response variable
                chains = 4,                                    # Number of Markov chains
                cores = 4,                                     # Number of cores for parallel computation
                iter = 4000,                                   # Number of iterations per chain
                prior = class_priors,
                seed = 123
)
summary(modlinsla)

modquadsla <- brm(fitness_rel ~ Population * Year + #fitness ~ pop*year
                   sla_pp_z + sla_pp2_z + #flowering plasticity and quadratic term
                   Population:sla_pp_z + Population:sla_pp2_z, #+ #interaction with population
                 #mean_doy_z, #mean flowering time
                 data = full_z,
                 family = gaussian(),                           # Gaussian distribution for response variable
                 chains = 4,                                    # Number of Markov chains
                 cores = 4,                                     # Number of cores for parallel computation
                 iter = 4000,                                   # Number of iterations per chain
                 prior = class_priors,
                 seed = 123
)
summary(modquadsla)
plot(modquadsla)


#saveRDS(modlinsla, file = "modlinsla.rds")
#modlinsla <- readRDS("modlinsla.rds")
#saveRDS(modquadsla, file = "modquadsla.rds")
#modquadsla <- readRDS("modquadsla.rds")





#across garden LDMC plasticity####

modlinldmc <- brm(fitness_rel ~ Population * Year + #fitness ~ pop*year
                   ldmc_pp_z + #LDMC plasticity
                   Population:ldmc_pp_z, #+  #interaction with population
                 #mean_doy_z, #mean flowering time
                 data = full_z,
                 family = gaussian(),                           # Gaussian distribution for response variable
                 chains = 4,                                    # Number of Markov chains
                 cores = 4,                                     # Number of cores for parallel computation
                 iter = 4000,                                   # Number of iterations per chain
                 prior = class_priors,
                 seed = 123
)
summary(modlinldmc)

modquadldmc <- brm(fitness_rel ~ Population * Year + #fitness ~ pop*year
                    ldmc_pp_z + ldmc_pp2_z + #LDMC plasticity and quadratic term
                    Population:ldmc_pp_z + Population:ldmc_pp2_z, #+ #interaction with population
                  #mean_doy_z, #mean flowering time
                  data = full_z,
                  family = gaussian(),                           # Gaussian distribution for response variable
                  chains = 4,                                    # Number of Markov chains
                  cores = 4,                                     # Number of cores for parallel computation
                  iter = 4000,                                   # Number of iterations per chain
                  prior = class_priors,
                  seed = 123
)
summary(modquadldmc)
plot(modquadldmc)


#saveRDS(modlinldmc, file = "modlinldmc.rds")
#modlinldmc <- readRDS("modlinldmc.rds")
#saveRDS(modquadldmc, file = "modquadldmc.rds")
#modquadldmc <- readRDS("modquadldmc.rds")





#2. Supplemental Table S2####
#Loo compare####

#Run loo_compare and convert to tidy data frames
#first flower plasticity
ff_tbl <- loo_compare(loo(modlinff), loo(modquadff)) %>%
  as.data.frame() %>%
  rownames_to_column("model") %>%
  mutate(trait = "First flower date plasticity")
#sla plasticity
sla_tbl <- loo_compare(loo(modlinsla), loo(modquadsla)) %>%
  as.data.frame() %>%
  rownames_to_column("model") %>%
  mutate(trait = "SLA pasticity")
#ldmc plasticity
ldmc_tbl <- loo_compare(loo(modlinldmc), loo(modquadldmc)) %>%
  as.data.frame() %>%
  rownames_to_column("model") %>%
  mutate(trait = "LDMC plasticity")

#Combine into one table- Supplemental Table S2
loo_model_comparison <- bind_rows(ff_tbl, sla_tbl, ldmc_tbl)
loo_model_comparison





#Supplemental Table S3####
#Coefficients table####

#first flower date plasticity
# Convert posterior draws to a dataframe
draws <- as_draws_df(modlinff)

# Define populations
populations <- c("N1", "N2", "C1", "C2", "S1", "S2")

# Identify reference population
pop_columns <- colnames(draws)[grepl("b_Population", colnames(draws))]
pop_reference <- setdiff(populations, gsub("b_Population|:.*", "", pop_columns))
pop_reference  # usually the reference population (e.g., C1)

# Loop over populations and summarize linear effects
plasticity_pop_summary <- map_dfr(populations, function(pop) {
  
  if (pop %in% pop_reference) {
    # Reference population: only main effect
    linear_samples <- draws[["b_doy_pp_z"]]
  } else {
    # Other populations: main effect + population interaction
    linear_samples <- draws[["b_doy_pp_z"]] + draws[[paste0("b_Population", pop, ":doy_pp_z")]]
  }
  
  # Summarize posterior samples
  tibble(
    Population = pop,
    linear_median = median(linear_samples),
    linear_lower  = quantile(linear_samples, 0.025),
    linear_upper  = quantile(linear_samples, 0.975)
  )
})

plasticity_pop_summary <- plasticity_pop_summary %>%
  mutate(trait = "First flower date plasticity")
plasticity_pop_summary 



#sla plasticity

# Convert posterior draws to a dataframe
draws <- as_draws_df(modlinsla)

# Define populations
populations <- c("N1", "N2", "C1", "C2", "S1", "S2")

# Identify reference population
pop_columns <- colnames(draws)[grepl("b_Population", colnames(draws))]
pop_reference <- setdiff(populations, gsub("b_Population|:.*", "", pop_columns))
pop_reference  # usually the reference population (e.g., C1)

# Loop over populations and summarize linear effects
SLAplasticity_pop_summary <- map_dfr(populations, function(pop) {
  
  if (pop %in% pop_reference) {
    # Reference population: only main effect
    linear_samples <- draws[["b_sla_pp_z"]]
  } else {
    # Other populations: main effect + population interaction
    linear_samples <- draws[["b_sla_pp_z"]] + draws[[paste0("b_Population", pop, ":sla_pp_z")]]
  }
  
  # Summarize posterior samples
  tibble(
    Population = pop,
    linear_median = median(linear_samples),
    linear_lower  = quantile(linear_samples, 0.025),
    linear_upper  = quantile(linear_samples, 0.975)
  )
})

SLAplasticity_pop_summary <- SLAplasticity_pop_summary %>%
  mutate(trait = "SLA plasticity")



#ldmc plasticity
# Convert posterior draws to a dataframe
draws <- as_draws_df(modlinldmc)

# Define populations
populations <- c("N1", "N2", "C1", "C2", "S1", "S2")

# Identify reference population
pop_columns <- colnames(draws)[grepl("b_Population", colnames(draws))]
pop_reference <- setdiff(populations, gsub("b_Population|:.*", "", pop_columns))
pop_reference  # usually the reference population (e.g., C1)

# Loop over populations and summarize linear effects
LDMCplasticity_pop_summary <- map_dfr(populations, function(pop) {
  
  if (pop %in% pop_reference) {
    # Reference population: only main effect
    linear_samples <- draws[["b_ldmc_pp_z"]]
  } else {
    # Other populations: main effect + population interaction
    linear_samples <- draws[["b_ldmc_pp_z"]] + draws[[paste0("b_Population", pop, ":ldmc_pp_z")]]
  }
  
  # Summarize posterior samples
  tibble(
    Population = pop,
    linear_median = median(linear_samples),
    linear_lower  = quantile(linear_samples, 0.025),
    linear_upper  = quantile(linear_samples, 0.975)
  )
})

LDMCplasticity_pop_summary <- LDMCplasticity_pop_summary %>%
  mutate(trait = "LDMC plasticity")

#combine tables: Supplemental Table S3
coefficients <- rbind(plasticity_pop_summary, SLAplasticity_pop_summary, LDMCplasticity_pop_summary)
coefficients





#3. Plot figure 3####
#First flower date plasticity####
#plot estimated slopes

full_z <- full_z %>%
  filter(!is.na(Population)) #one row is an NA 

#Get posterior distributions of slopes using emmeans and tidybayes
slope_draws <- modlinff %>% #linear model of first flower date plasticity
  emtrends(~ Population:Year, var = "doy_pp_z") %>%
  gather_emmeans_draws()

#Factor levels
slope_draws <-slope_draws %>%  mutate(Population = factor(Population, levels = c("S2","S1", "C2","C1","N2","N1")),
                                      Year = factor(Year, levels = c("2010", "2017")))

#Calculate probabilities that slopes are different from zero
slope_prob_df <- slope_draws %>%
  group_by(Population, Year) %>%
  summarize(
    mean_slope = mean(.value),
    prob_greater_than_zero = mean(.value > 0),
    prob_less_than_zero = mean(.value < 0),
    prob_diff_from_zero = max(prob_greater_than_zero, prob_less_than_zero)
  )
slope_prob_df

#The model doesnt allow slopes across cohorts to be different
#Supplemental Figure S7####
# Plot the slopes with uncertainty using tidybayes
ggplot() +
  # 2010 (darker)
  stat_halfeye(data = filter(slope_draws, Year == 2010),
               aes(x = .value, y = Population, fill = Population
                   #, alpha = factor(Year)
                   )) +
  
  # 2017 (lighter)
  #stat_halfeye(data = filter(slope_draws, Year == 2017),
  #            aes(x = .value, y = Population, fill = Population, alpha = factor(Year))) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  
  # Text for 2010
  geom_text(data = filter(slope_prob_df, Year == 2010), color = "black",
            aes(x = 0.15, y = Population,
                label = sprintf("PD(≠0) = %.2f", prob_diff_from_zero)),
            hjust = 0, vjust = -2.5, size = 3.5) +
  
  # Text for 2017
  #geom_text(data = filter(slope_prob_df, Year == 2017), color = "black",
  #          aes(x = 0.2, y = Population,
  #              label = sprintf("p(≠0) = %.3f", prob_diff_from_zero)),
  #          hjust = 0, vjust = -1, size = 3.5, fontface = "bold") +
  
  labs(#title = "Effect of plasticity on fitness",
       #subtitle = "2010 (light) and 2017 (dark) overlaid for each Population",
       x = "Slope coefficient (effect of adjusted first flower plasticity)",
       y = "Source Population",
       #alpha = "Year"
       ) +  # this label now applies
  scale_fill_manual(values = cols) +
  #scale_alpha_manual(values = c("2010" = 0.35, "2017" = 0.75))+
  xlim(-0.15, 0.3)+ 
  guides(fill = guide_legend(reverse = TRUE))




#Plot slopes with restricted ranges
#Find observed range of doy_pp_z for each Population-Year combo
restricted_ranges <- full_z %>%
  group_by(Population, Year) %>%
  summarize(
    min_doy = min(doy_pp_z, na.rm = TRUE),
    max_doy = max(doy_pp_z, na.rm = TRUE),
    .groups = "drop"
  )

#Build grid using those observed ranges
plasticity_grid_restricted <- do.call(rbind, lapply(1:nrow(restricted_ranges), function(i) {
  r <- restricted_ranges[i, ]
  expand.grid(
    doy_pp_z = seq(r$min_doy, r$max_doy, length.out = 100),
    Population = r$Population,
    Year = r$Year,
    mean_doy_z = mean(full_z$mean_doy_z, na.rm = TRUE)
  )
})) #%>%
  #mutate(doy_pp2_z = doy_pp_z^2)#quad terms if included in mod

conditional_effects <- fitted(
  modlinff, #flowering time pp
  newdata = plasticity_grid_restricted,
  re_formula = NA,
  summary = FALSE
)

drawslinff <- as_draws_df(modlinff) #convert posterior draws to df

conditional_df <- tibble(
  doy_pp_z = rep(plasticity_grid_restricted$doy_pp_z, each = nrow(drawslinff)), 
  Population = rep(plasticity_grid_restricted$Population, each = nrow(drawslinff)), 
  Year = rep(plasticity_grid_restricted$Year, each = nrow(drawslinff)),
  .draw = rep(1:nrow(drawslinff), times = nrow(plasticity_grid_restricted)),
  .value = as.vector(conditional_effects)
)

#slope_prob_df
slope_prob_df <- slope_prob_df %>%
  mutate(pop_year=paste(Population,"-",Year, sep=''))

conditional_df1 <- conditional_df %>%
  mutate(pop_year=paste(Population,"-",Year, sep=''))

conditional_df1 <- conditional_df %>%
  mutate(pop_year = paste(Population, "-", Year, sep = '')) %>%
  left_join(slope_prob_df %>% dplyr::select(pop_year, prob_diff_from_zero, Population), by = "pop_year") %>%
  dplyr::select(-Population.y) %>%  # Keep Population.x, drop duplicate
  rename(Population = Population.x) %>%
  mutate(sig = ifelse(prob_diff_from_zero >= 0.95, "Significant", "Not Significant")) #add in significant/not significant at 0.95

#Factor levels
full_z <-full_z %>%  mutate(Population = factor(Population, levels = c("N1","N2", "C2","C1","S2","S1")),
                                      Year = factor(Year, levels = c("2010", "2017")))

#Fig. 3 first flower date plasticity plot####
#with raw data points
ffplot <- ggplot() +
  # Points
  geom_point(data = full_z, 
             aes(x = doy_pp_z, y = fitness_rel, 
                 fill = Population, 
                 alpha = factor(Year)),
             color = "black",  
             shape = 21,
             size = 1.5) +
  
  # Lines + ribbons
  stat_lineribbon(#data = conditional_df1,
    data = dplyr::filter(conditional_df1, sig =="Significant"), 
                  aes(x = doy_pp_z, y = .value, 
                      fill = Population, 
                      alpha = factor(Year), #could do color within aes if you want diff colors per line
                      group = interaction(Population, Year), 
                      linetype = sig),
                  color = "black", #color of predicted line
                  .width = 0.50, #width of credible interval
                  linewidth = 1) +  
  
  facet_wrap(~ Population, ncol =1
             #, labeller = labeller(Population = pop_labels)
             )+
  
  labs(#title = "Estimated relationship between plasticity and fitness",
       #subtitle = "With raw data points, 50% credible intervals",
       x = "Adjusted First Flower Plasticity",
       y = "Relative Fitness (No. Reproductive Structures)", 
       fill = "Population",
       alpha = "Year", 
       linetype = "Significance") +
  
  scale_fill_manual(values = cols) +
  scale_alpha_manual(values = c("2010" = 0.35, "2017" = 0.75))+
  scale_linetype_manual(values = c("Significant" = "solid", "Not Significant" = "11"))+ #line type
  
 # guides(
  #  fill = "none", #remove population legend since colors match, and there is a header for each box
    #fill = guide_legend(order = 1),       # Population legend first
   # alpha = guide_legend(order = 1),      # Year legend second
   # linetype = guide_legend(order = 2)    # Significance legend last
  #)+
  theme_bw()+
  coord_cartesian(ylim = c(0,3.5))+ #zoom in just a little
  theme(legend.position = "none")

ffplot


#SLA plasticity ####

drawslinsla <- as_draws_df(modlinsla) #convert posterior draws to df 

#Get posterior distributions of slopes using emmeans and tidybayes
slope_draws <- modlinsla %>% #modJ
  emtrends(~ Population:Year, var = "sla_pp_z") %>%
  gather_emmeans_draws()

#Reorder levels
slope_draws <-slope_draws %>%  mutate(Population = factor(Population, levels = c("S2","S1", "C2","C1","N2","N1")),
                                      Year = factor(Year, levels = c("2010", "2017")))

#Calculate probabilities that slopes are different from zero
slope_prob_df <- slope_draws %>%
  group_by(Population, Year) %>%
  summarize(
    mean_slope = mean(.value),
    prob_greater_than_zero = mean(.value > 0),
    prob_less_than_zero = mean(.value < 0),
    prob_diff_from_zero = max(prob_greater_than_zero, prob_less_than_zero)
  )
slope_prob_df

#model doesnt allow slopes across cohorts to be different
#Supplemental Figure S8####
#Plot the slopes with uncertainty using tidybayes
ggplot() +
  # 2010 (darker)
  stat_halfeye(data = filter(slope_draws, Year == 2010),
               aes(x = .value, y = Population, fill = Population
                   #, alpha = factor(Year)
               )) +
  
  # 2017 (lighter)
  #stat_halfeye(data = filter(slope_draws, Year == 2017),
  #            aes(x = .value, y = Population, fill = Population, alpha = factor(Year))) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  
  # Text for 2010
  geom_text(data = filter(slope_prob_df, Year == 2010), color = "black",
            aes(x = 0.1, y = Population,
                label = sprintf("PD(≠0) = %.2f", prob_diff_from_zero)),
            hjust = 0, vjust = -2.5, size = 3.5) +
  
  # Text for 2017
  #geom_text(data = filter(slope_prob_df, Year == 2017), color = "black",
  #          aes(x = 0.2, y = Population,
  #              label = sprintf("p(≠0) = %.3f", prob_diff_from_zero)),
  #          hjust = 0, vjust = -1, size = 3.5, fontface = "bold") +
  
  labs(#title = "Effect of plasticity on fitness",
    #subtitle = "2010 (light) and 2017 (dark) overlaid for each Population",
    x = "Slope coefficient (effect of specific leaf area plasticity)",
    y = "Source Population",
    #alpha = "Year"
  ) +  # this label now applies
  scale_fill_manual(values = cols) +
  #scale_alpha_manual(values = c("2010" = 0.35, "2017" = 0.75))+
  xlim(-0.2, 0.25)+ 
  guides(fill = guide_legend(reverse = TRUE))



#Plot slopes with restricted ranges
#Find observed range of sla_pp_z for each Population-Year combo
restricted_ranges <- full_z %>%
  group_by(Population, Year) %>%
  summarize(
    min_sla = min(sla_pp_z, na.rm = TRUE),
    max_sla = max(sla_pp_z, na.rm = TRUE),
    .groups = "drop"
  )

#Build grid using those observed ranges
plasticity_grid_restricted <- do.call(rbind, lapply(1:nrow(restricted_ranges), function(i) {
  r <- restricted_ranges[i, ]
  expand.grid(
    sla_pp_z = seq(r$min_sla, r$max_sla, length.out = 100),
    Population = r$Population,
    Year = r$Year,
    mean_sla_z = mean(full_z$mean_sla_z, na.rm = TRUE)
  )
})) %>%
  mutate(sla_pp2_z = sla_pp_z^2)

conditional_effects <- fitted(
  modlinsla, #linear model for sla plasticity
  newdata = plasticity_grid_restricted,
  re_formula = NA,
  summary = FALSE
)

conditional_df <- tibble(
  sla_pp_z = rep(plasticity_grid_restricted$sla_pp_z, each = nrow(drawslinsla)), 
  Population = rep(plasticity_grid_restricted$Population, each = nrow(drawslinsla)), 
  Year = rep(plasticity_grid_restricted$Year, each = nrow(drawslinsla)),
  .draw = rep(1:nrow(drawslinsla), times = nrow(plasticity_grid_restricted)),
  .value = as.vector(conditional_effects)
)

#slope_prob_df
slope_prob_df <- slope_prob_df %>%
  mutate(pop_year=paste(Population,"-",Year, sep=''))

conditional_df1 <- conditional_df %>%
  mutate(pop_year=paste(Population,"-",Year, sep=''))


conditional_df1 <- conditional_df %>%
  mutate(pop_year = paste(Population, "-", Year, sep = '')) %>%
  left_join(slope_prob_df %>% dplyr::select(pop_year, prob_diff_from_zero, Population), by = "pop_year") %>%
  dplyr::select(-Population.y) %>%  # Keep Population.x, drop duplicate
  rename(Population = Population.x) %>%
  mutate(sig = ifelse(prob_diff_from_zero >= 0.95, "Significant", "Not Significant"))


#Fig. 3 SLA plasticity plot####
#with raw data points
slaplot <- 
  ggplot() +
  # Points
  geom_point(data = full_z, 
             aes(x = sla_pp_z, y = fitness_rel, 
                 fill = Population, 
                 alpha = factor(Year)),
             color = "black",  
             shape = 21,
             size = 1.5) +
  
  # Lines + ribbons
  stat_lineribbon(#data = conditional_df1,
    data = dplyr::filter(conditional_df1, sig =="Significant"), 
    aes(x = sla_pp_z, y = .value, 
        fill = Population, 
        alpha = factor(Year), #could do color within aes if you want diff colors per line
        group = interaction(Population, Year), 
        linetype = sig),
    color = "black", #color of predicted line
    .width = 0.50, #width of credible interval
    linewidth = 1) +  
  
  facet_wrap(~ Population, ncol =1
             #, labeller = labeller(Population = pop_labels)
  )+
  
  labs(#title = "Estimated relationship between plasticity and fitness",
    #subtitle = "With raw data points, 50% credible intervals",
    x = "Specific Leaf Area Plasticity",
    y = "Relative Fitness (No. Reproductive Structures)", 
    fill = "Population",
    alpha = "Year", 
    linetype = "Significance") +
  
  scale_fill_manual(values = cols) +
  scale_alpha_manual(values = c("2010" = 0.35, "2017" = 0.75))+
  scale_linetype_manual(values = c("Significant" = "solid", "Not Significant" = "11"))+ #line type
  
 #guides(
  #  fill = "none", #remove population legend since colors match, and there is a header for each box
    #fill = guide_legend(order = 1),       # Population legend first
  #  alpha = guide_legend(order = 1),      # Year legend second
   # linetype = guide_legend(order = 2)    # Significance legend last
  #)+
  theme_bw()+
  coord_cartesian(ylim = c(0,3.5))+ #zoom in just a little
  theme(legend.position = "none")

slaplot


#LDMC plasticity####

drawslinldmc <- as_draws_df(modlinldmc) #convert posterior draws to df 

#Get posterior distributions of slopes using emmeans and tidybayes
slope_draws <- modlinldmc %>% 
  emtrends(~ Population:Year, var = "ldmc_pp_z") %>%
  gather_emmeans_draws()

#Reorder levels
slope_draws <-slope_draws %>%  mutate(Population = factor(Population, levels = c("S2","S1", "C2","C1","N2","N1")),
                                      Year = factor(Year, levels = c("2010", "2017")))

#Calculate probabilities that slopes are different from zero
slope_prob_df <- slope_draws %>%
  group_by(Population, Year) %>%
  summarize(
    mean_slope = mean(.value),
    prob_greater_than_zero = mean(.value > 0),
    prob_less_than_zero = mean(.value < 0),
    prob_diff_from_zero = max(prob_greater_than_zero, prob_less_than_zero)
  )
slope_prob_df

#model doesnt allow slopes across cohorts to be different
#Supplemental Figure S9####
#Plot the slopes with uncertainty using tidybayes
ggplot() +
  # 2010 (darker)
  stat_halfeye(data = filter(slope_draws, Year == 2010),
               aes(x = .value, y = Population, fill = Population
                   #, alpha = factor(Year)
               )) +
  
  # 2017 (lighter)
  #stat_halfeye(data = filter(slope_draws, Year == 2017),
  #            aes(x = .value, y = Population, fill = Population, alpha = factor(Year))) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  
  # Text for 2010
  geom_text(data = filter(slope_prob_df, Year == 2010), color = "black",
            aes(x = 0.1, y = Population,
                label = sprintf("PD(≠0) = %.2f", prob_diff_from_zero)),
            hjust = 0, vjust = -2.5, size = 3.5) +
  
  # Text for 2017
  #geom_text(data = filter(slope_prob_df, Year == 2017), color = "black",
  #          aes(x = 0.2, y = Population,
  #              label = sprintf("p(≠0) = %.3f", prob_diff_from_zero)),
  #          hjust = 0, vjust = -1, size = 3.5, fontface = "bold") +
  
  labs(#title = "Effect of plasticity on fitness",
    #subtitle = "2010 (light) and 2017 (dark) overlaid for each Population",
    x = "Slope coefficient (effect of leaf dry matter content plasticity)",
    y = "Source Population",
    #alpha = "Year"
  ) +  # this label now applies
  scale_fill_manual(values = cols) +
  #scale_alpha_manual(values = c("2010" = 0.35, "2017" = 0.75))+
  xlim(-0.2, 0.25)+ 
  guides(fill = guide_legend(reverse = TRUE))


#Plot slopes with restricted ranges
#Find observed range of ldmc_pp_z for each Population-Year combo
restricted_ranges <- full_z %>%
  group_by(Population, Year) %>%
  summarize(
    min_ldmc = min(ldmc_pp_z, na.rm = TRUE),
    max_ldmc = max(ldmc_pp_z, na.rm = TRUE),
    .groups = "drop"
  )

#Build grid using those observed ranges
plasticity_grid_restricted <- do.call(rbind, lapply(1:nrow(restricted_ranges), function(i) {
  r <- restricted_ranges[i, ]
  expand.grid(
    ldmc_pp_z = seq(r$min_ldmc, r$max_ldmc, length.out = 100),
    Population = r$Population,
    Year = r$Year,
    mean_ldmc_z = mean(full_z$mean_ldmc_z, na.rm = TRUE)
  )
}))# %>%
  #mutate(ldmc_pp2_z = ldmc_pp_z^2)

conditional_effects <- fitted(
 modlinldmc, 
  newdata = plasticity_grid_restricted,
  re_formula = NA,
  summary = FALSE
)

conditional_df <- tibble(
  ldmc_pp_z = rep(plasticity_grid_restricted$ldmc_pp_z, each = nrow(drawslinldmc)), 
  Population = rep(plasticity_grid_restricted$Population, each = nrow(drawslinldmc)), 
  Year = rep(plasticity_grid_restricted$Year, each = nrow(drawslinldmc)),
  .draw = rep(1:nrow(drawslinldmc), times = nrow(plasticity_grid_restricted)),
  .value = as.vector(conditional_effects)
)

#slope_prob_df
slope_prob_df <- slope_prob_df %>%
  mutate(pop_year=paste(Population,"-",Year, sep=''))

conditional_df1 <- conditional_df %>%
  mutate(pop_year=paste(Population,"-",Year, sep=''))

conditional_df1 <- conditional_df %>%
  mutate(pop_year = paste(Population, "-", Year, sep = '')) %>%
  left_join(slope_prob_df %>% dplyr::select(pop_year, prob_diff_from_zero, Population), by = "pop_year") %>%
  dplyr::select(-Population.y) %>%  # Keep Population.x, drop duplicate
  rename(Population = Population.x) %>%
  mutate(sig = ifelse(prob_diff_from_zero >= 0.95, "Significant", "Not Significant"))

#Fig. 3 LDMC plasticity plot####
ldmcplot <- 
  ggplot() +
  # Points
  geom_point(data = full_z, 
             aes(x = ldmc_pp_z, y = fitness_rel, 
                 fill = Population, 
                 alpha = factor(Year)),
             color = "black",  
             shape = 21,
             size = 1.5) +
  
  # Lines + ribbons
  stat_lineribbon(#data = conditional_df1,
    data = dplyr::filter(conditional_df1, sig =="Significant"), 
    aes(x = ldmc_pp_z, y = .value, 
        fill = Population, 
        alpha = factor(Year), #could do color within aes if you want diff colors per line
        group = interaction(Population, Year), 
        linetype = sig),
    color = "black", #color of predicted line
    .width = 0.50, #width of credible interval
    linewidth = 1) +  
  
  facet_wrap(~ Population, ncol =1
             #, labeller = labeller(Population = pop_labels)
  )+
  
  labs(#title = "Estimated relationship between plasticity and fitness",
    #subtitle = "With raw data points, 50% credible intervals",
    x = "Leaf Dry Matter Content Plasticity",
    y = "Relative Fitness (No. Reproductive Structures)", 
    fill = "Population",
    alpha = "Cohort", 
    linetype = "Significance") +
  
  scale_fill_manual(values = cols) +
  scale_alpha_manual(values = c("2010" = 0.35, "2017" = 0.75))+
  scale_linetype_manual(values = c("Significant" = "solid", "Not Significant" = "11"))+ #line type
  
  guides(
    #fill = "none", #remove population legend since colors match, and there is a header for each box
    fill = guide_legend(order = 1),       # Population legend first
    alpha = guide_legend(order = 2),      # Year legend second
    #linetype = guide_legend(order = 2)    # Significance legend last
  )+
  theme_bw()+
  coord_cartesian(ylim = c(0,3.5))#+ #zoom in just a little
#theme(legend.position = "none")
ldmcplot



#join fitness plots####
library(patchwork)
#ffplot <- ffplot + theme(legend.position = "none") #remove legend
#slaplot <- slaplot + theme(legend.position = "none")

ffplot_A = ffplot +
  ggtitle('A') + theme(plot.title=element_text(hjust=0,size=20)) 
slaplot_B = slaplot +
  ggtitle('B') + theme(plot.title=element_text(hjust=0,size=20)) 
ldmcplot_C = ldmcplot +
  ggtitle('C') + theme(plot.title=element_text(hjust=0,size=20)) 




fitplot <- ffplot_A + slaplot_B + ldmcplot_C #side by side
fitplot

#END ACROSS GARDEN####

