#Author: Emily Cook
#Project: PERSIST plasticity
#Purpose: 
  #(Q1) Does the magnitude of plasticity differ among populations?
  #(Q2) Did the magnitude of plasticity evolve following a historic drought? 
#Using brms, estimate uncertainty for differences among populations (Q1), cohorts (Q2) in trait plasticity
#Traits: first flower date plasticity, SLA plasticity, and LDMC plasticity

#Workflow:
#Get Posterior predictions of the mean response for each group (cohort, population, cohort x population)
#Compute marginal predictions by averaging over pop/cohort or both

#load packages
library(brms)
library(tidybayes)
library(tictoc)
library(tidyverse) 
library(ggridges)
library(bayesplot)
#library(stat_pointinterval)


#load dataframes from 2_calculate_plasticity.R:
doy_pp_mn <- read.csv("doy_pp.csv") #read
sla_pp_mn <- read.csv("sla_pp.csv") #read
ldmc_pp_mn <- read.csv("ldmc_pp.csv") #read


#First flower date plasticity brms model####
#check factors
str(doy_pp_mn)
doy_pp_mn$Population <- as.factor(doy_pp_mn$Population)
doy_pp_mn$Year <- as.factor(doy_pp_mn$Year)

#run brm model:
mod11 <- brm(
  doy_plasticity ~ Population * Year # Fixed effects with interactions
  ,                                 
  data = doy_pp_mn,
  family = gaussian(),                           # Gaussian distribution for response variable
  chains = 4,                                    # Number of Markov chains
  cores = 4,                                     # Number of cores for parallel computation
  iter = 4000,                                   # Number of iterations per chain
)

summary(mod11) 
plot(mod11)
prior_summary(mod11)

#saveRDS(mod11, file = "mod11brms.rds")
#mod11 <- readRDS("mod11brms.rds")

#Create empty dataframe
cohort_data <- expand.grid( 
  Population = unique(doy_pp_mn$Population)
  ,Year=unique(doy_pp_mn$Year)
)

#Make posterior predictions of the mean response for each group
post_preds <-as.data.frame(posterior_epred(mod11, newdata = cohort_data, allow_new_levels=TRUE, ndraws=1000) )
colnames(post_preds) <- paste(cohort_data$Population, cohort_data$Year, sep="_") #this should assign based on order of new_data$garden

head(post_preds);dim(post_preds)

# Convert to long format
post_preds_long <- post_preds %>%
  mutate(draw = row_number()) %>%  # Add index for posterior draws
  pivot_longer(cols = -draw, names_to = "group", values_to = "prediction") %>%
  separate(group, into = c("Population", "Year"), sep = "_", remove = FALSE)
head(post_preds_long); dim(post_preds_long)

#Cohort effect (Q2)####
#Compute marginal predictions by averaging over both population
post_preds_marginal <- post_preds_long %>%
  group_by(draw, Year) %>%
  summarise(marginal_prediction = mean(prediction), .groups = "drop") 
head(post_preds_marginal); dim(post_preds_marginal)

difs_long <- post_preds_marginal %>%
  group_by(draw) %>% 
  summarise(diff_prediction = marginal_prediction[Year == 2010] - marginal_prediction[Year == 2017], .groups = "drop")

head(difs_long)

#Use the posterior distribution of each pairwise comparison to calculate the probability that cohorts were different (posterior probability of difference, PD)
cohort_probs <- difs_long %>%
  summarise(
    prob_greater_than_zero = mean(diff_prediction > 0),
    prob_less_than_zero    = mean(diff_prediction < 0),
    probability            = max(prob_greater_than_zero,
                                 prob_less_than_zero)
  )

difs_long <- difs_long %>%
  mutate(probability = cohort_probs$probability)

#plot overall cohort effect####
ggplot(difs_long, aes(x = diff_prediction)) +
  geom_density(fill = "gray",color = "darkgray", scale = 1) +  # Density plot
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Line at 0
  labs(
   # title = "Posterior Distribution of Differences (2010 - 2017)",
    x = "Difference in predicted posteriors (2010 - 2017)",
    y = "Density"
  ) +
  theme_minimal(base_size = 20)+
  labs(x="Difference in Predicted Posteriors")+
  # Add the probability as text only once per facet
  geom_text(data = difs_long, aes(x = 0, y = 1.5, label = paste0("PD(≠0) = ", round(probability, 2))), 
            hjust = 0, vjust = 1, size = 6, color = "black")

#Population effect (Q1)####
#Compute marginal predictions
post_preds_marginal <- post_preds_long %>%
  group_by(draw, Population) %>% 
  summarise(marginal_prediction = mean(prediction), .groups = "drop") 
head(post_preds_marginal); dim(post_preds_marginal)

#Take subtractions using a "wider" df
population_difs <- post_preds_marginal %>%
  pivot_wider(names_from = Population, values_from = marginal_prediction) %>%
  mutate(
    N1_C1 = N1 - C1, 
    N1_C2 = N1 - C2, 
    N2_C1 = N2 - C1, 
    N2_C2 = N2 - C2, 
    N1_S1 = N1 - S1, 
    N1_S2 = N1 - S2,
    N2_S1 = N2 - S1,
    N2_S2 = N2 - S2,
    C1_S1 = C1 - S1,
    C1_S2 = C1 - S2,
    C2_S1 = C2 - S1,
    C2_S2 = C2 - S2,
    N1_N2 = N1 - N2, 
    C1_C2 = C1 - C2, 
    S1_S2 = S1 - S2
  )

head(population_difs)

#pivot longer
population_difs_long <- population_difs %>%
  pivot_longer(
    cols = c(N1_C1, N1_C2, N2_C1, N2_C2,
             N1_S1, N1_S2, N2_S1, N2_S2,
             C1_S1, C1_S2, C2_S1, C2_S2,
             N1_N2, C1_C2, S1_S2),  # Columns to pivot
    names_to = "difference",   # New column for difference labels
    values_to = "pred_diff"    # New column for values
  )

#Posterior mean difference
#Use the posterior distribution of each pairwise comparison to calculate the probability that populations were different (posterior probability of difference, PD)
diffs <- population_difs_long %>%
  group_by(difference)%>%
  summarise(
    prob_greater_than_zero = mean(pred_diff > 0),
    prob_less_than_zero    = mean(pred_diff < 0),
    probability            = max(prob_greater_than_zero,
                                 prob_less_than_zero)
  )

population_difs_long <- population_difs_long %>%
  left_join(diffs, by = "difference")

#Supplemental Figure S1####
ggplot(population_difs_long, aes(x = pred_diff, y = difference, fill = difference)) +
  geom_density_ridges(fill = "gray", color = "darkgray", scale = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  
  theme_minimal(base_size = 20) +
  theme(legend.position = "none") +
  labs(x="Difference in Predicted Posteriors", y="Group Difference Density")+
  scale_y_discrete(limits = rev(c(  # Use `rev()` to flip order if needed
    "N1_C1", "N1_C2", "N2_C1", "N2_C2",
    "N1_S1", "N1_S2", "N2_S1", "N2_S2",
    "C1_S1", "C1_S2", "C2_S1", "C2_S2",
    "N1_N2", "C1_C2", "S1_S2"
  ))) + 
  geom_text(data = population_difs_long, aes(x = 0.18, y = difference, 
                                             label = paste0("PD(≠0) = ", round(probability, 2))),
            hjust = 0.7, vjust = 0.05, size = 6, color = "black", inherit.aes = FALSE)



#Pop-cohort effect (Q2)####
#Compute marginal predictions
post_preds_marginal <- post_preds_long %>%
  group_by(draw, Population, Year) %>%
  summarise(marginal_prediction = mean(prediction), .groups = "drop") 
head(post_preds_marginal); dim(post_preds_marginal)

#Specify the comparisons (i.e. is N1 2010 greater than N2 2010?)
comparisons <- list(
  #2010 comparisons within and across populations
  c("N1_2010", "N2_2010"), c("N1_2010", "C1_2010"), c("N1_2010", "C2_2010"),
  c("N1_2010", "S1_2010"), c("N1_2010", "S2_2010"), c("N2_2010", "C1_2010"),
  c("N2_2010", "C2_2010"), c("N2_2010", "S1_2010"), c("N2_2010", "S2_2010"),
  c("C1_2010", "C2_2010"), c("C1_2010", "S1_2010"), c("C1_2010", "S2_2010"),
  c("C2_2010", "S1_2010"), c("C2_2010", "S2_2010"), c("S1_2010", "S2_2010"),
  
  #2017 comparisons within and across populations
  c("N1_2017", "N2_2017"), c("N1_2017", "C1_2017"), c("N1_2017", "C2_2017"),
  c("N1_2017", "S1_2017"), c("N1_2017", "S2_2017"), c("N2_2017", "C1_2017"),
  c("N2_2017", "C2_2017"), c("N2_2017", "S1_2017"), c("N2_2017", "S2_2017"),
  c("C1_2017", "C2_2017"), c("C1_2017", "S1_2017"), c("C1_2017", "S2_2017"),
  c("C2_2017", "S1_2017"), c("C2_2017", "S2_2017"), c("S1_2017", "S2_2017"),
  
  #Cross-year comparisons
  c("N1_2010", "N1_2017"), c("N2_2010", "N2_2017"), c("C1_2010", "C1_2017"),
  c("C2_2010", "C2_2017"), c("S1_2010", "S1_2017"), c("S2_2010", "S2_2017")
)

#Pivot wider to create separate columns for each Population-Year group
post_preds_wide <- post_preds_marginal %>%
  pivot_wider(names_from = c(Population, Year), values_from = marginal_prediction)
head(post_preds_wide); dim(post_preds_wide)

#Compute row-wise differences for specified comparisons
df_differences <- post_preds_wide %>%
  rowwise() %>%
  mutate(
    !!!setNames(
      lapply(comparisons, function(cols) .[[cols[1]]] - .[[cols[2]]]),
      sapply(comparisons, function(cols) paste(cols, collapse = "-"))
    )
  ) %>%
  ungroup()
head(df_differences)

#probabilities - only keep columns 14- 49
df_differences <- df_differences[, 14:49]

difs_long <- df_differences %>%
  #select(14:49) %>%
  pivot_longer(cols = everything(), names_to = "var", values_to = "dif")

#Posterior mean difference
#Use the posterior distribution of each pairwise comparison to calculate the probability that populations were different (posterior probability of difference, PD)
diffs <- difs_long %>% group_by(var) %>%
  summarise(
    prob_greater_than_zero = mean(dif > 0),
    prob_less_than_zero    = mean(dif < 0),
    probability            = max(prob_greater_than_zero,
                                 prob_less_than_zero)
  )
diffs

difs_long <- difs_long %>%
  left_join(diffs, by = "var") #add "probability" back in to raw data for graph

#Supplemental Figure S2####
ggplot(difs_long, aes(x = dif, y = var, fill = var)) +
  geom_density_ridges(fill = "gray", color = "darkgray", scale = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  
  theme_minimal(base_size = 18) +
  theme(legend.position = "none") +
  labs(x="Difference in Predicted Posteriors", y="Group Difference Density")+
  scale_y_discrete(limits = rev(c(  # Use `rev()` to flip order if needed
    "N1_2010-N2_2010", "C1_2010-C2_2010", "S1_2010-S2_2010",
    "N1_2017-N2_2017", "C1_2017-C2_2017", "S1_2017-S2_2017",
    "N1_2010-C1_2010", "N1_2010-C2_2010", "N2_2010-C1_2010", "N2_2010-C2_2010",
    "N1_2010-S1_2010", "N1_2010-S2_2010", "N2_2010-S1_2010", "N2_2010-S2_2010",
    "C1_2010-S1_2010", "C1_2010-S2_2010", "C2_2010-S1_2010", "C2_2010-S2_2010",
    "N1_2017-C1_2017", "N1_2017-C2_2017", "N2_2017-C1_2017", "N2_2017-C2_2017",
    "N1_2017-S1_2017", "N1_2017-S2_2017", "N2_2017-S1_2017", "N2_2017-S2_2017",
    "C1_2017-S1_2017", "C1_2017-S2_2017", "C2_2017-S1_2017", "C2_2017-S2_2017",
    "N1_2010-N1_2017", "N2_2010-N2_2017", "C1_2010-C1_2017", "C2_2010-C2_2017",
    "S1_2010-S1_2017", "S2_2010-S2_2017"
  ))) + 
  #geom_text(data = difs_long, aes(x = 3, y = as.numeric(var) + 1, 
  #                                label = paste0("P(pred_diff < 0) = ", round(probability, 3))),
  #         hjust = 0, vjust = 1, size = 3, color = "black", inherit.aes = FALSE)
  geom_text(data = difs_long, aes(x = 0.2, y = var, 
                                  label = paste0("PD(≠0) = ", round(probability, 2))),
            hjust = 0, vjust = 0.05, size = 4, color = "black", inherit.aes = FALSE)+
  coord_cartesian(clip = "off")



#SLA plasticity brms model####
#check factors
sla_pp_mn$Population <- as.factor(sla_pp_mn$Population)
sla_pp_mn$Year <- as.factor(sla_pp_mn$Year)

mod12 <- brm(
  sla_plasticity ~ Population * Year # Fixed effects with interactions
  ,                                 
  data = sla_pp_mn,
  family = gaussian(),                           # Gaussian distribution for response variable
  chains = 4,                                    # Number of Markov chains
  cores = 4,                                     # Number of cores for parallel computation
  iter = 4000,                                   # Number of iterations per chain
)

summary(mod12) 
plot(mod12)

#saveRDS(mod12, file = "mod12brms.rds")
#mod12 <- readRDS("mod12brms.rds")

#Create empty dataframe
cohort_data_sla <- expand.grid( 
  Population = unique(doy_pp_mn$Population)
  ,Year=unique(doy_pp_mn$Year)
)

#Posterior predictions of the mean response for each group
post_preds_sla <-as.data.frame(posterior_epred(mod12, newdata = cohort_data_sla, allow_new_levels=TRUE, ndraws=1000) )
colnames(post_preds_sla) <- paste(cohort_data_sla$Population, cohort_data_sla$Year, sep="_")
head(post_preds_sla);dim(post_preds_sla)

#Convert to long format
post_preds_long_sla <- post_preds_sla %>%
  mutate(draw = row_number()) %>%  # Add index for posterior draws
  pivot_longer(cols = -draw, names_to = "group", values_to = "prediction") %>%
  separate(group, into = c("Population", "Year"), sep = "_", remove = FALSE)
head(post_preds_long_sla); dim(post_preds_long_sla)

#Cohort effect (Q2)####
#Compute marginal predictions 
post_preds_marginal <- post_preds_long_sla %>%
  group_by(draw, Year) %>%  
  summarise(marginal_prediction = mean(prediction), .groups = "drop") 
head(post_preds_marginal); dim(post_preds_marginal)

difs_long <- post_preds_marginal %>%
  group_by(draw) %>% 
  summarise(diff_prediction = marginal_prediction[Year == 2010] - marginal_prediction[Year == 2017], .groups = "drop")

head(difs_long)

#probabilities

#Use the posterior distribution of each pairwise comparison to calculate the probability that cohorts were different (posterior probability of difference, PD)
cohort_probs <- difs_long %>%
  summarise(
    prob_greater_than_zero = mean(diff_prediction > 0),
    prob_less_than_zero    = mean(diff_prediction < 0),
    probability            = max(prob_greater_than_zero,
                                 prob_less_than_zero)
  )

difs_long <- difs_long %>%
  mutate(probability = cohort_probs$probability)

#plot overall cohort effect####
ggplot(difs_long, aes(x = diff_prediction)) +
  geom_density(fill = "gray",color = "darkgray", scale = 1) +  # Density plot
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Line at 0
  labs(
   #title = "Posterior Distribution of Differences (2010 - 2017)",
    x = "Difference (2010 - 2017)",
    y = "Density"
  ) +
  theme_minimal(base_size = 20)+
  labs(x="Difference in Predicted Posteriors")+
  # Add the probability as text only once per facet
  geom_text(data = difs_long, aes(x = 0, y = 1.5, label = paste0("PD(≠0) = ", round(probability, 2))), 
            hjust = 0, vjust = 1, size = 6, color = "black")


#Population effect (Q1)####
#Compute marginal predictions
post_preds_marginal <- post_preds_long_sla %>%
  group_by(draw, Population) %>% 
  summarise(marginal_prediction = mean(prediction), .groups = "drop") 
head(post_preds_marginal); dim(post_preds_marginal)

#Take subtractions using a "wider" df
population_difs <- post_preds_marginal %>%
  pivot_wider(names_from = Population, values_from = marginal_prediction) %>%
  mutate(
    N1_C1 = N1 - C1, 
    N1_C2 = N1 - C2, 
    N2_C1 = N2 - C1, 
    N2_C2 = N2 - C2, 
    N1_S1 = N1 - S1, 
    N1_S2 = N1 - S2,
    N2_S1 = N2 - S1,
    N2_S2 = N2 - S2,
    C1_S1 = C1 - S1,
    C1_S2 = C1 - S2,
    C2_S1 = C2 - S1,
    C2_S2 = C2 - S2,
    N1_N2 = N1 - N2, 
    C1_C2 = C1 - C2, 
    S1_S2 = S1 - S2
  )
head(population_difs)

#pivot longer
population_difs_long <- population_difs %>%
  pivot_longer(
    cols = c(N1_C1, N1_C2, N2_C1, N2_C2,
             N1_S1, N1_S2, N2_S1, N2_S2,
             C1_S1, C1_S2, C2_S1, C2_S2,
             N1_N2, C1_C2, S1_S2),  # Columns to pivot
    names_to = "difference",   # New column for difference labels
    values_to = "pred_diff"    # New column for values
  )

#Posterior mean difference
#Use the posterior distribution of each pairwise comparison to calculate the probability that populations were different (posterior probability of difference, PD)
diffs <- population_difs_long %>%
  group_by(difference)%>%
  summarise(
    prob_greater_than_zero = mean(pred_diff > 0),
    prob_less_than_zero    = mean(pred_diff < 0),
    probability            = max(prob_greater_than_zero,
                                 prob_less_than_zero)
  )

population_difs_long <- population_difs_long %>%
  left_join(diffs, by = "difference")

#Supplemental Figure S3####
ggplot(population_difs_long, aes(x = pred_diff, y = difference, fill = difference)) +
  geom_density_ridges(fill = "gray", color = "darkgray", scale = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  
  theme_minimal(base_size = 20) +
  theme(legend.position = "none") +
  labs(x="Difference in Predicted Posteriors", y="Group Difference Density")+
  scale_y_discrete(limits = rev(c(  # Use `rev()` to flip order if needed
    "N1_C1", "N1_C2", "N2_C1", "N2_C2",
    "N1_S1", "N1_S2", "N2_S1", "N2_S2",
    "C1_S1", "C1_S2", "C2_S1", "C2_S2",
    "N1_N2", "C1_C2", "S1_S2"
  ))) + 
  geom_text(data = population_difs_long, aes(x = 0.12, y = difference, 
                                             label = paste0("PD(≠0) = ", round(probability, 2))),
            hjust = 0.7, vjust = 0.05, size = 6, color = "black", inherit.aes = FALSE)


#Pop-cohort effect (Q2)####
#Compute marginal predictions
post_preds_marginal <- post_preds_long_sla %>%
  group_by(draw, Population, Year) %>%  
  summarise(marginal_prediction = mean(prediction), .groups = "drop") 
head(post_preds_marginal); dim(post_preds_marginal)

#Specify the comparisons (i.e. is N1 2010 greater than N2 2010?)
comparisons <- list(
  #2010 comparisons within and across populations
  c("N1_2010", "N2_2010"), c("N1_2010", "C1_2010"), c("N1_2010", "C2_2010"),
  c("N1_2010", "S1_2010"), c("N1_2010", "S2_2010"), c("N2_2010", "C1_2010"),
  c("N2_2010", "C2_2010"), c("N2_2010", "S1_2010"), c("N2_2010", "S2_2010"),
  c("C1_2010", "C2_2010"), c("C1_2010", "S1_2010"), c("C1_2010", "S2_2010"),
  c("C2_2010", "S1_2010"), c("C2_2010", "S2_2010"), c("S1_2010", "S2_2010"),
  
  #2017 comparisons within and across populations
  c("N1_2017", "N2_2017"), c("N1_2017", "C1_2017"), c("N1_2017", "C2_2017"),
  c("N1_2017", "S1_2017"), c("N1_2017", "S2_2017"), c("N2_2017", "C1_2017"),
  c("N2_2017", "C2_2017"), c("N2_2017", "S1_2017"), c("N2_2017", "S2_2017"),
  c("C1_2017", "C2_2017"), c("C1_2017", "S1_2017"), c("C1_2017", "S2_2017"),
  c("C2_2017", "S1_2017"), c("C2_2017", "S2_2017"), c("S1_2017", "S2_2017"),
  
  #Cross-year comparisons
  c("N1_2010", "N1_2017"), c("N2_2010", "N2_2017"), c("C1_2010", "C1_2017"),
  c("C2_2010", "C2_2017"), c("S1_2010", "S1_2017"), c("S2_2010", "S2_2017")
)

#Pivot wider to create separate columns for each Population-Year group
post_preds_wide <- post_preds_marginal %>%
  pivot_wider(names_from = c(Population, Year), values_from = marginal_prediction)
head(post_preds_wide); dim(post_preds_wide)

#Compute row-wise differences for specified comparisons
df_differences <- post_preds_wide %>%
  rowwise() %>%
  mutate(
    !!!setNames(
      lapply(comparisons, function(cols) .[[cols[1]]] - .[[cols[2]]]),
      sapply(comparisons, function(cols) paste(cols, collapse = "-"))
    )
  ) %>%
  ungroup()
head(df_differences)

#probabilities - only keep columns 14- 49
df_differences <- df_differences[, 14:49]

difs_long <- df_differences %>%
  #select(14:49) %>%
  pivot_longer(cols = everything(), names_to = "var", values_to = "dif")

#Posterior mean difference
#Use the posterior distribution of each pairwise comparison to calculate the probability that populations were different (posterior probability of difference, PD
diffs <- difs_long %>% group_by(var) %>%
  summarise(
    prob_greater_than_zero = mean(dif > 0),
    prob_less_than_zero    = mean(dif < 0),
    probability            = max(prob_greater_than_zero,
                                 prob_less_than_zero)
  )
diffs

difs_long <- difs_long %>%
  left_join(diffs, by = "var") 

#Supplemental Figure S4####
ggplot(difs_long, aes(x = dif, y = var, fill = var)) +
  geom_density_ridges(fill = "gray", color = "darkgray", scale = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  
  theme_minimal(base_size = 16) +
  theme(legend.position = "none") +
  labs(x="Difference in Predicted Posteriors", y="Group Difference Density")+
  scale_y_discrete(limits = rev(c(  # Use `rev()` to flip order if needed
    "N1_2010-N2_2010", "C1_2010-C2_2010", "S1_2010-S2_2010",
    "N1_2017-N2_2017", "C1_2017-C2_2017", "S1_2017-S2_2017",
    "N1_2010-C1_2010", "N1_2010-C2_2010", "N2_2010-C1_2010", "N2_2010-C2_2010",
    "N1_2010-S1_2010", "N1_2010-S2_2010", "N2_2010-S1_2010", "N2_2010-S2_2010",
    "C1_2010-S1_2010", "C1_2010-S2_2010", "C2_2010-S1_2010", "C2_2010-S2_2010",
    "N1_2017-C1_2017", "N1_2017-C2_2017", "N2_2017-C1_2017", "N2_2017-C2_2017",
    "N1_2017-S1_2017", "N1_2017-S2_2017", "N2_2017-S1_2017", "N2_2017-S2_2017",
    "C1_2017-S1_2017", "C1_2017-S2_2017", "C2_2017-S1_2017", "C2_2017-S2_2017",
    "N1_2010-N1_2017", "N2_2010-N2_2017", "C1_2010-C1_2017", "C2_2010-C2_2017",
    "S1_2010-S1_2017", "S2_2010-S2_2017"
  ))) + 
  #geom_text(data = difs_long, aes(x = 3, y = as.numeric(var) + 1, 
  #                                label = paste0("P(pred_diff < 0) = ", round(probability, 3))),
  #         hjust = 0, vjust = 1, size = 3, color = "black", inherit.aes = FALSE)
  geom_text(data = difs_long, aes(x = 0.1, y = var, 
                                  label = paste0("PD(≠0) = ", round(probability, 2))),
            hjust = 0, vjust = 0.05, size = 5, color = "black", inherit.aes = FALSE)+
  coord_cartesian(clip = "off")




#LDMC plasticity brms model####
#check factors
ldmc_pp_mn$Population <- as.factor(ldmc_pp_mn$Population)
ldmc_pp_mn$Year <- as.factor(ldmc_pp_mn$Year)


mod13 <- brm(
  ldmc_plasticity ~ Population * Year # Fixed effects with interactions
  ,                                 
  data = ldmc_pp_mn,
  family = gaussian(),                           # Gaussian distribution for response variable
  chains = 4,                                    # Number of Markov chains
  cores = 4,                                     # Number of cores for parallel computation
  iter = 4000,                                   # Number of iterations per chain
)

summary(mod13) 
#trace plots
plot(mod13)

#saveRDS(mod13, file = "mod13brms.rds")
#mod13 <- readRDS("mod13brms.rds")

#Create empty dataframe
cohort_data_ldmc_pp <- expand.grid( 
  Population = unique(doy_pp_mn$Population)
  ,Year=unique(doy_pp_mn$Year)
)

#Make posterior predictions of the mean response for each group
post_preds_ldmc_pp <-as.data.frame(posterior_epred(mod13, newdata = cohort_data_ldmc_pp, allow_new_levels=TRUE, ndraws=1000) )
colnames(post_preds_ldmc_pp) <- paste(cohort_data_ldmc_pp$Population, cohort_data_ldmc_pp$Year, sep="_")
head(post_preds_ldmc_pp);dim(post_preds_ldmc_pp)

# Convert to long format
post_preds_long_ldmc_pp <- post_preds_ldmc_pp %>%
  mutate(draw = row_number()) %>%  # Add index for posterior draws
  pivot_longer(cols = -draw, names_to = "group", values_to = "prediction") %>%
  separate(group, into = c("Population", "Year"), sep = "_", remove = FALSE)
head(post_preds_long_ldmc_pp); dim(post_preds_long_ldmc_pp)

#Cohort effect (Q2)####
#Compute marginal predictions
post_preds_marginal <- post_preds_long_ldmc_pp %>%
  group_by(draw, Year) %>% 
  summarise(marginal_prediction = mean(prediction), .groups = "drop") 
head(post_preds_marginal); dim(post_preds_marginal)

difs_long <- post_preds_marginal %>%
  group_by(draw) %>% 
  summarise(diff_prediction = marginal_prediction[Year == 2010] - marginal_prediction[Year == 2017], .groups = "drop")

head(difs_long)

#Use the posterior distribution of each pairwise comparison to calculate the probability that cohorts were different (posterior probability of difference, PD)
cohort_probs <- difs_long %>%
  summarise(
    prob_greater_than_zero = mean(diff_prediction > 0),
    prob_less_than_zero    = mean(diff_prediction < 0),
    probability            = max(prob_greater_than_zero,
                                 prob_less_than_zero)
  )

difs_long <- difs_long %>%
  mutate(probability = cohort_probs$probability)

#plot overall cohort effect####
ggplot(difs_long, aes(x = diff_prediction)) +
  geom_density(fill = "gray",color = "darkgray", scale = 1) +  # Density plot
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Line at 0
  labs(
    #title = "Posterior Distribution of Differences (2010 - 2017)",
    #x = "Difference in predicted posteriors (2010 - 2017)",
    y = "Density"
  ) +
  theme_minimal(base_size = 20)+
  labs(x="Difference in Predicted Posteriors")+
  # Add the probability as text only once per facet
  geom_text(data = difs_long, aes(x = 0, y = 1.5, label = paste0("PD(≠0) = ", round(probability, 2))), 
            hjust = 0, vjust = 1, size = 6, color = "black")


#Population effect (Q1)####
#Compute marginal predictions
post_preds_marginal <- post_preds_long_ldmc_pp %>%
  group_by(draw, Population) %>%
  summarise(marginal_prediction = mean(prediction), .groups = "drop") 
head(post_preds_marginal); dim(post_preds_marginal)

#Take subtractions using a "wider" df
population_difs <- post_preds_marginal %>%
  pivot_wider(names_from = Population, values_from = marginal_prediction) %>%
  mutate(
    N1_C1 = N1 - C1, 
    N1_C2 = N1 - C2, 
    N2_C1 = N2 - C1, 
    N2_C2 = N2 - C2, 
    N1_S1 = N1 - S1, 
    N1_S2 = N1 - S2,
    N2_S1 = N2 - S1,
    N2_S2 = N2 - S2,
    C1_S1 = C1 - S1,
    C1_S2 = C1 - S2,
    C2_S1 = C2 - S1,
    C2_S2 = C2 - S2,
    N1_N2 = N1 - N2, 
    C1_C2 = C1 - C2, 
    S1_S2 = S1 - S2
  )
head(population_difs)

#pivot longer
population_difs_long <- population_difs %>%
  pivot_longer(
    cols = c(N1_C1, N1_C2, N2_C1, N2_C2,
             N1_S1, N1_S2, N2_S1, N2_S2,
             C1_S1, C1_S2, C2_S1, C2_S2,
             N1_N2, C1_C2, S1_S2),  # Columns to pivot
    names_to = "difference",   # New column for difference labels
    values_to = "pred_diff"    # New column for values
  )

#Posterior mean difference
#Use the posterior distribution of each pairwise comparison to calculate the probability that populations were different (posterior probability of difference, PD)
diffs <- population_difs_long %>%
  group_by(difference)%>%
  summarise(
    prob_greater_than_zero = mean(pred_diff > 0),
    prob_less_than_zero    = mean(pred_diff < 0),
    probability            = max(prob_greater_than_zero,
                                 prob_less_than_zero)
  )

population_difs_long <- population_difs_long %>%
  left_join(diffs, by = "difference")

#Supplemental Figure S5####
ggplot(population_difs_long, aes(x = pred_diff, y = difference, fill = difference)) +
  geom_density_ridges(fill = "gray", color = "darkgray", scale = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  
  theme_minimal(base_size = 20) +
  theme(legend.position = "none") +
  labs(x="Difference in Predicted Posteriors", y="Group Difference Density")+
  scale_y_discrete(limits = rev(c(  # Use `rev()` to flip order if needed
    "N1_C1", "N1_C2", "N2_C1", "N2_C2",
    "N1_S1", "N1_S2", "N2_S1", "N2_S2",
    "C1_S1", "C1_S2", "C2_S1", "C2_S2",
    "N1_N2", "C1_C2", "S1_S2"
  ))) + 
  geom_text(data = population_difs_long, aes(x = 0.18, y = difference, 
                                             label = paste0("PD(≠0) = ", round(probability, 2))),
            hjust = 0.7, vjust = 0.05, size = 6, color = "black", inherit.aes = FALSE)



#Pop-cohort effect (Q2)####
# Compute row-wise differences for specified comparisons
df_differences_ldmc_pp <- post_preds_ldmc_pp%>%
  rowwise() %>%
  mutate(
    !!!setNames(
      lapply(comparisons, function(cols) .[[cols[1]]] - .[[cols[2]]]),
      sapply(comparisons, function(cols) paste(cols, collapse = "-"))
    )
  ) %>%
  ungroup()

# View result
head(df_differences_ldmc_pp)

#probabilities
df_differences_ldmc_pp <- df_differences_ldmc_pp[, 13:48]

difs_long_ldmc_pp <- df_differences_ldmc_pp %>%
  pivot_longer(cols = everything(), names_to = "var", values_to = "dif")

#Posterior mean difference
#Use the posterior distribution of each pairwise comparison to calculate the probability that populations were different (posterior probability of difference, PD
diffs <- difs_long_ldmc_pp %>% group_by(var) %>%
  summarise(
    prob_greater_than_zero = mean(dif > 0),
    prob_less_than_zero    = mean(dif < 0),
    probability            = max(prob_greater_than_zero,
                                 prob_less_than_zero)
  )
diffs

difs_long_ldmc_pp <- difs_long_ldmc_pp %>%
  left_join(diffs, by = "var") 

#Supplemental figure S6####
ggplot(difs_long_ldmc_pp, aes(x = dif, y = var)) +
  geom_density_ridges(fill = "gray", color = "darkgray", scale = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  
  theme_minimal(base_size = 18) +
  theme(legend.position = "none") +
  labs(x="Difference in Predicted Posteriors", y="Group Difference Density")+
  scale_y_discrete(limits = rev(c(  # Use `rev()` to flip order if needed
    "N1_2010-N2_2010", "C1_2010-C2_2010", "S1_2010-S2_2010",
    "N1_2017-N2_2017", "C1_2017-C2_2017", "S1_2017-S2_2017",
    "N1_2010-C1_2010", "N1_2010-C2_2010", "N2_2010-C1_2010", "N2_2010-C2_2010",
    "N1_2010-S1_2010", "N1_2010-S2_2010", "N2_2010-S1_2010", "N2_2010-S2_2010",
    "C1_2010-S1_2010", "C1_2010-S2_2010", "C2_2010-S1_2010", "C2_2010-S2_2010",
    "N1_2017-C1_2017", "N1_2017-C2_2017", "N2_2017-C1_2017", "N2_2017-C2_2017",
    "N1_2017-S1_2017", "N1_2017-S2_2017", "N2_2017-S1_2017", "N2_2017-S2_2017",
    "C1_2017-S1_2017", "C1_2017-S2_2017", "C2_2017-S1_2017", "C2_2017-S2_2017",
    "N1_2010-N1_2017", "N2_2010-N2_2017", "C1_2010-C1_2017", "C2_2010-C2_2017",
    "S1_2010-S1_2017", "S2_2010-S2_2017"
  ))) + 
  xlim(-0.15,0.15)+
  #geom_text(data = difs_long, aes(x = 3, y = as.numeric(var) + 1, 
  #                                label = paste0("P(pred_diff < 0) = ", round(probability, 3))),
  #         hjust = 0, vjust = 1, size = 3, color = "black", inherit.aes = FALSE)
  geom_text(data = difs_long_ldmc_pp, aes(x = 0.12, y = var, 
                                  label = paste0("PD(≠0) = ", round(probability, 2))),
            hjust = 0.75, vjust = 0.03, size = 5, color = "black", inherit.aes = FALSE)+
  coord_cartesian(clip = "off")




