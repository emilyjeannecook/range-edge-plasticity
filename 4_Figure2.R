#Author: Emily Cook
#Project: PERSIST plasticity
#Purpose: Figure 2, plasticity in first flower date, sla, and ldmc across populations and cohorts

#packages
library(tidyverse)
library(colorspace)
library(ggdist)
library(patchwork)


# Assign colors per population
cols <- c("N1" = "#5E4FA2", 
          "N2" = "dodgerblue4", 
          "C1" = "gold", 
          "C2" = "goldenrod1", 
          "S1" = "red2", 
          "S2" = "#9E0142")



#load dataframes from 2_calculate_plasticity.R
#doy_pp_mn <- read.csv("doy_pp.csv") #read
#sla_pp_mn <- read.csv("sla_pp.csv") #read
#ldmc_pp_mn <- read.csv("ldmc_pp.csv") #read

#set factors
doy_pp_mn$Year <- as.factor(doy_pp_mn$Year)
sla_pp_mn$Year <- as.factor(sla_pp_mn$Year)
ldmc_pp_mn$Year <- as.factor(ldmc_pp_mn$Year)



#First flower date plasticity####

#Load in needed models to assess significant cohort differences from 3_Q1_Q2_brms_FiguresS1-S6.R
mod11 <- readRDS("mod11brms.rds")

#Get differences for each pop-cohort
#create empty dataframe
cohort_data <- expand.grid( 
  Population = unique(doy_pp_mn$Population)
  ,Year=unique(doy_pp_mn$Year)
)

#Get Posterior predictions of the mean response for each group
post_preds <-as.data.frame(posterior_epred(mod11, newdata = cohort_data, allow_new_levels=TRUE, ndraws=1000) )
colnames(post_preds) <- paste(cohort_data$Population, cohort_data$Year, sep="_")

head(post_preds);dim(post_preds)

#Convert to long format
post_preds_long <- post_preds %>%
  mutate(draw = row_number()) %>%  # Add index for posterior draws
  pivot_longer(cols = -draw, names_to = "group", values_to = "prediction") %>%
  separate(group, into = c("Population", "Year"), sep = "_", remove = FALSE)

#pop-cohort effect: Compute marginal predictions
post_preds_marginal <- post_preds_long %>%
  group_by(draw, Population, Year) %>%  # Group by draw, cohort
  summarise(marginal_prediction = mean(prediction), .groups = "drop") 

# Specify the comparisons (i.e. is N1 2010 greater than N2 2010?)
comparisons <- list(
  #2010 comparisons within and across populations
  c("N1_2010", "N2_2010"), c("N1_2010", "C1_2010"), c("N1_2010", "C2_2010"),
  c("N1_2010", "S1_2010"), c("N1_2010", "S2_2010"), c("N2_2010", "C1_2010"),
  c("N2_2010", "C2_2010"), c("N2_2010", "S1_2010"), c("N2_2010", "S2_2010"),
  c("C1_2010", "C2_2010"), c("C1_2010", "S1_2010"), c("C1_2010", "S2_2010"),
  c("C2_2010", "S1_2010"), c("C2_2010", "S2_2010"), c("S1_2010", "S2_2010"),
  
  # 2017 comparisons within and across populations
  c("N1_2017", "N2_2017"), c("N1_2017", "C1_2017"), c("N1_2017", "C2_2017"),
  c("N1_2017", "S1_2017"), c("N1_2017", "S2_2017"), c("N2_2017", "C1_2017"),
  c("N2_2017", "C2_2017"), c("N2_2017", "S1_2017"), c("N2_2017", "S2_2017"),
  c("C1_2017", "C2_2017"), c("C1_2017", "S1_2017"), c("C1_2017", "S2_2017"),
  c("C2_2017", "S1_2017"), c("C2_2017", "S2_2017"), c("S1_2017", "S2_2017"),
  
  # Cross-year comparisons
  c("N1_2010", "N1_2017"), c("N2_2010", "N2_2017"), c("C1_2010", "C1_2017"),
  c("C2_2010", "C2_2017"), c("S1_2010", "S1_2017"), c("S2_2010", "S2_2017")
)

# Pivot wider to create separate columns for each Population-Year group
post_preds_wide <- post_preds_marginal %>%
  pivot_wider(names_from = c(Population, Year), values_from = marginal_prediction)

# Compute row-wise differences for specified comparisons
df_differences <- post_preds_wide %>%
  rowwise() %>%
  mutate(
    !!!setNames(
      lapply(comparisons, function(cols) .[[cols[1]]] - .[[cols[2]]]),
      sapply(comparisons, function(cols) paste(cols, collapse = "-"))
    )
  ) %>%
  ungroup()

# View result
head(df_differences)

#probabilities - only keep columns 14- 49
df_differences <- df_differences[, 14:49]

difs_long <- df_differences %>%
  #select(14:49) %>%
  pivot_longer(cols = everything(), names_to = "var", values_to = "dif")

#compute number of differences that are diff than zero: the probabilities
diffs <- difs_long %>% group_by(var) %>%
  summarise(
    prob_greater_than_zero = mean(dif > 0),
    prob_less_than_zero    = mean(dif < 0),
    probability            = max(prob_greater_than_zero, prob_less_than_zero)
  )
diffs

difs_long <- difs_long %>%
  left_join(diffs, by = "var") #add "probability" back in to raw data for graph

#make labels for Figure 1 with asterisks for significant cohort differences: 
#get a probability that posterior difference is different than zero (not less or greater than zero)
diffs_labels <- difs_long %>%
  group_by(var) %>%
  summarise(
    mean_diff = mean(dif),                          # average difference - dont need this for asterisks
    prob_greater_than_zero = mean(dif > 0),         # % of draws > 0
    prob_less_than_zero  = mean(dif < 0),           # % of draws < 0
    prob_diff_from_zero  = max(prob_greater_than_zero, prob_less_than_zero),
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      prob_diff_from_zero >= 0.999 ~ "***",
      prob_diff_from_zero >= 0.99  ~ "**",
      prob_diff_from_zero >= 0.95  ~ "*",
      TRUE                         ~ ""
    )
  )
#keep within population contrasts only:
contrasts_keep <- c("N1_2010-N1_2017", "N2_2010-N2_2017",
                    "C1_2010-C1_2017", "C2_2010-C2_2017",
                    "S1_2010-S1_2017", "S2_2010-S2_2017")

pvals_df_sub <- diffs_labels |>
  dplyr::filter(var %in% contrasts_keep) #filter to just cohort differences 
pvals_df_sub

pvals_df_sub <- pvals_df_sub |>
  dplyr::mutate(Population = stringr::str_extract(var, "^[A-Z]\\d")) #add in a population column to be able to map to ggplot

#get posterior 95% credible intervals
# summarise posterior draws from brms
posterior_summary <- post_preds_marginal %>%
  group_by(Population, Year) %>%
  summarise(
    median = median(marginal_prediction),
    lower  = quantile(marginal_prediction, 0.025),
    upper  = quantile(marginal_prediction, 0.975)
  )


#Make letters dataframe showing significant diffs between populations
letters_df <- tibble(
  Population = factor(c("N1", "N2", "C1", "C2", "S1", "S2"),
                      levels = c("N1","N2","C1","C2","S1","S2")),
  letters = c("a", "a", "b", "c", "c", "c"),
  ypos = 1.15  # adjust as needed
)

# Y-position for stars showing significant cohort differences
ypos_stars <- 1.05

#rank from low to high CV
doy_pp_mn <- doy_pp_mn %>% mutate(Population=fct_relevel(Population, "S1", "S2", "C1", "C2", "N2", "N1"))
doy_pp_mn <- doy_pp_mn %>% mutate(Year=fct_relevel(Year, "2017", "2010")) #ancestors "on top"

#Use custom labels (only label the first occurrence)
pop_labels <- c("N1", "", "N2", "", "C1", "", "C2", "", "S1", "", "S2", "")

#First flower date plot####
doy_plot <- ggplot(doy_pp_mn, aes(x = Population, y = doy_plasticity,
                                  fill = Population, alpha = factor(Year))) +
  geom_boxplot(outlier.shape = NA, width = 0.5,
               position = position_dodge(width = 0.65)) +
  geom_jitter(shape = 21, color = "black", size = 1,
              position = position_jitterdodge(dodge.width = 0.65,
                                              jitter.width = 0.65)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Source Population", y = "Adjusted First Flower Date Plasticity",
       alpha = "Year") +
  scale_fill_manual(values = cols) +
  scale_alpha_manual(values = c("2010" = 0.35, "2017" = 0.75)) +
  scale_y_continuous(limits = c(-0.25, 1.2),
                     breaks = c(-0.25, 0, 0.5, 1)) +
  theme_classic() +
  theme(legend.position = "none") +
  # add stars
  geom_text(data = pvals_df_sub,
            aes(x = Population, y = ypos_stars, label = significance),
            inherit.aes = FALSE,
            size = 6,            # make text bigger
            fontface = "bold")+
  # Letters for population differences
  geom_text(data = letters_df,
            aes(x = Population, y = ypos, label = letters),
            inherit.aes = FALSE,
            size = 4,)+
  theme(
    axis.text = element_text(size = 12),       # axis tick labels
    axis.title = element_text(size = 14),      # x and y titles
    #plot.title = element_text(size = 16, face = "bold"), # plot title
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  )+
  coord_flip()  # flips x and y axes

doy_plot

#Add plot label
doy_plot_A = doy_plot +
  ggtitle('A') + theme(plot.title=element_text(hjust=0,size=20)) 
doy_plot_A


#SLA plasticity####
#no cohort effects

#Set the desired order of x-axis (factor levels)
sla_pp_mn <- sla_pp_mn %>% mutate(Population=fct_relevel(Population, "N1", "N2", "C1", "C2", "S1", "S2"))

sla_pp_mn$pop_year <- factor(sla_pp_mn$pop_year, levels = c(
  "N1-2010", "N1-2017",
  "N2-2010", "N2-2017",
  "C1-2010", "C1-2017",
  "C2-2010", "C2-2017",
  "S1-2010", "S1-2017",
  "S2-2010", "S2-2017"
))

# Letters dataframe showing significant population differences
letters_df <- tibble(
  Population = factor(c("N1", "N2", "C1", "C2", "S1", "S2"),
                      levels = c("N1","N2","C1","C2","S1","S2")),
  letters = c("a", "a", "b", "b", "b", "c"),
  ypos = 1.15  # adjust as needed
)

#rank from low to high CV
sla_pp_mn <- sla_pp_mn %>% mutate(Population=fct_relevel(Population, "S1", "S2", "C1", "C2", "N2", "N1"))
sla_pp_mn <- sla_pp_mn %>% mutate(Year=fct_relevel(Year, "2017", "2010")) #ancestors "on top"

#SLA plot####
plot_sla <- ggplot(sla_pp_mn, aes(x = Population, y = sla_plasticity, fill = Population, alpha = factor(Year))) +
  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.65)) +
  geom_jitter(shape = 21, color = "black", size = 1, position = position_jitterdodge(dodge.width = 0.65, jitter.width = 0.65,)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Source Population", y = "Specific Leaf Area Plasticity", alpha = "Year") +
  scale_fill_manual(values = cols) +
  scale_alpha_manual(values = c("2010" = 0.35, "2017" = 0.75)) +
  #scale_x_discrete(labels = pop_labels) +
  #scale_y_continuous(breaks = seq(floor(min(sla_pp_mn$sla_plasticity)),
  #                                ceiling(max(sla_pp_mn$sla_plasticity)), 
  #                               by = 0.25)) +
  theme_classic()+
  scale_y_continuous(limits = c(-0.25, 1.2),
                     breaks = c(-0.25, 0, 0.5, 1)) +
  theme(legend.position = "none")+
  theme(
    axis.text = element_text(size = 12),       # axis tick labels
    axis.title = element_text(size = 14),      # x and y titles
    #plot.title = element_text(size = 16, face = "bold"), # plot title
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  )+
  # Letters for population differences
  geom_text(data = letters_df,
            aes(x = Population, y = ypos, label = letters),
            inherit.aes = FALSE,
            size = 4,)+
  coord_flip()  # flips x and y axes

plot_sla 
#add plot label
sla_plot_B = plot_sla +
  ggtitle('B') + theme(plot.title=element_text(hjust=0,size=20)) 
sla_plot_B


#LDMC plasticity####
#no cohort effect

#Set the desired order of x-axis (factor levels)
ldmc_pp_mn$pop_year <- factor(ldmc_pp_mn$pop_year, levels = c(
  "N1-2010", "N1-2017",
  "N2-2010", "N2-2017",
  "C1-2010", "C1-2017",
  "C2-2010", "C2-2017",
  "S1-2010", "S1-2017",
  "S2-2010", "S2-2017"
))

ldmc_pp_mn <- ldmc_pp_mn %>% mutate(Year=fct_relevel(Year, "2017", "2010")) #ancestors "on top"


#Letters dataframe showing significant population differences
letters_df <- tibble(
  Population = factor(c("N1", "N2", "C1", "C2", "S1", "S2"),
                      levels = c("N1","N2","C1","C2","S1","S2")),
  letters = c("ab", "a", "ab", "b", "ab", "c"),
  ypos = 1.15  # adjust as needed
)

#rank from low to high CV
ldmc_pp_mn <- ldmc_pp_mn %>% mutate(Population=fct_relevel(Population, "S1", "S2", "C1", "C2", "N2", "N1"))

#LDMC plot####
plot_ldmc <- ggplot(ldmc_pp_mn, aes(x = Population, y = ldmc_plasticity, fill = Population, alpha = factor(Year))) +
  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.65)) +
  geom_jitter(shape = 21, color = "black", size = 1,
              position = position_jitterdodge(dodge.width = 0.65, jitter.width = 0.65)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Source Population", y = "Leaf Dry Matter Content Plasticity", alpha = "Year") +
  scale_fill_manual(values = cols, guide = guide_legend(reverse = TRUE, order = 1)) + #flip order of legend, N1 on top
  #scale_alpha_manual(values = c("2010" = 0.35, "2017" = 0.75), #alpha = year
  #                  guide = guide_legend(order = 2, override.aes = list(fill = "black", shape = 21))) +
  scale_alpha_manual(
    values = c("2010" = 0.35, "2017" = 0.75),
    breaks = c("2010", "2017"),   # legend order ancestors on top
    name = "Cohort",
    guide = guide_legend(
      order = 2,
      override.aes = list(fill = "black", shape = 21)
    )
  )+
  scale_y_continuous(limits = c(-0.25, 1.2), breaks = c(-0.25, 0, 0.5, 1)) + #set y axis ticks
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  ) +
  geom_text(data = letters_df, #add in letters for significant differencces
            aes(x = Population, y = ypos, label = letters),
            inherit.aes = FALSE,
            size = 4) +
  coord_flip() #flip x and yaxis

plot_ldmc 
#add plot label
ldmc_plot_C = plot_ldmc +
  ggtitle('C') + theme(plot.title=element_text(hjust=0,size=20)) 
ldmc_plot_C


#join plasticity plots####
plasticity <- doy_plot_A / sla_plot_B / ldmc_plot_C #vertical
plasticity












