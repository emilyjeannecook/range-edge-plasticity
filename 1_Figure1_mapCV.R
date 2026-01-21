#Author: Emily Cook
#Project: PERSIST plasticity
#Purpose: Plot map and climate variables for M. cardinalis study populations (Fig. 1)


# Easy code for installing packages in R (if not installed) and calling their libraries
# From: https://gist.github.com/DrK-Lo/a945a29d6606b899022d0f03109b9483

#make vector of packages needed
packages_needed <- c("tidyverse","cowplot","RColorBrewer","ggmap","maps","mapdata","mapproj","raster","rnaturalearth","rnaturalearthdata","extrafont")

#install packages needed (if not already installed)
#for (i in 1:length(packages_needed)){
#  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
#}

#load packages needed
for (i in 1:length(packages_needed)){
  library( packages_needed[i], character.only = TRUE)
}

#data####
#read in cardinalis focal pops and common garden locations from DRYAD: https://doi.org/10.5061/dryad.18931zd7g
focal_cardinalis=read_csv("PERSIST_populations_gardens_1901-2021SY.csv") |>
 distinct(ID1,Latitude,Longitude) |>
 mutate(ID = ifelse(ID1=="N_garden"|ID1=="C_garden"|ID1=="S_garden","common garden","focal population")) # add a column to denote if location if focal population or common garden

#reorder populations
focal_cardinalis$ID1 = factor(focal_cardinalis$ID1,levels = c("S2","S1","S_garden","C2","C1","C_garden","N2","N1","N_garden")) 

#2010-2017 data
# read in .csv file of annual climate NA data (1951-2021) for each population and garden 
annual_clim=read_csv("PERSIST_populations_gardens_1901-2021SY.csv") |>
  mutate(ID1=factor(ID1),ID2=factor(ID2)) |>
  filter(Year > 2009 & Year < 2018, ID1 != "S_garden" & ID1 != "C_garden" & ID1 != "N_garden" ) |> # remove garden locations
  mutate(Year=factor(Year)) 
unique(annual_clim$Year)

#1981-2010 data for 30 year CV & anomalies
# read in .csv file of annual climate NA data (1951-2021) for each population and garden -- 30 year mean for anomaly
annual_clim30=read_csv("PERSIST_populations_gardens_1901-2021SY.csv") |>
  mutate(ID1=factor(ID1),ID2=factor(ID2)) |>
  filter(Year > 1980 & Year < 2011, ID1 != "S_garden" & ID1 != "C_garden" & ID1 != "N_garden" ) |> # remove garden locations
  mutate(Year=factor(Year)) 
unique(annual_clim30$Year)


#calculate anomalies####
# calculate mean winter precipitation 30 year
ppt30 <-  annual_clim30 %>%
  group_by(ID1)%>%
  summarize(mn_ppt30 = mean(PPT_wt))

# calculate mean winter ppt for 2012–2015
ppt1215 <- annual_clim %>%
  filter(Year %in% 2012:2015) %>%
  group_by(ID1, Year) %>%
  summarize(mn_ppt1215 = mean(PPT_wt), .groups="drop") #winter precip

# join and compute anomaly
ppt_anom <- left_join(ppt30, ppt1215, by="ID1") %>%
  mutate(anomaly = mn_ppt1215 - mn_ppt30)

#summarize
ppt_anom_summary <- ppt_anom %>%
  group_by(ID1) %>%
  summarise(mean_anomaly = mean(anomaly, na.rm = TRUE), .groups = "drop")

ppt_anom_summary


#compute CV####
# calculate 30-year means (1981 - 2010) for MAP and winter precipitation so that we can calculate climate CV
hist_clim <- annual_clim30 %>%
  group_by(ID1) %>%
  summarize(hist_MAP = mean(MAP), #mean annual precip
            hist_PPT_wt = mean(PPT_wt),#mean winter precip
            hist_MAT = mean(MAT), #mean annual temperature
            hist_Tmax_sm = mean(Tmax_sm),#max temp summer
            hist_CMD = mean(CMD), #climate moisture deficit
            cv_MAP = sd(MAP, na.rm = TRUE) / mean(MAP, na.rm = TRUE), # CV for MAP
            cv_PPTwt = sd(PPT_wt, na.rm = TRUE) / mean(PPT_wt, na.rm = TRUE), #CV for winter precip
            cv_MAT = sd(MAT, na.rm = TRUE) / mean(MAT, na.rm = TRUE))  # CV for MAT

#contemporary CV 2011 - 2017
contemp_clim <- annual_clim %>%
  group_by(ID1) %>%
  filter(Year!=2010)%>% #remove 2010
  summarize(hist_MAP = mean(MAP), #mean annual precip
            hist_PPT_wt = mean(PPT_wt),#mean winter precip
            hist_MAT = mean(MAT), #mean annual temperature
            hist_Tmax_sm = mean(Tmax_sm),#max temp summer
            hist_CMD = mean(CMD), #climate moisture deficit
            cv_MAP = sd(MAP, na.rm = TRUE) / mean(MAP, na.rm = TRUE), # CV for MAP
            cv_PPTwt = sd(PPT_wt, na.rm = TRUE) / mean(PPT_wt, na.rm = TRUE), #CV for winter precip
            cv_MAT = sd(MAT, na.rm = TRUE) / mean(MAT, na.rm = TRUE))  # CV for MAT

#join contemporary and historical dataframes: 
#Add period labels
hist_long <- hist_clim %>%
  mutate(period = "historical") %>%
  dplyr::select(ID1, period, cv_MAP, cv_PPTwt, cv_MAT)

contemp_long <- contemp_clim %>%
  mutate(period = "contemporary") %>%
  dplyr::select(ID1, period, cv_MAP, cv_PPTwt, cv_MAT)

clim_cv_long <- bind_rows(hist_long, contemp_long)

#2023 garden averages: 
#create data frame of 2023 common garden data (specifically mean annual temperature and winter precipitation)
#these were downloaded using climateNA v4.60 via web interface https://climatena.ca/mapversion on 20250509
garden = c("2023 north","2023 center","2023 south") 
garden_rank = c(3, 2, 1)
MAT=c(12,14.8,17.3)
PPT_wt=c(336,1119,490)
CMD=c(574,920,878)

Year = c(2016,2016,2016) # only for plotting purposes
clim_garden = data.frame(garden, garden_rank, MAT, PPT_wt,Year) |>
  mutate(Year = factor(Year))


#reorder populations
annual_clim$ID1 = factor(annual_clim$ID1,levels = c("S2","S1","C2","C1","N2","N1"))

#note: 2023 CMD for north garden was 574 (from climateNA)


#map of focal populations and common gardens
library(rnaturalearth)
library(terra)
library(raster)
# Assign colors per population
cols <- c("N1" = "#5E4FA2", 
          "N2" = "dodgerblue4", 
          "C1" = "gold", 
          "C2" = "goldenrod1", 
          "S1" = "red2", 
          "S2" = "#9E0142")
# get world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# create map of states
states <- map_data("state")

# make bounding box
ymin = 29.3499
ymax = 45.206755
xmin = -125.547702 
xmax = -114.05
e = extent(xmin, xmax,ymin, ymax)

#A. map####
#create a map of focal populations and occurrences without a legend  
card_map = ggplot(data=world,fill="lightgrey",col="black",size=0.3) + 
  geom_sf() +
  geom_polygon(aes(x = long, y = lat, group = group), data=states,fill="transparent",col="black",size=0.2) +
  theme_minimal() +
  coord_sf(xlim = c(xmin,xmax), ylim = c(ymin+0.5,ymax), expand = FALSE) +
  geom_point(aes(x=Longitude,y=Latitude, ,fill=ID1),data=filter(focal_cardinalis,ID=="focal population"),alpha=0.95,shape=21,size=4) +
  #geom_point(aes(x=Longitude,y=Latitude),data=filter(focal_cardinalis,ID=="common garden" & ID1=="N_garden"),alpha=1,shape=23,size=5,fill="black") +    
  labs(x="Longitude",y="Latitude") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16),legend.position = c(.99, .80),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        plot.title=element_text(hjust=0,size=18),
        legend.title = element_blank()) +
      scale_fill_manual(values = cols)+
  geom_text(aes(label = ID1, x = Longitude, y = Latitude), data = filter(focal_cardinalis, ID1 == "N1"), size = 4, fontface = "bold", nudge_x = 1, nudge_y = 0.2) +
  geom_text(aes(label = ID1, x = Longitude, y = Latitude), data = filter(focal_cardinalis, ID1 == "N2"), size = 4, fontface = "bold", nudge_x = -0.18, nudge_y = 0.5) +
  geom_text(aes(label = ID1, x = Longitude, y = Latitude), data = filter(focal_cardinalis, ID1 == "C1"), size = 4, fontface = "bold", nudge_x = 0.6, nudge_y = 0.2) +
  geom_text(aes(label = ID1, x = Longitude, y = Latitude), data = filter(focal_cardinalis, ID1 == "C2"), size = 4, fontface = "bold", nudge_x = 0.8, nudge_y = -0.1) +
  geom_text(aes(label = ID1, x = Longitude, y = Latitude), data = filter(focal_cardinalis, ID1 == "S1"), size = 4, fontface = "bold", nudge_x = .7, nudge_y = 0.3) +
  geom_text(aes(label = ID1, x = Longitude, y = Latitude), data = filter(focal_cardinalis, ID1 == "S2"), size = 4, fontface = "bold", nudge_x = 0.7, nudge_y = -0.2) +
  #geom_text(aes(label = "common garden", x = Longitude, y = Latitude), data = filter(focal_cardinalis, ID1 == "N_garden"), size = 4, fontface = "bold", nudge_x = 3.1, nudge_y = 0.2)+
  theme(legend.position = "none")
card_map

card_map_gards = card_map +
  geom_point(aes(x=Longitude,y=Latitude),
             data=filter(focal_cardinalis,ID=="common garden"),
             alpha=0.75,shape=24,size=5,fill="black") #add in triangles for common gardens

card_map_gards


#FIGURE of climatic moisture deficit of focal populations and northern common garden


#B. winter precip####
#plot Year vs. winter precipitation for multipanel figure 1
PPT_wt_plot=
  ggplot() + 
  geom_line(stat="smooth",data = annual_clim, aes(y = PPT_wt, x = Year, col=ID1, group=ID1), alpha = 0.8,linewidth = 2,se=FALSE) + 
  geom_hline(data=hist_clim,aes(yintercept = hist_PPT_wt, color = factor(ID1)), linewidth = 1.5, alpha = 0.5, linetype=2) +
  geom_text(data=clim_garden,aes(x=Year,y=PPT_wt,label=garden),size=4.5, fontface="bold") +
  xlab("Year") + 
  ylab("Winter precipitation (mm)") +
    ylim(0, 1600)+
  theme_bw(base_size=20)  +
  scale_color_manual(values = cols) +
  scale_x_discrete(breaks=c(2010,2012,2014,2016)) +
  theme(legend.position="none",axis.text=element_text(size=12),
        axis.title=element_text(size=16),legend.title=element_blank(),legend.text = element_text(size = 20)) 
PPT_wt_plot


#plot precipitation anomaly
ggplot(ppt_anom_summary, aes(x = ID1, y = mean_anomaly, col = ID1, group = ID1)) + 
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = cols, guide = "none") +  # <- hide legend here
    ylab("Winter\nprecipitation anomaly (mm)") +
    xlab("Population") +
    ylim(-300, 0) +
    theme_bw(base_size = 16)

#C. historical vs current CV####
#plot historical and contemporary side by side
#order from highest to lowest CV
clim_cv_long$ID1 = factor(clim_cv_long$ID1,levels = c("N1","N2","C2","C1","S2","S1"))
clim_cv_long$period <- factor(
  clim_cv_long$period,
  levels = c("historical", "contemporary"),
  labels = c("Historical", "Contemporary")
)

precip_CV_plot <- ggplot(clim_cv_long, aes(x = ID1,
                         y = cv_PPTwt,
                         color = ID1,
                         alpha = period)) +
  geom_point(size = 4,
             position = position_nudge(x = c(-0.07, 0.07))) +  # separates hist/contemp a bit
  scale_color_manual(values = cols) +
  labs(x = "Source Population",
       y = "CV in winter precipitation") +
  scale_alpha_manual(
    name = "Period",  #capitalize legend
    values = c(Historical = 0.55, Contemporary = 1),
    labels = c(Historical = "1981–2010", #change labels to actual years on legend
               Contemporary = "2011–2017")
  ) +
  ylim(0.2, 0.75) +
  theme_bw(base_size = 14)+
  theme( # show legend
    legend.position = c(0.84, 0.15),   # inside plot, bottom-right
    legend.background = element_rect(fill = alpha("white", 0.7), colour = "black"), #fill legend so no grid behind
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10)
  )+      
  guides(color = "none")+                   # hide color legend, keep alpha legend
  facet_wrap(~ period)
precip_CV_plot 




#join plots####
#add panel labels
card_map_A = card_map_gards +
  ggtitle('A') + theme(plot.title=element_text(hjust=0,size=20)) 

#winter precip
PPT_wt_plot_B = PPT_wt_plot +
  ggtitle('B') + theme(plot.title=element_text(hjust=0,size=20)) 

#historical vs contemporary CV
precip_CV_plot_C = precip_CV_plot +
  ggtitle('C') + theme(plot.title=element_text(hjust=0,size=20)) 


#font size & layout
card_map_A        <- card_map_A        + theme_bw(base_size = 14)
PPT_wt_plot_B     <- PPT_wt_plot_B     + theme_bw(base_size = 14)
precip_CV_plot_C  <- precip_CV_plot_C  + theme_bw(base_size = 14)
#legends
card_map_A <- card_map_A + theme(legend.position = "none")
PPT_wt_plot_B <- PPT_wt_plot_B + theme(legend.position = "none")
precip_CV_plot_C <- precip_CV_plot_C + theme(legend.position = "none")

#stack B and C on the right
right_panel <- plot_grid(
  PPT_wt_plot_B,
  precip_CV_plot_C,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1)
)

#combine A with the stacked right panels
final_figure <- plot_grid(
  card_map_A,
  right_panel,
  ncol = 2,
  rel_widths = c(0.9,1)  # tweak to match visual balance
)

final_figure

#add white rectangle
final_figure_white <- ggdraw() +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf,
    fill = "white",
    color = NA
  ) +
  draw_plot(final_figure)

#Figure 1####
final_figure_white




