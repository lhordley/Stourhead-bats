##########################
#### user: Lisbeth Hordley
#### date: November 2021
#### info: Stourhead bats: calculate richness and activity to use in models

rm(list = ls())
options(scipen=999)

# Packages
library(plyr)
library(dplyr)
library(ggplot2)
library(devtools)
library(ggiraphExtra)
library(ggeffects)
library(tidyr)
library(MASS)       
library(tidyverse)

## read in data
bat_data <- read.csv("Data/bat_data_cleaned_final.csv", header=TRUE)
hab_data <- read.csv("Data/Stourhead_hab_structure_draft.csv", header=TRUE)

#####################################################################
## Data Prep

# for some reason if time is at midnight (00:xx), the date is taken as the next day (but not if time is 1am onwards)
# need to fix this so there are only two dates (i.e. two visits) for each plot

# if first two characters of time = 00, then minus one number from day
bat_data$day <- ifelse(grepl('^00', bat_data$time), bat_data$day-1, bat_data$day)

## create visit column: 1 = first visit, 2 = second visit for each plot
bat_data <- bat_data %>% group_by(plot) %>%
  mutate(visit=ifelse(day==min(day), 1, 2))

## calculate total species richness & richness of guilds
library(dplyr)
bat_data_rich <- bat_data %>% 
  group_by(plot, visit, year) %>%
  summarise(richness = n_distinct(scientific_name),
            sre_rich = n_distinct(scientific_name[guild=="SRE"], na.rm=TRUE),
            mre_rich = n_distinct(scientific_name[guild=="MRE"], na.rm=TRUE),
            lre_rich = n_distinct(scientific_name[guild=="LRE"], na.rm=TRUE))

## merge hab data with bat data
bat_hab_data <- merge(bat_data_rich, hab_data, by.x=c("plot","year"), by.y=c("plot", "survey_year"))

bat_hab_data[,c("location","dbh_tree_1","dbh_tree_2","dbh_tree_3","dbh_tree_4","dbh_tree_5",
                "tree_species_names", "canopy_.cover", "sub_canopy_.cover", "tall_understorey_cover",
                "low_shrub_cover", "ground_layer_cover", "canopy_layer_score", "sub_canopy_layer_score",
                "tall_understorey_score", "low_shrub_score", "ground_layer_score", "photo_number", "comments")] <- list(NULL)
unique(bat_hab_data$treatment) ## 3 treatments: early transitioning irregular, irregular, and clear fell
## 64 plots
## change habitat codes
bat_hab_data$treatment <- recode(bat_hab_data$treatment, "Clear Fell" = "PS", "Early Transitioning Irregular" = "RIS", "Irregular" = "SDS")

bat_hab_data$visit <- as.factor(bat_hab_data$visit)
bat_hab_data$subcompartment <- as.factor(bat_hab_data$subcompartment)

res <- bartlett.test(richness ~ treatment, data = bat_hab_data)
res ## non-significant: variance in richness is the same between treatments

## simple plot of richness between treatments
ggplot(bat_hab_data, aes(x=treatment, y=richness, fill=visit))+
  geom_boxplot()+ 
  theme_classic() + 
  ylab("Species richness")

ggplot(bat_hab_data, aes(treatment, sre_rich, fill=visit))+
  geom_boxplot()+ 
  theme_classic() + 
  ylab("Species richness")

ggplot(bat_hab_data, aes(treatment, mre_rich, fill=visit))+
  geom_boxplot()+ 
  theme_classic() + 
  ylab("Species richness")

ggplot(bat_hab_data, aes(treatment, lre_rich, fill=visit))+
  geom_boxplot()+ 
  theme_classic() + 
  ylab("Species richness")

## save file
write.csv(bat_hab_data, file="Data/Stourhead_bats_richness_final.csv", row.names=FALSE)

######################################################################
## Calculate bat activity for each species

## 1. Occupancy rate: proportion of 10-minute periods during a night, within which the bat species was recorded at least once
## this accounts for night length differences

## number of 10minute periods each night (e.g. 8 hours: 48)
## number of 10minute periods where, e.g. common pip, is recorded at least once (e.g. 30)
## 30/48: 0.625 (62.5%)
library(base)
bat_data$time2 <- as.POSIXct(bat_data$time, format = "%H:%M:%S")
hrs <- 12 * 60 * 60
bat_data$time3 <- bat_data$time2 + hrs
library(data.table)
bat_data$time4 <- as.ITime(bat_data$time3) 

bat_data$time4 <- strftime(bat_data$time4, format="%H:%M:%S")
bat_data$time4 <- as.POSIXct(bat_data$time4, format="%H:%M:%S")
bat_data$warnings <- NULL

plot <- unique(bat_data$plot)
visit <- unique(bat_data$visit)
bat_activity_rate <- NULL

for (g in plot){ # loop through each plot
  print(g) 
  plot_data <- bat_data[bat_data$plot==g,]
  
  for (i in visit){ # loop through each visit
    
    sub_data <- plot_data[plot_data$visit==i,]
    sub_data$ten_min_group <- cut(sub_data$time4, breaks="10 mins") ## match each row with 10 minute categories
    MyDatesTable <- table(cut(sub_data$time4, breaks="10 mins")) ## final ALL 10 min categories 
    MyDates <- data.frame(MyDatesTable)
    sub_data <- merge(sub_data, MyDates, by.x="ten_min_group", by.y="Var1", all.y=TRUE) ## merge together
    tot_ten_min_intervals <- length(unique(MyDates$Var1)) ## 57
    
    sub_data <- na.omit(sub_data)
    species <- unique(sub_data$scientific_name)
    
    for(h in species){
      
      sub_data_sp <- subset(sub_data, scientific_name == h)
      sp_ten_min_intervals <- length(unique(sub_data_sp$ten_min_group)) ## 6
      activity_rate <- sp_ten_min_intervals / tot_ten_min_intervals
      
      tmp_df <- data.frame(plot=g, visit=i, scientific_name=h, tot_intervals = tot_ten_min_intervals,
                           sp_intervals = sp_ten_min_intervals, activity_rate=activity_rate)
      bat_activity_rate <- rbind(tmp_df, bat_activity_rate)
      
    }
  }
}

## merge with hab data
## merge hab data with moth data
bat_activity_rate$year <- 2021
bat_act_hab <- merge(bat_activity_rate, hab_data, by.x=c("plot","year"), by.y=c("plot", "survey_year"))

bat_act_hab[,c("location","dbh_tree_1","dbh_tree_2","dbh_tree_3","dbh_tree_4","dbh_tree_5",
               "tree_species_names", "canopy_.cover", "sub_canopy_.cover", "tall_understorey_cover",
               "low_shrub_cover", "ground_layer_cover", "canopy_layer_score", "sub_canopy_layer_score",
               "tall_understorey_score", "low_shrub_score", "ground_layer_score", "photo_number", "comments")] <- list(NULL)
unique(bat_act_hab$treatment) ## 3 treatments: early transitioning irregular, irregular, and clear fell
## 64 plots

bat_act_hab$visit <- as.factor(bat_act_hab$visit)
bat_act_hab$subcompartment <- as.factor(bat_act_hab$subcompartment)

## run model for each species

## remove species with only very few records in plots - Nycatalus leisleri, pipistrellus nathusii, rhinolophus ferrumequinum and rhinolophus hipposideros
bat_act_hab3 <- subset(bat_act_hab, !(scientific_name %in% c("Nyctalus leisleri", "Pipistrellus nathusii",
                                                             "Rhinolophus ferrumequinum", "Rhinolophus hipposideros")))
# sp_results <- lapply(split(bat_act_hab3, bat_act_hab3$scientific_name), function(x) summary(glmer(
#   activity_rate ~ treatment + (1|visit) + (1|subcompartment),family=Gamma(link="log"),data=x)))
## quick way to get results by species - but need to check model assumptions & do post-hoc testing on each separately 

bat_act_hab3 <- bat_act_hab3 %>% mutate_all(na_if,"n/a")
bat_act_hab3$basal_area <- as.numeric(bat_act_hab3$basal_area)
bat_act_hab3$canopy_openess <- as.numeric(bat_act_hab3$canopy_openess)
bat_act_hab3$treatment <- recode(bat_act_hab3$treatment, "Clear Fell" = "PS", "Early Transitioning Irregular" = "RIS", "Irregular" = "SDS")

## save file
write.csv(bat_act_hab3, file="Data/Stourhead_bats_activity_final.csv", row.names=FALSE)
