##########################
#### user: Lisbeth Hordley
#### date: November 2021
#### info: Stourhead bats: run models testing effect of habitat and surrounding landscape

rm(list = ls())
options(scipen=999)

# Packages
library(MuMIn)
library(DHARMa)
library(lme4)
library(StatisticalModels)
library(glmmTMB)
library(buildmer)
library(dplyr)
library(ggeffects)
library(easystats)

# read in data
bat_rich <- read.csv("Data/Stourhead_bats_richness_final.csv", header=TRUE)
bat_act <- read.csv("Data/Stourhead_bats_activity_final.csv", header=TRUE)
land_cover <- read.csv("Data/Land_cover_stourhead_2015.csv", header=TRUE)
# change orientation of land cover data
colnames(land_cover) <- c("site","scale", "prop_arable","prop_broadleaf","prop_conifer","prop_grassland")
land_cover <- reshape(land_cover, idvar = "site", timevar = "scale", direction = "wide")

# merge in land cover data
bat_rich <- merge(bat_rich, land_cover, by.x="gridref", by.y="site", all.x=TRUE)
bat_act <- merge(bat_act, land_cover, by.x="gridref", by.y="site", all.x=TRUE)

# Make sure variables are in right format
bat_rich <- bat_rich %>% mutate_all(na_if,"n/a")
bat_rich$basal_area <- as.numeric(bat_rich$basal_area)
bat_rich$canopy_openess <- as.numeric(bat_rich$canopy_openess)
summary(bat_rich) # no NAs - 32 plots with two visits each

########## SURROUNDING LANDSCAPE AND HABITAT MODELS
# This is to determine which parameters (radius of each land cover type) will be included in the habitat variable models
# Following Fuentes-Montemayor et al. 2012 paper methods
# Then including the relevant land cover type/radius in models with habitat variables and compare AIC

round(cor(bat_rich[,c("dbh_average","basal_area","per_broadleaf_canopy","complexity_score","canopy_openess",
                      "deadwood_snags","fallen_deadwood", "dis_to_edge", "prop_grassland.1500")]),3)
cor.test(bat_rich$per_broadleaf_canopy, bat_rich$deadwood_snags) ## significant and r>0.8
## include the variable most strongly correlated with response (sp richness)
cor.test(bat_rich$richness, bat_rich$per_broadleaf_canopy) # r=0.064, p=0.61
cor.test(bat_rich$richness, bat_rich$deadwood_snags) # r=0.061, p=0.63
## remove deadwood snags from analysis

# Total richness
tot_rich_arable1 <- glmer(richness ~ prop_arable.500 + visit + (1|subcompartment),
                          family=poisson(), data=bat_rich)
summary(tot_rich_arable1)
r.squaredGLMM(tot_rich_arable1) # 0.003728972    
tot_rich_arable2 <- glmer(richness ~ prop_arable.1500 + visit + (1|subcompartment),
                          family=poisson(), data=bat_rich)
summary(tot_rich_arable2) 
r.squaredGLMM(tot_rich_arable2) # 0.0001509244   
tot_rich_arable3 <- glmer(richness ~ prop_arable.3000 + visit + (1|subcompartment),
                          family=poisson(), data=bat_rich)
summary(tot_rich_arable3)    
r.squaredGLMM(tot_rich_arable3) # 0.003820198   
#### ARABLE 3000M #### (previously arable 500m)

tot_rich_broad1 <- glmer(richness ~ prop_broadleaf.500 + visit + (1|subcompartment),
                         family=poisson(), data=bat_rich)
summary(tot_rich_broad1) 
r.squaredGLMM(tot_rich_broad1) # 0.009587220     
tot_rich_broad2 <- glmer(richness ~ prop_broadleaf.1500 + visit + (1|subcompartment),
                         family=poisson(), data=bat_rich)
summary(tot_rich_broad2)
r.squaredGLMM(tot_rich_broad2) # 0.005287890     
tot_rich_broad3 <- glmer(richness ~ prop_broadleaf.3000 + visit + (1|subcompartment),
                         family=poisson(), data=bat_rich)
summary(tot_rich_broad3) 
r.squaredGLMM(tot_rich_broad3) # 0.01426713    
#### BROADLEAF 3000M ####  (same as before)


tot_rich_conifer1 <- glmer(richness ~ prop_conifer.500 + visit + (1|subcompartment),
                           family=poisson(), data=bat_rich)
summary(tot_rich_conifer1)   
r.squaredGLMM(tot_rich_conifer1) # 0.002762292      
tot_rich_conifer2 <- glmer(richness ~ prop_conifer.1500 + visit + (1|subcompartment),
                           family=poisson(), data=bat_rich)
summary(tot_rich_conifer2)       
r.squaredGLMM(tot_rich_conifer2) # 0.03986401      
tot_rich_conifer3 <- glmer(richness ~ prop_conifer.3000 + visit + (1|subcompartment),
                           family=poisson(), data=bat_rich)
summary(tot_rich_conifer3)         
r.squaredGLMM(tot_rich_conifer3) # 0.03763621     
#### CONIFER 1500M #### (previously conifer 3000m) 

tot_rich_grass1 <- glmer(richness ~ prop_grassland.500 + visit + (1|subcompartment),
                         family=poisson(), data=bat_rich)
summary(tot_rich_grass1)        
r.squaredGLMM(tot_rich_grass1) # 0.0007923292       
tot_rich_grass2 <- glmer(richness ~ prop_grassland.1500 + visit + (1|subcompartment),
                         family=poisson(), data=bat_rich)
summary(tot_rich_grass2)                   
r.squaredGLMM(tot_rich_grass2) # 0.03962581       
tot_rich_grass3 <- glmer(richness ~ prop_grassland.3000 + visit + (1|subcompartment),
                         family=poisson(), data=bat_rich)
summary(tot_rich_grass3)          
r.squaredGLMM(tot_rich_grass3) # 0.01488594      
#### GRASSLAND 1500M #### (same as before)
# same buffers selected when visit moved from random to fixed

summary(tot_rich_arable3)$coefficients[1,1]
tot_rich_arable3$terms
names(model.frame(tot_rich_arable3))[2]

# save results of lowest R2 model for each land use type
tot_rich_land <- data.frame(response=rep("species richness"),
                      var=c(names(model.frame(tot_rich_arable3))[2], names(model.frame(tot_rich_broad3))[2], 
                            names(model.frame(tot_rich_conifer2))[2], names(model.frame(tot_rich_grass2))[2]),
                      estimate=c(summary(tot_rich_arable3)$coefficients[2,1], summary(tot_rich_broad3)$coefficients[2,1],
                                 summary(tot_rich_conifer2)$coefficients[2,1], summary(tot_rich_grass2)$coefficients[2,1]), 
                      se=c(summary(tot_rich_arable3)$coefficients[2,2], summary(tot_rich_broad3)$coefficients[2,2],
                           summary(tot_rich_conifer2)$coefficients[2,2], summary(tot_rich_grass2)$coefficients[2,2]), 
                      p_value=c(summary(tot_rich_arable3)$coefficients[2,4], summary(tot_rich_broad3)$coefficients[2,4],
                                summary(tot_rich_conifer2)$coefficients[2,4], summary(tot_rich_grass2)$coefficients[2,4]),
                      r2=c(r.squaredGLMM(tot_rich_arable3)[1,1], r.squaredGLMM(tot_rich_broad3)[1,1], r.squaredGLMM(tot_rich_conifer2)[1,1],
                           r.squaredGLMM(tot_rich_grass2)[1,1]))           
                                 
write.csv(tot_rich_land, file="Results/Bats/Total_richness_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
tot_rich_hab <- glmmTMB(richness ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + 
                        scale(fallen_deadwood) + scale(dis_to_edge) + visit + (1|subcompartment),
                        data = bat_rich, family = "genpois", na.action = "na.fail",
                        control=glmmTMBControl(optimizer=optim,
                                               optArgs=list(method="BFGS")))
summary(tot_rich_hab)
AIC(tot_rich_hab)

pred <- ggpredict(tot_rich_hab, terms = "dbh_average")
richness_dbh <- ggplot(pred, aes(x, predicted)) +
        geom_line() +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
        geom_point(data=bat_rich, aes(x=dbh_average, y=richness))+
        labs(x="Mean DBH", y="Total species \nrichness")+
        scale_y_continuous(breaks=c(2,4,6,8,10))+
        theme_classic()
richness_dbh
ggsave(richness_dbh, file="Graphs/Bats/Total_richness_dbh.png")

final_mod2 <- glmmTMB(richness ~ scale(dbh_average) + scale(basal_area) +
                          scale(per_broadleaf_canopy) + scale(complexity_score) + 
                          scale(canopy_openess) + 
                          scale(fallen_deadwood) + scale(dis_to_edge) + 
                          scale(prop_arable.3000) + visit + (1|subcompartment),
                          data = bat_rich, family = "genpois")
AIC(final_mod2)
summary(final_mod2)

final_mod3 <- glmmTMB(richness ~ scale(dbh_average) + scale(basal_area) +
                          scale(per_broadleaf_canopy) + scale(complexity_score) + 
                          scale(canopy_openess) + 
                          scale(fallen_deadwood) + scale(dis_to_edge) + 
                          scale(prop_broadleaf.3000) + visit + (1|subcompartment),
                          data = bat_rich, family = "genpois", na.action = "na.fail")
AIC(final_mod3)
summary(final_mod3)

final_mod4 <- glmmTMB(richness ~ scale(dbh_average) + scale(basal_area) +
                          scale(per_broadleaf_canopy) + scale(complexity_score) + 
                          scale(canopy_openess) + 
                          scale(fallen_deadwood) + scale(dis_to_edge) + 
                          scale(prop_conifer.1500) + visit + (1|subcompartment),
                          data = bat_rich, family = "genpois", na.action = "na.fail")
AIC(final_mod4)
summary(final_mod4)

final_mod5 <- glmmTMB(richness ~ scale(dbh_average) + scale(basal_area) +
                          scale(per_broadleaf_canopy) + scale(complexity_score) + 
                          scale(canopy_openess) + 
                          scale(fallen_deadwood) + scale(dis_to_edge) + 
                          scale(prop_grassland.1500) + visit + (1|subcompartment),
                          data = bat_rich, family = "genpois", na.action = "na.fail")
AIC(final_mod5)
summary(final_mod5)

# same model selected random vs fixed (no model with AIC lower than 2 values, so habitat model selected)
# same significant habitat variables too

## check model assumptions
testDispersion(final_mod5) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = final_mod5, plot = F)
plot(simulationOutput) ## no assumptions violated!
## Check for multicolinearity
car::vif(tot_rich_hab) ## looking at GIF^(1/(2*Df)) ==> all look good (under 3)
check_collinearity(tot_rich_hab,component = "conditional")

## put all 5 model summaries into one data frame and save
tot_rich_hab_sum <- as.data.frame(summary(tot_rich_hab)$coefficients$cond) #selecting full model coefficient averages
tot_rich_hab_sum$parameters <- row.names(tot_rich_hab_sum)
row.names(tot_rich_hab_sum) <- 1:nrow(tot_rich_hab_sum)
tot_rich_hab_sum$model <- "habitat"
tot_rich_hab_sum$AIC <- AIC(tot_rich_hab)

final_mod2_sum <- as.data.frame(summary(final_mod2)$coefficients$cond) #selecting full model coefficient averages
final_mod2_sum$parameters <- row.names(final_mod2_sum)
row.names(final_mod2_sum) <- 1:nrow(final_mod2_sum)
final_mod2_sum$model <- "arable"
final_mod2_sum$AIC <- AIC(final_mod2)

final_mod3_sum <- as.data.frame(summary(final_mod3)$coefficients$cond) #selecting full model coefficient averages
final_mod3_sum$parameters <- row.names(final_mod3_sum)
row.names(final_mod3_sum) <- 1:nrow(final_mod3_sum)
final_mod3_sum$model <- "broadleaf"
final_mod3_sum$AIC <- AIC(final_mod3)

final_mod4_sum <- as.data.frame(summary(final_mod4)$coefficients$cond) #selecting full model coefficient averages
final_mod4_sum$parameters <- row.names(final_mod4_sum)
row.names(final_mod4_sum) <- 1:nrow(final_mod4_sum)
final_mod4_sum$model <- "conifer"
final_mod4_sum$AIC <- AIC(final_mod4)

final_mod5_sum <- as.data.frame(summary(final_mod5)$coefficients$cond) #selecting full model coefficient averages
final_mod5_sum$parameters <- row.names(final_mod5_sum)
row.names(final_mod5_sum) <- 1:nrow(final_mod5_sum)
final_mod5_sum$model <- "grassland"
final_mod5_sum$AIC <- AIC(final_mod5)

tot_rich_hab_mods <- rbind(tot_rich_hab_sum, final_mod2_sum, final_mod3_sum, final_mod4_sum, final_mod5_sum)
write.csv(tot_rich_hab_mods, file="Results/Bats/Total_richness_habitat.csv", row.names=FALSE)

############################################################################################################################
## Sp1. Nyctalus noctula
noctula <- bat_act[bat_act$scientific_name=="Nyctalus noctula",]
tot_intervals <- noctula$tot_intervals

round(cor(noctula[,c("dbh_average","basal_area","per_broadleaf_canopy","complexity_score","canopy_openess",
                      "deadwood_snags","fallen_deadwood", "dis_to_edge", "prop_grassland.1500")]),3)
cor.test(noctula$per_broadleaf_canopy, noctula$deadwood_snags) ## significant and r>0.8
## include the variable most strongly correlated with response (sp richness)
cor.test(noctula$activity_rate, noctula$per_broadleaf_canopy) # r=-0.0396, p=0.77
cor.test(noctula$activity_rate, noctula$deadwood_snags) # r=-0.03377, p=0.8
## remove deadwood snags from analysis

tot_act_arable1 <- glmer(activity_rate ~ prop_arable.500 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=noctula)
summary(tot_act_arable1)
r.squaredGLMM(tot_act_arable1) # 0.012609211     
tot_act_arable2 <- glmer(activity_rate ~ prop_arable.1500 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=noctula)
summary(tot_act_arable2) 
r.squaredGLMM(tot_act_arable2) # 0.005775566     
tot_act_arable3 <- glmer(activity_rate ~ prop_arable.3000 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=noctula)
summary(tot_act_arable3)    
r.squaredGLMM(tot_act_arable3) # 0.00034868295    
#### ARABLE 500M #### 

tot_act_broad1 <- glmer(activity_rate ~ prop_broadleaf.500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=noctula)
summary(tot_act_broad1) 
r.squaredGLMM(tot_act_broad1) # 0.0008915553      
tot_act_broad2 <- glmer(activity_rate ~ prop_broadleaf.1500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=noctula)
summary(tot_act_broad2)
r.squaredGLMM(tot_act_broad2) # 0.006999837      
tot_act_broad3 <- glmer(activity_rate ~ prop_broadleaf.3000 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=noctula)
summary(tot_act_broad3) 
r.squaredGLMM(tot_act_broad3) # 0.009129797     
#### BROADLEAF 1500M #### (3000m previously)
# 3000M now when visit is fixed variable


tot_act_conifer1 <- glmer(activity_rate ~ prop_conifer.500 + visit + (1|subcompartment),
                           family=binomial, weights=tot_intervals, data=noctula)
summary(tot_act_conifer1)   
r.squaredGLMM(tot_act_conifer1) # 0.000033127905       
tot_act_conifer2 <- glmer(activity_rate ~ prop_conifer.1500 + visit + (1|subcompartment),
                           family=binomial, weights=tot_intervals, data=noctula)
summary(tot_act_conifer2)       
r.squaredGLMM(tot_act_conifer2) # 0.020423808       
tot_act_conifer3 <- glmer(activity_rate ~ prop_conifer.3000 + visit + (1|subcompartment),
                           family=binomial, weights=tot_intervals, data=noctula)
summary(tot_act_conifer3)         
r.squaredGLMM(tot_act_conifer3) # 0.0043533300      
#### CONIFER 1500M #### 

tot_act_grass1 <- glmer(activity_rate ~ prop_grassland.500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=noctula)
summary(tot_act_grass1)        
r.squaredGLMM(tot_act_grass1) # 0.00012180625        
tot_act_grass2 <- glmer(activity_rate ~ prop_grassland.1500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=noctula)
summary(tot_act_grass2)                   
r.squaredGLMM(tot_act_grass2) # 0.007494143        
tot_act_grass3 <- glmer(activity_rate ~ prop_grassland.3000 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=noctula)
summary(tot_act_grass3)          
r.squaredGLMM(tot_act_grass3) # 0.0038252987       
#### GRASSLAND 1500M #### (500m previously)

# save results
tot_act_land <- data.frame(response=rep("Nyctalus noctula"),
                            var=c(names(model.frame(tot_act_arable1))[2], names(model.frame(tot_act_broad2))[2], 
                                  names(model.frame(tot_act_conifer2))[2], names(model.frame(tot_act_grass2))[2]),
                            estimate=c(summary(tot_act_arable1)$coefficients[2,1], summary(tot_act_broad2)$coefficients[2,1],
                                       summary(tot_act_conifer2)$coefficients[2,1], summary(tot_act_grass2)$coefficients[2,1]), 
                            se=c(summary(tot_act_arable1)$coefficients[2,2], summary(tot_act_broad2)$coefficients[2,2],
                                 summary(tot_act_conifer2)$coefficients[2,2], summary(tot_act_grass2)$coefficients[2,2]), 
                            p_value=c(summary(tot_act_arable1)$coefficients[2,4], summary(tot_act_broad2)$coefficients[2,4],
                                      summary(tot_act_conifer2)$coefficients[2,4], summary(tot_act_grass2)$coefficients[2,4]),
                            r2=c(r.squaredGLMM(tot_act_arable1)[1,1], r.squaredGLMM(tot_act_broad2)[1,1], r.squaredGLMM(tot_act_conifer2)[1,1],
                                 r.squaredGLMM(tot_act_grass2)[1,1]))           

write.csv(tot_act_land, file="Results/Bats/Nyctalus_noctula_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
sp_act_hab <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                          scale(per_broadleaf_canopy) + scale(complexity_score) + 
                          scale(canopy_openess) + 
                          scale(fallen_deadwood) + scale(dis_to_edge) + visit + (1|subcompartment),
                        data = noctula, family = binomial, weights=tot_intervals)
summary(sp_act_hab)
AIC(sp_act_hab)

final_mod2 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + 
                        scale(fallen_deadwood) + scale(dis_to_edge) + 
                        scale(prop_arable.500) + visit + (1|subcompartment),
                      data = noctula, family = binomial, weights=tot_intervals)
AIC(final_mod2)

final_mod3 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + 
                        scale(fallen_deadwood) + scale(dis_to_edge) + 
                        scale(prop_broadleaf.3000) + visit + (1|subcompartment),
                      data = noctula, family = binomial, weights=tot_intervals)
AIC(final_mod3)

final_mod4 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) +
                        scale(fallen_deadwood) + scale(dis_to_edge) + 
                        scale(prop_conifer.1500) + visit + (1|subcompartment),
                      data = noctula, family = binomial, weights=tot_intervals)
AIC(final_mod4)
summary(final_mod4)
plot(ggpredict(final_mod4, terms = "prop_conifer.1500"))

final_mod5 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                            scale(per_broadleaf_canopy) + scale(complexity_score) + 
                            scale(canopy_openess) +
                            scale(fallen_deadwood) + scale(dis_to_edge) + 
                            scale(prop_grassland.1500) + visit + (1|subcompartment),
                    data = noctula, family = binomial, weights=tot_intervals)
AIC(final_mod5)
summary(final_mod5)

# same model selected (grassland 1500m) when visit is random vs fixed 
# but conifer 1500m also lower - same as before but we picked grassland
# same significant habitat variables

library(rnrfa)
library(stringr)
lon_lat <- as.data.frame(osg_parse(grid_refs=noctula$gridref, coord_system = "WGS84"))
noctula <- cbind(noctula, lon_lat)
noctula$coords <- paste(noctula$lon,", ",noctula$lat)
coords <- c(unique(noctula$coords))
x_unique <- c(str_extract(coords, "^.+(?=,)"))
y_unique <- c(str_extract(coords, "(?<=, ).+$"))
sims<-simulateResiduals(final_mod4)
simsrecalc<-recalculateResiduals(sims,group = noctula$coords) # recalculate residuals based on unique coordinates
testSpatialAutocorrelation(simsrecalc, x = x_unique, y = y_unique) ## non significant - no evidence of spatial autocorrelation

pred <- ggpredict(final_mod4, terms = "prop_conifer.1500")
noctula_conifer <- ggplot(pred, aes(x, predicted)) +
        geom_line() +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
        geom_point(data=noctula, aes(x=prop_conifer.1500, y=activity_rate))+
        labs(x="Proportion of conifer within 1500m", y="N. noctula \nactivity rate")+
        theme_classic()
noctula_conifer
ggsave(noctula_conifer, file="Graphs/Bats/Noctula_conifer_1500m.png")

final_mod5 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                        scale(per_broadleaf_canopy) + scale(complexity_score) + 
                        scale(canopy_openess) + 
                        scale(fallen_deadwood) + scale(dis_to_edge) + 
                        scale(prop_grassland.1500) + (1|visit) + (1|subcompartment),
                      data = noctula, family = binomial, weights=tot_intervals)
AIC(final_mod5)
summary(final_mod5)
plot(ggpredict(final_mod5, terms = "prop_grassland.1500"))

pred <- ggpredict(final_mod5, terms = "prop_grassland.1500 [all]")
noctula_grass <- ggplot(pred, aes(x, predicted)) +
        geom_line() +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
        geom_point(data=noctula, aes(x=prop_grassland.1500, y=activity_rate))+
        labs(x="Proportion of semi-natural \ngrasssland within 1500m", y="N. noctula \nactivity rate")+
        theme_classic()
noctula_grass
ggsave(noctula_grass, file="Graphs/Bats/Noctula_grassland_1500m.png")


## check model assumptions
testDispersion(sp_act_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = sp_act_hab, plot = F)
plot(simulationOutput) ## no assumptions violated!
## Check for multicolinearity
car::vif(final_mod4) ## looking at GIF^(1/(2*Df)) ==> all look good (under 3)

## put all 5 model summaries into one data frame and save
sp_act_hab_sum <- as.data.frame(coef(summary(sp_act_hab))) #selecting full model coefficient averages
sp_act_hab_sum$parameters <- row.names(sp_act_hab_sum)
row.names(sp_act_hab_sum) <- 1:nrow(sp_act_hab_sum)
sp_act_hab_sum$model <- "habitat"
sp_act_hab_sum$AIC <- AIC(sp_act_hab)

final_mod2_sum <- as.data.frame(coef(summary(final_mod2))) #selecting full model coefficient averages
final_mod2_sum$parameters <- row.names(final_mod2_sum)
row.names(final_mod2_sum) <- 1:nrow(final_mod2_sum)
final_mod2_sum$model <- "arable"
final_mod2_sum$AIC <- AIC(final_mod2)

final_mod3_sum <- as.data.frame(coef(summary(final_mod3))) #selecting full model coefficient averages
final_mod3_sum$parameters <- row.names(final_mod3_sum)
row.names(final_mod3_sum) <- 1:nrow(final_mod3_sum)
final_mod3_sum$model <- "broadleaf"
final_mod3_sum$AIC <- AIC(final_mod3)

final_mod4_sum <- as.data.frame(coef(summary(final_mod4))) #selecting full model coefficient averages
final_mod4_sum$parameters <- row.names(final_mod4_sum)
row.names(final_mod4_sum) <- 1:nrow(final_mod4_sum)
final_mod4_sum$model <- "conifer"
final_mod4_sum$AIC <- AIC(final_mod4)

final_mod5_sum <- as.data.frame(coef(summary(final_mod5))) #selecting full model coefficient averages
final_mod5_sum$parameters <- row.names(final_mod5_sum)
row.names(final_mod5_sum) <- 1:nrow(final_mod5_sum)
final_mod5_sum$model <- "grassland"
final_mod5_sum$AIC <- AIC(final_mod5)

sp_act_hab_mods <- rbind(sp_act_hab_sum, final_mod2_sum, final_mod3_sum, final_mod4_sum, final_mod5_sum)
write.csv(sp_act_hab_mods, file="Results/Bats/Nyctalus_noctula_habitat.csv", row.names=FALSE)

#######################################################################################################################
## Sp2. Plecotus auritus
auritus <- bat_act[bat_act$scientific_name=="Plecotus auritus",]
tot_intervals <- auritus$tot_intervals

round(cor(auritus[,c("dbh_average","basal_area","per_broadleaf_canopy","complexity_score","canopy_openess",
                     "deadwood_snags","fallen_deadwood", "dis_to_edge", "prop_grassland.3000")]),3)
cor.test(auritus$per_broadleaf_canopy, auritus$deadwood_snags) ## significant and r>0.8
## include the variable most strongly correlated with response (sp richness)
cor.test(auritus$activity_rate, auritus$per_broadleaf_canopy) # r=0.1783, p=0.25
cor.test(auritus$activity_rate, auritus$deadwood_snags) # r=0.1122, p=0.47
## remove deadwood snags from analysis

tot_act_arable1 <- glmer(activity_rate ~ prop_arable.500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=auritus)
summary(tot_act_arable1)
r.squaredGLMM(tot_act_arable1) # 0.009987801      
tot_act_arable2 <- glmer(activity_rate ~ prop_arable.1500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=auritus)
summary(tot_act_arable2) 
r.squaredGLMM(tot_act_arable2) # 0.0025652355      
tot_act_arable3 <- glmer(activity_rate ~ prop_arable.3000 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=auritus)
summary(tot_act_arable3)    
r.squaredGLMM(tot_act_arable3) # 0.0045587240     
#### ARABLE 500M #### 

tot_act_broad1 <- glmer(activity_rate ~ prop_broadleaf.500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=auritus)
summary(tot_act_broad1) 
r.squaredGLMM(tot_act_broad1) # 0.0047742281       
tot_act_broad2 <- glmer(activity_rate ~ prop_broadleaf.1500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=auritus)
summary(tot_act_broad2)
r.squaredGLMM(tot_act_broad2) # 0.005598168       
tot_act_broad3 <- glmer(activity_rate ~ prop_broadleaf.3000 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=auritus)
summary(tot_act_broad3) 
r.squaredGLMM(tot_act_broad3) # 0.010089343      
#### BROADLEAF 3000M #### 


tot_act_conifer1 <- glmer(activity_rate ~ prop_conifer.500 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=auritus)
summary(tot_act_conifer1)   
r.squaredGLMM(tot_act_conifer1) # 0.008230157        
tot_act_conifer2 <- glmer(activity_rate ~ prop_conifer.1500 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=auritus)
summary(tot_act_conifer2)       
r.squaredGLMM(tot_act_conifer2) # 0.000053214210        
tot_act_conifer3 <- glmer(activity_rate ~ prop_conifer.3000 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=auritus)
summary(tot_act_conifer3)         
r.squaredGLMM(tot_act_conifer3) # 0.0000029746877       
#### CONIFER 500M #### (1500m previously)

tot_act_grass1 <- glmer(activity_rate ~ prop_grassland.500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=auritus)
summary(tot_act_grass1)        
r.squaredGLMM(tot_act_grass1) # 0.0015524649         
tot_act_grass2 <- glmer(activity_rate ~ prop_grassland.1500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=auritus)
summary(tot_act_grass2)                   
r.squaredGLMM(tot_act_grass2) # 0.00043515734         
tot_act_grass3 <- glmer(activity_rate ~ prop_grassland.3000 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=auritus)
summary(tot_act_grass3)          
r.squaredGLMM(tot_act_grass3) # 0.006864375        
#### GRASSLAND 3000M #### (500m previously)

# same buffer sizes selected random vs fixed

# save results
tot_act_land <- data.frame(response=rep("Plecotus auritus"),
                           var=c(names(model.frame(tot_act_arable1))[2], names(model.frame(tot_act_broad3))[2], 
                                 names(model.frame(tot_act_conifer1))[2], names(model.frame(tot_act_grass3))[2]),
                           estimate=c(summary(tot_act_arable1)$coefficients[2,1], summary(tot_act_broad3)$coefficients[2,1],
                                      summary(tot_act_conifer1)$coefficients[2,1], summary(tot_act_grass3)$coefficients[2,1]), 
                           se=c(summary(tot_act_arable1)$coefficients[2,2], summary(tot_act_broad3)$coefficients[2,2],
                                summary(tot_act_conifer1)$coefficients[2,2], summary(tot_act_grass3)$coefficients[2,2]), 
                           p_value=c(summary(tot_act_arable1)$coefficients[2,4], summary(tot_act_broad3)$coefficients[2,4],
                                     summary(tot_act_conifer1)$coefficients[2,4], summary(tot_act_grass3)$coefficients[2,4]),
                           r2=c(r.squaredGLMM(tot_act_arable1)[1,1], r.squaredGLMM(tot_act_broad3)[1,1], 
                                r.squaredGLMM(tot_act_conifer1)[1,1], r.squaredGLMM(tot_act_grass3)[1,1]))           

write.csv(tot_act_land, file="Results/Bats/Plecotus_auritus_land_cover.csv", row.names=FALSE)


## Now run models with habitat variables only 
sp_act_hab <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + visit + (1|subcompartment),
                    data = auritus, family = binomial, weights=tot_intervals)
summary(sp_act_hab)
AIC(sp_act_hab)

final_mod2 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_arable.500) + visit + (1|subcompartment),
                    data = auritus, family = binomial, weights=tot_intervals)
AIC(final_mod2)

final_mod3 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_broadleaf.3000) + visit + (1|subcompartment),
                    data = auritus, family = binomial, weights=tot_intervals)
AIC(final_mod3)
summary(final_mod3)
plot(ggpredict(final_mod3, terms = "prop_broadleaf.3000"))

final_mod4 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) +  
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_conifer.500) + visit + (1|subcompartment),
                    data = auritus, family = binomial, weights=tot_intervals)
AIC(final_mod4)
summary(final_mod4)
plot(ggpredict(final_mod4, terms = "prop_conifer.1500"))

final_mod5 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_grassland.3000) + visit + (1|subcompartment),
                    data = auritus, family = binomial, weights=tot_intervals)
AIC(final_mod5)

# same model selected random vs fixed (habitat model)
# same significant habitat variables

## put all 5 model summaries into one data frame and save
sp_act_hab_sum <- as.data.frame(coef(summary(sp_act_hab))) #selecting full model coefficient averages
sp_act_hab_sum$parameters <- row.names(sp_act_hab_sum)
row.names(sp_act_hab_sum) <- 1:nrow(sp_act_hab_sum)
sp_act_hab_sum$model <- "habitat"
sp_act_hab_sum$AIC <- AIC(sp_act_hab)

final_mod2_sum <- as.data.frame(coef(summary(final_mod2))) #selecting full model coefficient averages
final_mod2_sum$parameters <- row.names(final_mod2_sum)
row.names(final_mod2_sum) <- 1:nrow(final_mod2_sum)
final_mod2_sum$model <- "arable"
final_mod2_sum$AIC <- AIC(final_mod2)

final_mod3_sum <- as.data.frame(coef(summary(final_mod3))) #selecting full model coefficient averages
final_mod3_sum$parameters <- row.names(final_mod3_sum)
row.names(final_mod3_sum) <- 1:nrow(final_mod3_sum)
final_mod3_sum$model <- "broadleaf"
final_mod3_sum$AIC <- AIC(final_mod3)

final_mod4_sum <- as.data.frame(coef(summary(final_mod4))) #selecting full model coefficient averages
final_mod4_sum$parameters <- row.names(final_mod4_sum)
row.names(final_mod4_sum) <- 1:nrow(final_mod4_sum)
final_mod4_sum$model <- "conifer"
final_mod4_sum$AIC <- AIC(final_mod4)

final_mod5_sum <- as.data.frame(coef(summary(final_mod5))) #selecting full model coefficient averages
final_mod5_sum$parameters <- row.names(final_mod5_sum)
row.names(final_mod5_sum) <- 1:nrow(final_mod5_sum)
final_mod5_sum$model <- "grassland"
final_mod5_sum$AIC <- AIC(final_mod5)

sp_act_hab_mods <- rbind(sp_act_hab_sum, final_mod2_sum, final_mod3_sum, final_mod4_sum, final_mod5_sum)
write.csv(sp_act_hab_mods, file="Results/Bats/Plecotus_auritus_habitat.csv", row.names=FALSE)


## check model assumptions
testDispersion(sp_act_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = sp_act_hab, plot = F)
plot(simulationOutput) ## no assumptions violated!
## Check for multicolinearity
car::vif(sp_act_hab) ## looking at GIF^(1/(2*Df)) ==> all look good (under 3)

################################################################################################################

## Sp3. Myotis nattereri
nattereri <- bat_act[bat_act$scientific_name=="Myotis nattereri",]
tot_intervals <- nattereri$tot_intervals

round(cor(nattereri[,c("dbh_average","basal_area","per_broadleaf_canopy","complexity_score","canopy_openess",
                     "deadwood_snags","fallen_deadwood", "dis_to_edge", "prop_grassland.3000")]),3)
cor.test(nattereri$per_broadleaf_canopy, nattereri$deadwood_snags) ## significant and r>0.8
## include the variable most strongly correlated with response (sp richness)
cor.test(nattereri$activity_rate, nattereri$per_broadleaf_canopy) # r=-0.103, p=0.56
cor.test(nattereri$activity_rate, nattereri$deadwood_snags) # r=-0.228, p=0.19
## remove % broadleaf canopy from analysis

tot_act_arable1 <- glmer(activity_rate ~ prop_arable.500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=nattereri)
summary(tot_act_arable1)
r.squaredGLMM(tot_act_arable1) # 0.0077813497       
tot_act_arable2 <- glmer(activity_rate ~ prop_arable.1500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=nattereri)
summary(tot_act_arable2) 
r.squaredGLMM(tot_act_arable2) # 0.000017139511       
tot_act_arable3 <- glmer(activity_rate ~ prop_arable.3000 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=nattereri)
summary(tot_act_arable3)    
r.squaredGLMM(tot_act_arable3) # 0.0026739779      
#### ARABLE 500M #### (3000m previously)

tot_act_broad1 <- glmer(activity_rate ~ prop_broadleaf.500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=nattereri)
summary(tot_act_broad1) 
r.squaredGLMM(tot_act_broad1) # 0.021919872        
tot_act_broad2 <- glmer(activity_rate ~ prop_broadleaf.1500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=nattereri)
summary(tot_act_broad2)
r.squaredGLMM(tot_act_broad2) # 0.0015728564        
tot_act_broad3 <- glmer(activity_rate ~ prop_broadleaf.3000 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=nattereri)
summary(tot_act_broad3) 
r.squaredGLMM(tot_act_broad3) # 0.0013583374       
#### BROADLEAF 500M #### (1500m previously)


tot_act_conifer1 <- glmer(activity_rate ~ prop_conifer.500 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=nattereri)
summary(tot_act_conifer1)   
r.squaredGLMM(tot_act_conifer1) # 0.004733775         
tot_act_conifer2 <- glmer(activity_rate ~ prop_conifer.1500 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=nattereri)
summary(tot_act_conifer2)       
r.squaredGLMM(tot_act_conifer2) # 0.0026404236         
tot_act_conifer3 <- glmer(activity_rate ~ prop_conifer.3000 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=nattereri)
summary(tot_act_conifer3)         
r.squaredGLMM(tot_act_conifer3) # 0.0020950230        
#### CONIFER 500M #### (1500m previously)

tot_act_grass1 <- glmer(activity_rate ~ prop_grassland.500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=nattereri)
summary(tot_act_grass1)        
r.squaredGLMM(tot_act_grass1) # 0.0018574377          
tot_act_grass2 <- glmer(activity_rate ~ prop_grassland.1500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=nattereri)
summary(tot_act_grass2)                   
r.squaredGLMM(tot_act_grass2) # 0.000070325546          
tot_act_grass3 <- glmer(activity_rate ~ prop_grassland.3000 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=nattereri)
summary(tot_act_grass3)          
r.squaredGLMM(tot_act_grass3) # 0.0024899517         
#### GRASSLAND 3000M #### 

# same buffer sizes selected random vs fixed

# save results
tot_act_land <- data.frame(response=rep("Myotis nattereri"),
                           var=c(names(model.frame(tot_act_arable1))[2], names(model.frame(tot_act_broad1))[2], 
                                 names(model.frame(tot_act_conifer1))[2], names(model.frame(tot_act_grass3))[2]),
                           estimate=c(summary(tot_act_arable1)$coefficients[2,1], summary(tot_act_broad1)$coefficients[2,1],
                                      summary(tot_act_conifer1)$coefficients[2,1], summary(tot_act_grass3)$coefficients[2,1]), 
                           se=c(summary(tot_act_arable1)$coefficients[2,2], summary(tot_act_broad1)$coefficients[2,2],
                                summary(tot_act_conifer1)$coefficients[2,2], summary(tot_act_grass3)$coefficients[2,2]), 
                           p_value=c(summary(tot_act_arable1)$coefficients[2,4], summary(tot_act_broad1)$coefficients[2,4],
                                     summary(tot_act_conifer1)$coefficients[2,4], summary(tot_act_grass3)$coefficients[2,4]),
                           r2=c(r.squaredGLMM(tot_act_arable1)[1,1], r.squaredGLMM(tot_act_broad1)[1,1], 
                                r.squaredGLMM(tot_act_conifer1)[1,1], r.squaredGLMM(tot_act_grass3)[1,1]))           

write.csv(tot_act_land, file="Results/Bats/Myotis_nattereri_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
sp_act_hab <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + visit + (1|subcompartment),
                    data = nattereri, family = binomial, weights=tot_intervals)
summary(sp_act_hab)
AIC(sp_act_hab)

final_mod2 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_arable.500) + visit + (1|subcompartment),
                    data = nattereri, family = binomial, weights=tot_intervals)
AIC(final_mod2)

final_mod3 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_broadleaf.500) + visit + (1|subcompartment),
                    data = nattereri, family = binomial, weights=tot_intervals)
AIC(final_mod3)

final_mod4 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_conifer.500) + visit + (1|subcompartment),
                    data = nattereri, family = binomial, weights=tot_intervals)
AIC(final_mod4)

final_mod5 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_grassland.3000) + visit + (1|subcompartment),
                    data = nattereri, family = binomial, weights=tot_intervals)
AIC(final_mod5)

# same model selected random vs fixed (habitat only)
# same significant habitat variables (none)

## check model assumptions
testDispersion(sp_act_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = sp_act_hab, plot = F)
plot(simulationOutput) ## no assumptions violated!
## Check for multicolinearity
car::vif(final_mod4) ## looking at GIF^(1/(2*Df)) ==> all look good (under 3)

## put all 5 model summaries into one data frame and save
sp_act_hab_sum <- as.data.frame(coef(summary(sp_act_hab))) #selecting full model coefficient averages
sp_act_hab_sum$parameters <- row.names(sp_act_hab_sum)
row.names(sp_act_hab_sum) <- 1:nrow(sp_act_hab_sum)
sp_act_hab_sum$model <- "habitat"
sp_act_hab_sum$AIC <- AIC(sp_act_hab)

final_mod2_sum <- as.data.frame(coef(summary(final_mod2))) #selecting full model coefficient averages
final_mod2_sum$parameters <- row.names(final_mod2_sum)
row.names(final_mod2_sum) <- 1:nrow(final_mod2_sum)
final_mod2_sum$model <- "arable"
final_mod2_sum$AIC <- AIC(final_mod2)

final_mod3_sum <- as.data.frame(coef(summary(final_mod3))) #selecting full model coefficient averages
final_mod3_sum$parameters <- row.names(final_mod3_sum)
row.names(final_mod3_sum) <- 1:nrow(final_mod3_sum)
final_mod3_sum$model <- "broadleaf"
final_mod3_sum$AIC <- AIC(final_mod3)

final_mod4_sum <- as.data.frame(coef(summary(final_mod4))) #selecting full model coefficient averages
final_mod4_sum$parameters <- row.names(final_mod4_sum)
row.names(final_mod4_sum) <- 1:nrow(final_mod4_sum)
final_mod4_sum$model <- "conifer"
final_mod4_sum$AIC <- AIC(final_mod4)

final_mod5_sum <- as.data.frame(coef(summary(final_mod5))) #selecting full model coefficient averages
final_mod5_sum$parameters <- row.names(final_mod5_sum)
row.names(final_mod5_sum) <- 1:nrow(final_mod5_sum)
final_mod5_sum$model <- "grassland"
final_mod5_sum$AIC <- AIC(final_mod5)

sp_act_hab_mods <- rbind(sp_act_hab_sum, final_mod2_sum, final_mod3_sum, final_mod4_sum, final_mod5_sum)
write.csv(sp_act_hab_mods, file="Results/Bats/Myotis_nattereri_habitat.csv", row.names=FALSE)

#################################################################################################################
## Sp4. Myotis mystacinus/brandtii
brandtii <- bat_act[bat_act$scientific_name=="Myotis mystacinus/brandtii",]
tot_intervals <- brandtii$tot_intervals

round(cor(brandtii[,c("dbh_average","basal_area","per_broadleaf_canopy","complexity_score","canopy_openess",
                       "deadwood_snags","fallen_deadwood", "dis_to_edge", "prop_grassland.3000")]),3)
cor.test(brandtii$per_broadleaf_canopy, brandtii$deadwood_snags) ## significant and r>0.8
## include the variable most strongly correlated with response (sp richness)
cor.test(brandtii$activity_rate, brandtii$per_broadleaf_canopy) # r=0.544, p<0.001
cor.test(brandtii$activity_rate, brandtii$deadwood_snags) # r=0.32, p=0.034
## remove deadwood snags from analysis

tot_act_arable1 <- glmer(activity_rate ~ prop_arable.500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=brandtii)
summary(tot_act_arable1)
r.squaredGLMM(tot_act_arable1) # 0.019636323        
tot_act_arable2 <- glmer(activity_rate ~ prop_arable.1500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=brandtii)
summary(tot_act_arable2) 
r.squaredGLMM(tot_act_arable2) # 0.05463143        
tot_act_arable3 <- glmer(activity_rate ~ prop_arable.3000 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=brandtii)
summary(tot_act_arable3)    
r.squaredGLMM(tot_act_arable3) # 0.015861915       
#### ARABLE 1500M #### 

tot_act_broad1 <- glmer(activity_rate ~ prop_broadleaf.500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=brandtii)
summary(tot_act_broad1) 
r.squaredGLMM(tot_act_broad1) # 0.04287694        
tot_act_broad2 <- glmer(activity_rate ~ prop_broadleaf.1500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=brandtii)
summary(tot_act_broad2)
r.squaredGLMM(tot_act_broad2) # 0.04136314         
tot_act_broad3 <- glmer(activity_rate ~ prop_broadleaf.3000 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=brandtii)
summary(tot_act_broad3) 
r.squaredGLMM(tot_act_broad3) # 0.05188494        
#### BROADLEAF 3000M #### (500m previously) 


tot_act_conifer1 <- glmer(activity_rate ~ prop_conifer.500 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=brandtii)
summary(tot_act_conifer1)   
r.squaredGLMM(tot_act_conifer1) # 0.05732421          
tot_act_conifer2 <- glmer(activity_rate ~ prop_conifer.1500 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=brandtii)
summary(tot_act_conifer2)       
r.squaredGLMM(tot_act_conifer2) # 0.02754497          
tot_act_conifer3 <- glmer(activity_rate ~ prop_conifer.3000 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=brandtii)
summary(tot_act_conifer3)         
r.squaredGLMM(tot_act_conifer3) # 0.08751967         
#### CONIFER 3000M #### 

tot_act_grass1 <- glmer(activity_rate ~ prop_grassland.500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=brandtii)
summary(tot_act_grass1)        
r.squaredGLMM(tot_act_grass1) # 0.008179273           
tot_act_grass2 <- glmer(activity_rate ~ prop_grassland.1500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=brandtii)
summary(tot_act_grass2)                   
r.squaredGLMM(tot_act_grass2) # 0.032429610           
tot_act_grass3 <- glmer(activity_rate ~ prop_grassland.3000 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=brandtii)
summary(tot_act_grass3)          
r.squaredGLMM(tot_act_grass3) # 0.019279408          
#### GRASSLAND 1500M #### 

# same buffer sizes selected random vs fixed

# save results
tot_act_land <- data.frame(response=rep("Myotis mystacinus/brandtii"),
                           var=c(names(model.frame(tot_act_arable2))[2], names(model.frame(tot_act_broad3))[2], 
                                 names(model.frame(tot_act_conifer3))[2], names(model.frame(tot_act_grass2))[2]),
                           estimate=c(summary(tot_act_arable2)$coefficients[2,1], summary(tot_act_broad3)$coefficients[2,1],
                                      summary(tot_act_conifer3)$coefficients[2,1], summary(tot_act_grass2)$coefficients[2,1]), 
                           se=c(summary(tot_act_arable2)$coefficients[2,2], summary(tot_act_broad3)$coefficients[2,2],
                                summary(tot_act_conifer3)$coefficients[2,2], summary(tot_act_grass2)$coefficients[2,2]), 
                           p_value=c(summary(tot_act_arable2)$coefficients[2,4], summary(tot_act_broad3)$coefficients[2,4],
                                     summary(tot_act_conifer3)$coefficients[2,4], summary(tot_act_grass2)$coefficients[2,4]),
                           r2=c(r.squaredGLMM(tot_act_arable2)[1,1], r.squaredGLMM(tot_act_broad3)[1,1], 
                                r.squaredGLMM(tot_act_conifer3)[1,1], r.squaredGLMM(tot_act_grass2)[1,1]))           

write.csv(tot_act_land, file="Results/Bats/Myotis_mystacinus_brandtii_land_cover.csv", row.names=FALSE)


## Now run models with habitat variables only 
sp_act_hab <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + visit + (1|subcompartment),
                    data = brandtii, family = binomial, weights=tot_intervals)
summary(sp_act_hab)
AIC(sp_act_hab)

pred <- ggpredict(sp_act_hab, terms = "per_broadleaf_canopy")
brandtii_per_broadleaf <- ggplot(pred, aes(x, predicted)) +
        geom_line() +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
        geom_point(data=brandtii, aes(x=per_broadleaf_canopy, y=activity_rate))+
        labs(x="Percentage \nbroadleaf canopy", y="M. mystacinus \n/brandtii activity \nrate")+
        #scale_y_continuous(breaks=c(2,4,6,8,10))+
        theme_classic()
brandtii_per_broadleaf
ggsave(brandtii_per_broadleaf, file="Graphs/Bats/Brandtii_per_broadleaf.png")

final_mod2 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_arable.1500) + visit + (1|subcompartment),
                    data = brandtii, family = binomial, weights=tot_intervals)
AIC(final_mod2)

final_mod3 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_broadleaf.3000) + visit + (1|subcompartment),
                    data = brandtii, family = binomial, weights=tot_intervals)
AIC(final_mod3)

final_mod4 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_conifer.3000) + visit + (1|subcompartment),
                    data = brandtii, family = binomial, weights=tot_intervals)
AIC(final_mod4)

final_mod5 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_grassland.1500) + visit + (1|subcompartment),
                    data = brandtii, family = binomial, weights=tot_intervals)
AIC(final_mod5)

# same model selected (habitat only) random vs fixed
# same significant habitat variables

## check model assumptions
testDispersion(sp_act_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = sp_act_hab, plot = F)
plot(simulationOutput) ## no assumptions violated!
## Check for multicolinearity
car::vif(final_mod4) ## looking at GIF^(1/(2*Df)) ==> all look good (under 3)

## put all 5 model summaries into one data frame and save
sp_act_hab_sum <- as.data.frame(coef(summary(sp_act_hab))) #selecting full model coefficient averages
sp_act_hab_sum$parameters <- row.names(sp_act_hab_sum)
row.names(sp_act_hab_sum) <- 1:nrow(sp_act_hab_sum)
sp_act_hab_sum$model <- "habitat"
sp_act_hab_sum$AIC <- AIC(sp_act_hab)

final_mod2_sum <- as.data.frame(coef(summary(final_mod2))) #selecting full model coefficient averages
final_mod2_sum$parameters <- row.names(final_mod2_sum)
row.names(final_mod2_sum) <- 1:nrow(final_mod2_sum)
final_mod2_sum$model <- "arable"
final_mod2_sum$AIC <- AIC(final_mod2)

final_mod3_sum <- as.data.frame(coef(summary(final_mod3))) #selecting full model coefficient averages
final_mod3_sum$parameters <- row.names(final_mod3_sum)
row.names(final_mod3_sum) <- 1:nrow(final_mod3_sum)
final_mod3_sum$model <- "broadleaf"
final_mod3_sum$AIC <- AIC(final_mod3)

final_mod4_sum <- as.data.frame(coef(summary(final_mod4))) #selecting full model coefficient averages
final_mod4_sum$parameters <- row.names(final_mod4_sum)
row.names(final_mod4_sum) <- 1:nrow(final_mod4_sum)
final_mod4_sum$model <- "conifer"
final_mod4_sum$AIC <- AIC(final_mod4)

final_mod5_sum <- as.data.frame(coef(summary(final_mod5))) #selecting full model coefficient averages
final_mod5_sum$parameters <- row.names(final_mod5_sum)
row.names(final_mod5_sum) <- 1:nrow(final_mod5_sum)
final_mod5_sum$model <- "grassland"
final_mod5_sum$AIC <- AIC(final_mod5)

sp_act_hab_mods <- rbind(sp_act_hab_sum, final_mod2_sum, final_mod3_sum, final_mod4_sum, final_mod5_sum)
write.csv(sp_act_hab_mods, file="Results/Bats/Myotis_mystacinus_brandtii_habitat.csv", row.names=FALSE)

############################################################################################################################
## Sp5. Pipistrellus pipistrellus
pipistrellus <- bat_act[bat_act$scientific_name=="Pipistrellus pipistrellus",]
tot_intervals <- pipistrellus$tot_intervals

round(cor(pipistrellus[,c("dbh_average","basal_area","per_broadleaf_canopy","complexity_score","canopy_openess",
                      "deadwood_snags","fallen_deadwood", "dis_to_edge", "prop_arable.1500")]),3)
cor.test(pipistrellus$per_broadleaf_canopy, pipistrellus$deadwood_snags) ## significant and r>0.8
## include the variable most strongly correlated with response (sp richness)
cor.test(pipistrellus$activity_rate, pipistrellus$per_broadleaf_canopy) 
cor.test(pipistrellus$activity_rate, pipistrellus$deadwood_snags) 
## remove % broadleaf canopy from analysis

tot_act_arable1 <- glmer(activity_rate ~ prop_arable.500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=pipistrellus)
summary(tot_act_arable1)
r.squaredGLMM(tot_act_arable1) # 0.0011111648         
tot_act_arable2 <- glmer(activity_rate ~ prop_arable.1500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=pipistrellus)
summary(tot_act_arable2) 
r.squaredGLMM(tot_act_arable2) # 0.03751759         
tot_act_arable3 <- glmer(activity_rate ~ prop_arable.3000 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=pipistrellus)
summary(tot_act_arable3)    
r.squaredGLMM(tot_act_arable3) # 0.00002733767        
#### ARABLE 1500M #### 

tot_act_broad1 <- glmer(activity_rate ~ prop_broadleaf.500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=pipistrellus)
summary(tot_act_broad1) 
r.squaredGLMM(tot_act_broad1) # 0.013676790         
tot_act_broad2 <- glmer(activity_rate ~ prop_broadleaf.1500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=pipistrellus)
summary(tot_act_broad2)
r.squaredGLMM(tot_act_broad2) # 0.009447555          
tot_act_broad3 <- glmer(activity_rate ~ prop_broadleaf.3000 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=pipistrellus)
summary(tot_act_broad3) 
r.squaredGLMM(tot_act_broad3) # 0.03349531         
#### BROADLEAF 3000M #### (500m previously)


tot_act_conifer1 <- glmer(activity_rate ~ prop_conifer.500 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=pipistrellus)
summary(tot_act_conifer1)   
r.squaredGLMM(tot_act_conifer1) # 0.002556930           
tot_act_conifer2 <- glmer(activity_rate ~ prop_conifer.1500 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=pipistrellus)
summary(tot_act_conifer2)       
r.squaredGLMM(tot_act_conifer2) # 0.03242995          
tot_act_conifer3 <- glmer(activity_rate ~ prop_conifer.3000 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=pipistrellus)
summary(tot_act_conifer3)         
r.squaredGLMM(tot_act_conifer3) # 0.0013140366          
#### CONIFER 1500M #### 

tot_act_grass1 <- glmer(activity_rate ~ prop_grassland.500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=pipistrellus)
summary(tot_act_grass1)        
r.squaredGLMM(tot_act_grass1) # 0.03542227            
tot_act_grass2 <- glmer(activity_rate ~ prop_grassland.1500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=pipistrellus)
summary(tot_act_grass2)                   
r.squaredGLMM(tot_act_grass2) # 0.04571914            
tot_act_grass3 <- glmer(activity_rate ~ prop_grassland.3000 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=pipistrellus)
summary(tot_act_grass3)          
r.squaredGLMM(tot_act_grass3) # 0.006147952           
#### GRASSLAND 1500M #### 

# same buffer sizes selected random vs fixed 

# save results
tot_act_land <- data.frame(response=rep("Pipistrellus pipistrellus"),
                           var=c(names(model.frame(tot_act_arable2))[2], names(model.frame(tot_act_broad3))[2], 
                                 names(model.frame(tot_act_conifer2))[2], names(model.frame(tot_act_grass2))[2]),
                           estimate=c(summary(tot_act_arable2)$coefficients[2,1], summary(tot_act_broad3)$coefficients[2,1],
                                      summary(tot_act_conifer2)$coefficients[2,1], summary(tot_act_grass2)$coefficients[2,1]), 
                           se=c(summary(tot_act_arable2)$coefficients[2,2], summary(tot_act_broad3)$coefficients[2,2],
                                summary(tot_act_conifer2)$coefficients[2,2], summary(tot_act_grass2)$coefficients[2,2]), 
                           p_value=c(summary(tot_act_arable2)$coefficients[2,4], summary(tot_act_broad3)$coefficients[2,4],
                                     summary(tot_act_conifer2)$coefficients[2,4], summary(tot_act_grass2)$coefficients[2,4]),
                           r2=c(r.squaredGLMM(tot_act_arable2)[1,1], r.squaredGLMM(tot_act_broad3)[1,1], 
                                r.squaredGLMM(tot_act_conifer2)[1,1], r.squaredGLMM(tot_act_grass2)[1,1]))           

write.csv(tot_act_land, file="Results/Bats/Pipistrellus_pipistrellus_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
sp_act_hab <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + visit + (1|subcompartment),
                    data = pipistrellus, family = binomial, weights=tot_intervals)
summary(sp_act_hab)
AIC(sp_act_hab)

final_mod2 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_arable.1500) + visit + (1|subcompartment),
                    data = pipistrellus, family = binomial, weights=tot_intervals)
AIC(final_mod2)

final_mod3 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_broadleaf.3000) + visit + (1|subcompartment),
                    data = pipistrellus, family = binomial, weights=tot_intervals)
AIC(final_mod3)
summary(final_mod3)

final_mod4 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_conifer.1500) + visit + (1|subcompartment),
                    data = pipistrellus, family = binomial, weights=tot_intervals)
AIC(final_mod4)
summary(final_mod4)
plot(ggpredict(final_mod4, terms = "prop_conifer.1500"))

pred <- ggpredict(final_mod3, terms = "canopy_openess [all]")
pipistrellus_openness <- ggplot(pred, aes(x, predicted)) +
        geom_line() +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
        geom_point(data=pipistrellus, aes(x=canopy_openess, y=activity_rate))+
        labs(x="Canopy openness", y="P. pipistrellus \nactivity rate")+
        theme_classic()
pipistrellus_openness
ggsave(pipistrellus_openness, file="Graphs/Bats/Pipistrellus_canopy_openness.png")

pred <- ggpredict(final_mod3, terms = "deadwood_snags [all]")
pipistrellus_snags <- ggplot(pred, aes(x, predicted)) +
        geom_line() +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
        geom_point(data=pipistrellus, aes(x=deadwood_snags, y=activity_rate))+
        labs(x="Deadwood snags", y="P. pipistrellus \nactivity rate")+
        theme_classic()
pipistrellus_snags
ggsave(pipistrellus_snags, file="Graphs/Bats/Pipistrellus_deadwood_snags.png")

pred <- ggpredict(final_mod3, terms = "fallen_deadwood [all]")
pipistrellus_fallen <- ggplot(pred, aes(x, predicted)) +
        geom_line() +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
        geom_point(data=pipistrellus, aes(x=fallen_deadwood, y=activity_rate))+
        labs(x="Fallen deadwood", y="P. pipistrellus \nactivity rate")+
        theme_classic()
pipistrellus_fallen
ggsave(pipistrellus_fallen, file="Graphs/Bats/Pipistrellus_fallen_deadwood.png")

pred <- ggpredict(final_mod3, terms = "complexity_score [all]")
pipistrellus_complex <- ggplot(pred, aes(x, predicted)) +
        geom_line() +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
        geom_point(data=pipistrellus, aes(x=complexity_score, y=activity_rate))+
        labs(x="Complexity score", y="P. pipistrellus \nactivity rate")+
        theme_classic()
pipistrellus_complex
ggsave(pipistrellus_complex, file="Graphs/Bats/Pipistrellus_complexity_score.png")


final_mod5 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_grassland.1500) + visit + (1|subcompartment),
                    data = pipistrellus, family = binomial, weights=tot_intervals)
AIC(final_mod5)

## check model assumptions
testDispersion(sp_act_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = sp_act_hab, plot = F)
plot(simulationOutput) ## no assumptions violated!
## Check for multicolinearity
car::vif(final_mod4) ## looking at GIF^(1/(2*Df)) ==> all look good (under 3)

## put all 5 model summaries into one data frame and save
sp_act_hab_sum <- as.data.frame(coef(summary(sp_act_hab))) #selecting full model coefficient averages
sp_act_hab_sum$parameters <- row.names(sp_act_hab_sum)
row.names(sp_act_hab_sum) <- 1:nrow(sp_act_hab_sum)
sp_act_hab_sum$model <- "habitat"
sp_act_hab_sum$AIC <- AIC(sp_act_hab)

final_mod2_sum <- as.data.frame(coef(summary(final_mod2))) #selecting full model coefficient averages
final_mod2_sum$parameters <- row.names(final_mod2_sum)
row.names(final_mod2_sum) <- 1:nrow(final_mod2_sum)
final_mod2_sum$model <- "arable"
final_mod2_sum$AIC <- AIC(final_mod2)

final_mod3_sum <- as.data.frame(coef(summary(final_mod3))) #selecting full model coefficient averages
final_mod3_sum$parameters <- row.names(final_mod3_sum)
row.names(final_mod3_sum) <- 1:nrow(final_mod3_sum)
final_mod3_sum$model <- "broadleaf"
final_mod3_sum$AIC <- AIC(final_mod3)

final_mod4_sum <- as.data.frame(coef(summary(final_mod4))) #selecting full model coefficient averages
final_mod4_sum$parameters <- row.names(final_mod4_sum)
row.names(final_mod4_sum) <- 1:nrow(final_mod4_sum)
final_mod4_sum$model <- "conifer"
final_mod4_sum$AIC <- AIC(final_mod4)

final_mod5_sum <- as.data.frame(coef(summary(final_mod5))) #selecting full model coefficient averages
final_mod5_sum$parameters <- row.names(final_mod5_sum)
row.names(final_mod5_sum) <- 1:nrow(final_mod5_sum)
final_mod5_sum$model <- "grassland"
final_mod5_sum$AIC <- AIC(final_mod5)

sp_act_hab_mods <- rbind(sp_act_hab_sum, final_mod2_sum, final_mod3_sum, final_mod4_sum, final_mod5_sum)
write.csv(sp_act_hab_mods, file="Results/Bats/Pipistrellus_pipistrellus_habitat.csv", row.names=FALSE)

############################################################################################################
## Sp6. Pipistrellus pygmaeus
pygmaeus <- bat_act[bat_act$scientific_name=="Pipistrellus pygmaeus",]
tot_intervals <- pygmaeus$tot_intervals

round(cor(pygmaeus[,c("dbh_average","basal_area","per_broadleaf_canopy","complexity_score","canopy_openess",
                           "deadwood_snags","fallen_deadwood", "dis_to_edge", "prop_grassland.500")]),3)
cor.test(pygmaeus$per_broadleaf_canopy, pygmaeus$deadwood_snags) ## significant and r>0.8
## include the variable most strongly correlated with response (sp richness)
cor.test(pygmaeus$activity_rate, pygmaeus$per_broadleaf_canopy) # r=0.3384, p=0.0071
cor.test(pygmaeus$activity_rate, pygmaeus$deadwood_snags) # r=0.433, p<0.001
## remove % broadleaf canopy from analysis

tot_act_arable1 <- glmer(activity_rate ~ prop_arable.500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=pygmaeus)
summary(tot_act_arable1)
r.squaredGLMM(tot_act_arable1) # 0.003162026          
tot_act_arable2 <- glmer(activity_rate ~ prop_arable.1500 + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=pygmaeus)
summary(tot_act_arable2) 
r.squaredGLMM(tot_act_arable2) # 0.00008031302          
tot_act_arable3 <- glmer(activity_rate ~ prop_arable.3000 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=pygmaeus)
summary(tot_act_arable3)    
r.squaredGLMM(tot_act_arable3) # 0.009573085        
#### ARABLE 3000M #### (500m previously) 

tot_act_broad1 <- glmer(activity_rate ~ prop_broadleaf.500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=pygmaeus)
summary(tot_act_broad1) 
r.squaredGLMM(tot_act_broad1) # 0.01997744          
tot_act_broad2 <- glmer(activity_rate ~ prop_broadleaf.1500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=pygmaeus)
summary(tot_act_broad2)
r.squaredGLMM(tot_act_broad2) # 0.013685592           
tot_act_broad3 <- glmer(activity_rate ~ prop_broadleaf.3000 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=pygmaeus)
summary(tot_act_broad3) 
r.squaredGLMM(tot_act_broad3) # 0.013702968          
#### BROADLEAF 500M #### 


tot_act_conifer1 <- glmer(activity_rate ~ prop_conifer.500 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=pygmaeus)
summary(tot_act_conifer1)   
r.squaredGLMM(tot_act_conifer1) # 0.002969737            
tot_act_conifer2 <- glmer(activity_rate ~ prop_conifer.1500 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=pygmaeus)
summary(tot_act_conifer2)       
r.squaredGLMM(tot_act_conifer2) # 0.004672719           
tot_act_conifer3 <- glmer(activity_rate ~ prop_conifer.3000 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=pygmaeus)
summary(tot_act_conifer3)         
r.squaredGLMM(tot_act_conifer3) # 0.007447870           
#### CONIFER 3000M #### (500m previously)

tot_act_grass1 <- glmer(activity_rate ~ prop_grassland.500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=pygmaeus)
summary(tot_act_grass1)        
r.squaredGLMM(tot_act_grass1) # 0.05533833             
tot_act_grass2 <- glmer(activity_rate ~ prop_grassland.1500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=pygmaeus)
summary(tot_act_grass2)                   
r.squaredGLMM(tot_act_grass2) # 0.00019908973             
tot_act_grass3 <- glmer(activity_rate ~ prop_grassland.3000 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=pygmaeus)
summary(tot_act_grass3)          
r.squaredGLMM(tot_act_grass3) # 0.016515591            
#### GRASSLAND 500M #### 

# same buffer sizes selected random vs fixed

# save results
tot_act_land <- data.frame(response=rep("Pipistrellus pygmaeus"),
                           var=c(names(model.frame(tot_act_arable3))[2], names(model.frame(tot_act_broad1))[2], 
                                 names(model.frame(tot_act_conifer3))[2], names(model.frame(tot_act_grass1))[2]),
                           estimate=c(summary(tot_act_arable3)$coefficients[2,1], summary(tot_act_broad1)$coefficients[2,1],
                                      summary(tot_act_conifer3)$coefficients[2,1], summary(tot_act_grass1)$coefficients[2,1]), 
                           se=c(summary(tot_act_arable3)$coefficients[2,2], summary(tot_act_broad1)$coefficients[2,2],
                                summary(tot_act_conifer3)$coefficients[2,2], summary(tot_act_grass1)$coefficients[2,2]), 
                           p_value=c(summary(tot_act_arable3)$coefficients[2,4], summary(tot_act_broad1)$coefficients[2,4],
                                     summary(tot_act_conifer3)$coefficients[2,4], summary(tot_act_grass1)$coefficients[2,4]),
                           r2=c(r.squaredGLMM(tot_act_arable3)[1,1], r.squaredGLMM(tot_act_broad1)[1,1], 
                                r.squaredGLMM(tot_act_conifer3)[1,1], r.squaredGLMM(tot_act_grass1)[1,1]))           

write.csv(tot_act_land, file="Results/Bats/Pipistrellus_pygmaeus_land_cover.csv", row.names=FALSE)


## Now run models with habitat variables only 
sp_act_hab <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + visit + (1|subcompartment),
                    data = pygmaeus, family = binomial, weights=tot_intervals)
summary(sp_act_hab)
AIC(sp_act_hab)

final_mod2 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_arable.3000) + visit + (1|subcompartment),
                    data = pygmaeus, family = binomial, weights=tot_intervals)
AIC(final_mod2)

final_mod3 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_broadleaf.500) + visit + (1|subcompartment),
                    data = pygmaeus, family = binomial, weights=tot_intervals)
AIC(final_mod3)
summary(final_mod3)
plot(ggpredict(final_mod3, terms = "prop_broadleaf.500"))

pred <- ggpredict(final_mod3, terms = "prop_broadleaf.500 [all]")
pygmaeus_broadleaf <- ggplot(pred, aes(x, predicted)) +
        geom_line() +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
        geom_point(data=pygmaeus, aes(x=prop_broadleaf.500, y=activity_rate))+
        labs(x="Proportion of broadleaf within 500m", y="P. pygmaeus \nactivity rate")+
        theme_classic()
pygmaeus_broadleaf
ggsave(pygmaeus_broadleaf, file="Graphs/Bats/Pygmaeus_broadleaf_500m.png")

final_mod4 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_conifer.3000) + visit + (1|subcompartment),
                    data = pygmaeus, family = binomial, weights=tot_intervals)
AIC(final_mod4)

final_mod5 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_grassland.500) + visit + (1|subcompartment),
                    data = pygmaeus, family = binomial, weights=tot_intervals)
AIC(final_mod5)
summary(final_mod5)

pred <- ggpredict(final_mod5, terms = "dbh_average [all]")
pygmaeus_dbh <- ggplot(pred, aes(x, predicted)) +
        geom_line() +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
        geom_point(data=pygmaeus, aes(x=dbh_average, y=activity_rate))+
        labs(x="Mean DBH", y="P. pygmaeus \nactivity rate")+
        theme_classic()
pygmaeus_dbh
ggsave(pygmaeus_dbh, file="Graphs/Bats/Pygmaeus_mean_dbh.png")

pred <- ggpredict(final_mod5, terms = "canopy_openess [all]")
pygmaeus_open <- ggplot(pred, aes(x, predicted)) +
        geom_line() +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
        geom_point(data=pygmaeus, aes(x=canopy_openess, y=activity_rate))+
        labs(x="Canopy openness", y="P. pygmaeus \nactivity rate")+
        theme_classic()
pygmaeus_open
ggsave(pygmaeus_open, file="Graphs/Bats/Pygmaeus_canopy_openness.png")

pred <- ggpredict(final_mod5, terms = "deadwood_snags [all]")
pygmaeus_snags <- ggplot(pred, aes(x, predicted)) +
        geom_line() +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
        geom_point(data=pygmaeus, aes(x=deadwood_snags, y=activity_rate))+
        labs(x="Deadwood snags", y="P. pygmaeus \nactivity rate")+
        theme_classic()
pygmaeus_snags
ggsave(pygmaeus_snags, file="Graphs/Bats/Pygmaeus_deadwood_snags.png")

pred <- ggpredict(final_mod5, terms = "fallen_deadwood [all]")
pygmaeus_fallen <- ggplot(pred, aes(x, predicted)) +
        geom_line() +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
        geom_point(data=pygmaeus, aes(x=fallen_deadwood, y=activity_rate))+
        labs(x="Fallen deadwood", y="P. pygmaeus \nactivity rate")+
        theme_classic()
pygmaeus_fallen
ggsave(pygmaeus_fallen, file="Graphs/Bats/Pygmaeus_mean_fallen_deadwood.png")

# same model selected (grassland) random vs fixed
# same habitat variables significant 

## check model assumptions
testDispersion(sp_act_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = sp_act_hab, plot = F)
plot(simulationOutput) ## no assumptions violated!
## Check for multicolinearity
car::vif(final_mod4) ## looking at GIF^(1/(2*Df)) ==> all look good (under 3)

## put all 5 model summaries into one data frame and save
sp_act_hab_sum <- as.data.frame(coef(summary(sp_act_hab))) #selecting full model coefficient averages
sp_act_hab_sum$parameters <- row.names(sp_act_hab_sum)
row.names(sp_act_hab_sum) <- 1:nrow(sp_act_hab_sum)
sp_act_hab_sum$model <- "habitat"
sp_act_hab_sum$AIC <- AIC(sp_act_hab)

final_mod2_sum <- as.data.frame(coef(summary(final_mod2))) #selecting full model coefficient averages
final_mod2_sum$parameters <- row.names(final_mod2_sum)
row.names(final_mod2_sum) <- 1:nrow(final_mod2_sum)
final_mod2_sum$model <- "arable"
final_mod2_sum$AIC <- AIC(final_mod2)

final_mod3_sum <- as.data.frame(coef(summary(final_mod3))) #selecting full model coefficient averages
final_mod3_sum$parameters <- row.names(final_mod3_sum)
row.names(final_mod3_sum) <- 1:nrow(final_mod3_sum)
final_mod3_sum$model <- "broadleaf"
final_mod3_sum$AIC <- AIC(final_mod3)

final_mod4_sum <- as.data.frame(coef(summary(final_mod4))) #selecting full model coefficient averages
final_mod4_sum$parameters <- row.names(final_mod4_sum)
row.names(final_mod4_sum) <- 1:nrow(final_mod4_sum)
final_mod4_sum$model <- "conifer"
final_mod4_sum$AIC <- AIC(final_mod4)

final_mod5_sum <- as.data.frame(coef(summary(final_mod5))) #selecting full model coefficient averages
final_mod5_sum$parameters <- row.names(final_mod5_sum)
row.names(final_mod5_sum) <- 1:nrow(final_mod5_sum)
final_mod5_sum$model <- "grassland"
final_mod5_sum$AIC <- AIC(final_mod5)

sp_act_hab_mods <- rbind(sp_act_hab_sum, final_mod2_sum, final_mod3_sum, final_mod4_sum, final_mod5_sum)
write.csv(sp_act_hab_mods, file="Results/Bats/Pipistrellus_pygmaeus_habitat.csv", row.names=FALSE)

#########################################################################################################
## Sp7. Barbastella barbastellus
barbastellus <- bat_act[bat_act$scientific_name=="Barbastella barbastellus",]
tot_intervals <- barbastellus$tot_intervals

round(cor(barbastellus[,c("dbh_average","basal_area","per_broadleaf_canopy","complexity_score","canopy_openess",
                      "deadwood_snags","fallen_deadwood", "dis_to_edge", "prop_conifer.500")]),3)
cor.test(barbastellus$per_broadleaf_canopy, barbastellus$deadwood_snags) ## significant and r>0.8
## include the variable most strongly correlated with response (sp richness)
cor.test(barbastellus$activity_rate, barbastellus$per_broadleaf_canopy)
cor.test(barbastellus$activity_rate, barbastellus$deadwood_snags) 
## remove % broadleaf canopy from analysis

tot_act_arable1 <- glmer(activity_rate ~ prop_arable.500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=barbastellus)
summary(tot_act_arable1)
r.squaredGLMM(tot_act_arable1) # 0.00039527800           
tot_act_arable2 <- glmer(activity_rate ~ prop_arable.1500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=barbastellus)
summary(tot_act_arable2) 
r.squaredGLMM(tot_act_arable2) # 0.006048960           
tot_act_arable3 <- glmer(activity_rate ~ prop_arable.3000 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=barbastellus)
summary(tot_act_arable3)    
r.squaredGLMM(tot_act_arable3) # 0.0017600405         
#### ARABLE 3000M #### (500m previously)
# 1500m with visit as fixed effect

tot_act_broad1 <- glmer(activity_rate ~ prop_broadleaf.500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=barbastellus)
summary(tot_act_broad1) 
r.squaredGLMM(tot_act_broad1) # 0.017485263           
tot_act_broad2 <- glmer(activity_rate ~ prop_broadleaf.1500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=barbastellus)
summary(tot_act_broad2)
r.squaredGLMM(tot_act_broad2) # 0.002081460            
tot_act_broad3 <- glmer(activity_rate ~ prop_broadleaf.3000 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=barbastellus)
summary(tot_act_broad3) 
r.squaredGLMM(tot_act_broad3) # 0.0020257022           
#### BROADLEAF 500M #### 


tot_act_conifer1 <- glmer(activity_rate ~ prop_conifer.500 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=barbastellus)
summary(tot_act_conifer1)   
r.squaredGLMM(tot_act_conifer1) # 0.033629389             
tot_act_conifer2 <- glmer(activity_rate ~ prop_conifer.1500 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=barbastellus)
summary(tot_act_conifer2)       
r.squaredGLMM(tot_act_conifer2) # 0.006865320            
tot_act_conifer3 <- glmer(activity_rate ~ prop_conifer.3000 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=barbastellus)
summary(tot_act_conifer3)         
r.squaredGLMM(tot_act_conifer3) # 0.0025175748            
#### CONIFER 500M #### 

tot_act_grass1 <- glmer(activity_rate ~ prop_grassland.500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=barbastellus)
summary(tot_act_grass1)        
r.squaredGLMM(tot_act_grass1) # 0.0013298201              
tot_act_grass2 <- glmer(activity_rate ~ prop_grassland.1500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=barbastellus)
summary(tot_act_grass2)                   
r.squaredGLMM(tot_act_grass2) # 0.034135776              
tot_act_grass3 <- glmer(activity_rate ~ prop_grassland.3000 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=barbastellus)
summary(tot_act_grass3)          
r.squaredGLMM(tot_act_grass3) # 0.000014157098             
#### GRASSLAND 1500M #### 

# save results
tot_act_land <- data.frame(response=rep("Barbastella barbastellus"),
                           var=c(names(model.frame(tot_act_arable3))[2], names(model.frame(tot_act_broad1))[2], 
                                 names(model.frame(tot_act_conifer1))[2], names(model.frame(tot_act_grass2))[2]),
                           estimate=c(summary(tot_act_arable3)$coefficients[2,1], summary(tot_act_broad1)$coefficients[2,1],
                                      summary(tot_act_conifer1)$coefficients[2,1], summary(tot_act_grass2)$coefficients[2,1]), 
                           se=c(summary(tot_act_arable3)$coefficients[2,2], summary(tot_act_broad1)$coefficients[2,2],
                                summary(tot_act_conifer1)$coefficients[2,2], summary(tot_act_grass2)$coefficients[2,2]), 
                           p_value=c(summary(tot_act_arable3)$coefficients[2,4], summary(tot_act_broad1)$coefficients[2,4],
                                     summary(tot_act_conifer1)$coefficients[2,4], summary(tot_act_grass2)$coefficients[2,4]),
                           r2=c(r.squaredGLMM(tot_act_arable3)[1,1], r.squaredGLMM(tot_act_broad1)[1,1], 
                                r.squaredGLMM(tot_act_conifer1)[1,1], r.squaredGLMM(tot_act_grass2)[1,1]))           

write.csv(tot_act_land, file="Results/Bats/Barbastella_barbastellus_land_cover.csv", row.names=FALSE)


## Now run models with habitat variables only 
sp_act_hab <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + visit + (1|subcompartment),
                    data = barbastellus, family = binomial, weights=tot_intervals)
summary(sp_act_hab)
AIC(sp_act_hab)

final_mod2 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_arable.1500) + visit + (1|subcompartment),
                    data = barbastellus, family = binomial, weights=tot_intervals)
AIC(final_mod2)

final_mod3 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_broadleaf.500) + visit + (1|subcompartment),
                    data = barbastellus, family = binomial, weights=tot_intervals)
AIC(final_mod3)

final_mod4 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_conifer.500) + visit + (1|subcompartment),
                    data = barbastellus, family = binomial, weights=tot_intervals)
AIC(final_mod4)
summary(final_mod4)
plot(ggpredict(final_mod4, terms = "prop_conifer.500"))

pred <- ggpredict(final_mod4, terms = "prop_conifer.500 [all]")
barbastellus_conifer <- ggplot(pred, aes(x, predicted)) +
        geom_line(linetype="dashed") +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
        geom_point(data=barbastellus, aes(x=prop_conifer.500, y=activity_rate))+
        labs(x="Proportion of conifer within 500m", y="B. barbastellus \nactivity rate")+
        theme_classic()
barbastellus_conifer
ggsave(barbastellus_conifer, file="Graphs/Bats/Barbastellus_conifer_500m.png")

final_mod5 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(complexity_score) + 
                      scale(canopy_openess) + scale(deadwood_snags) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_grassland.1500) + visit + (1|subcompartment),
                    data = barbastellus, family = binomial, weights=tot_intervals)
AIC(final_mod5)

# same model selected (conifer) random vs fixed
# only change was from 3000 to 1500m buffer for arable
# same habitat variables significant (none)

## check model assumptions
testDispersion(sp_act_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = sp_act_hab, plot = F)
plot(simulationOutput) ## no assumptions violated!
## Check for multicolinearity
car::vif(final_mod4) ## looking at GIF^(1/(2*Df)) ==> all look good (under 3)

## put all 5 model summaries into one data frame and save
sp_act_hab_sum <- as.data.frame(coef(summary(sp_act_hab))) #selecting full model coefficient averages
sp_act_hab_sum$parameters <- row.names(sp_act_hab_sum)
row.names(sp_act_hab_sum) <- 1:nrow(sp_act_hab_sum)
sp_act_hab_sum$model <- "habitat"
sp_act_hab_sum$AIC <- AIC(sp_act_hab)

final_mod2_sum <- as.data.frame(coef(summary(final_mod2))) #selecting full model coefficient averages
final_mod2_sum$parameters <- row.names(final_mod2_sum)
row.names(final_mod2_sum) <- 1:nrow(final_mod2_sum)
final_mod2_sum$model <- "arable"
final_mod2_sum$AIC <- AIC(final_mod2)

final_mod3_sum <- as.data.frame(coef(summary(final_mod3))) #selecting full model coefficient averages
final_mod3_sum$parameters <- row.names(final_mod3_sum)
row.names(final_mod3_sum) <- 1:nrow(final_mod3_sum)
final_mod3_sum$model <- "broadleaf"
final_mod3_sum$AIC <- AIC(final_mod3)

final_mod4_sum <- as.data.frame(coef(summary(final_mod4))) #selecting full model coefficient averages
final_mod4_sum$parameters <- row.names(final_mod4_sum)
row.names(final_mod4_sum) <- 1:nrow(final_mod4_sum)
final_mod4_sum$model <- "conifer"
final_mod4_sum$AIC <- AIC(final_mod4)

final_mod5_sum <- as.data.frame(coef(summary(final_mod5))) #selecting full model coefficient averages
final_mod5_sum$parameters <- row.names(final_mod5_sum)
row.names(final_mod5_sum) <- 1:nrow(final_mod5_sum)
final_mod5_sum$model <- "grassland"
final_mod5_sum$AIC <- AIC(final_mod5)

sp_act_hab_mods <- rbind(sp_act_hab_sum, final_mod2_sum, final_mod3_sum, final_mod4_sum, final_mod5_sum)
write.csv(sp_act_hab_mods, file="Results/Bats/Barbastella_barbastellus_habitat.csv", row.names=FALSE)

######################################################################################################################
## Sp8. Eptesicus serotinus
serotinus <- bat_act[bat_act$scientific_name=="Eptesicus serotinus",]
tot_intervals <- serotinus$tot_intervals

round(cor(serotinus[,c("dbh_average","basal_area","per_broadleaf_canopy","complexity_score","canopy_openess",
                          "deadwood_snags","fallen_deadwood", "dis_to_edge", "prop_grassland.3000")]),3)
cor.test(serotinus$basal_area, serotinus$prop_arable.1500) ## significant and r>0.7
## include the variable most strongly correlated with response (sp richness)
cor.test(serotinus$activity_rate, serotinus$per_broadleaf_canopy)
cor.test(serotinus$activity_rate, serotinus$deadwood_snags) 
## remove deadwood snags from analysis

tot_act_arable1 <- glmer(activity_rate ~ prop_arable.500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=serotinus)
summary(tot_act_arable1)
r.squaredGLMM(tot_act_arable1) # 0.0038962118            
tot_act_arable2 <- glmer(activity_rate ~ prop_arable.1500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=serotinus)
summary(tot_act_arable2) 
r.squaredGLMM(tot_act_arable2) # 0.031369852            
tot_act_arable3 <- glmer(activity_rate ~ prop_arable.3000 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=serotinus)
summary(tot_act_arable3)    
r.squaredGLMM(tot_act_arable3) # 0.017357811          
#### ARABLE 1500M #### 

tot_act_broad1 <- glmer(activity_rate ~ prop_broadleaf.500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=serotinus)
summary(tot_act_broad1) 
r.squaredGLMM(tot_act_broad1) # 0.001260761            
tot_act_broad2 <- glmer(activity_rate ~ prop_broadleaf.1500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=serotinus)
summary(tot_act_broad2)
r.squaredGLMM(tot_act_broad2) # 0.022914050             
tot_act_broad3 <- glmer(activity_rate ~ prop_broadleaf.3000 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=serotinus)
summary(tot_act_broad3) 
r.squaredGLMM(tot_act_broad3) # 0.019462188            
#### BROADLEAF 1500M #### 


tot_act_conifer1 <- glmer(activity_rate ~ prop_conifer.500 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=serotinus)
summary(tot_act_conifer1)   
r.squaredGLMM(tot_act_conifer1) # 0.018724410              
tot_act_conifer2 <- glmer(activity_rate ~ prop_conifer.1500 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=serotinus)
summary(tot_act_conifer2)       
r.squaredGLMM(tot_act_conifer2) # 0.024735391             
tot_act_conifer3 <- glmer(activity_rate ~ prop_conifer.3000 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=serotinus)
summary(tot_act_conifer3)         
r.squaredGLMM(tot_act_conifer3) # 0.000025387304             
#### CONIFER 1500M #### (500m previously)

tot_act_grass1 <- glmer(activity_rate ~ prop_grassland.500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=serotinus)
summary(tot_act_grass1)        
r.squaredGLMM(tot_act_grass1) # 0.0059950774               
tot_act_grass2 <- glmer(activity_rate ~ prop_grassland.1500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=serotinus)
summary(tot_act_grass2)                   
r.squaredGLMM(tot_act_grass2) # 0.0009307332               
tot_act_grass3 <- glmer(activity_rate ~ prop_grassland.3000 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=serotinus)
summary(tot_act_grass3)          
r.squaredGLMM(tot_act_grass3) # 0.0078805492              
#### GRASSLAND 3000M #### 

# same buffer sizes selected random vs fixed

# save results
tot_act_land <- data.frame(response=rep("Eptesicus serotinus"),
                           var=c(names(model.frame(tot_act_arable2))[2], names(model.frame(tot_act_broad2))[2], 
                                 names(model.frame(tot_act_conifer2))[2], names(model.frame(tot_act_grass3))[2]),
                           estimate=c(summary(tot_act_arable2)$coefficients[2,1], summary(tot_act_broad2)$coefficients[2,1],
                                      summary(tot_act_conifer2)$coefficients[2,1], summary(tot_act_grass3)$coefficients[2,1]), 
                           se=c(summary(tot_act_arable2)$coefficients[2,2], summary(tot_act_broad2)$coefficients[2,2],
                                summary(tot_act_conifer2)$coefficients[2,2], summary(tot_act_grass3)$coefficients[2,2]), 
                           p_value=c(summary(tot_act_arable2)$coefficients[2,4], summary(tot_act_broad2)$coefficients[2,4],
                                     summary(tot_act_conifer2)$coefficients[2,4], summary(tot_act_grass3)$coefficients[2,4]),
                           r2=c(r.squaredGLMM(tot_act_arable2)[1,1], r.squaredGLMM(tot_act_broad2)[1,1], 
                                r.squaredGLMM(tot_act_conifer2)[1,1], r.squaredGLMM(tot_act_grass3)[1,1]))           

write.csv(tot_act_land, file="Results/Bats/Eptesicus_serotinus_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
sp_act_hab <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + visit + (1|subcompartment),
                    data = serotinus, family = binomial, weights=tot_intervals)
summary(sp_act_hab)
AIC(sp_act_hab)

final_mod2 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_arable.1500) + visit + (1|subcompartment),
                    data = serotinus, family = binomial, weights=tot_intervals)
AIC(final_mod2)

final_mod3 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_broadleaf.1500) + visit + (1|subcompartment),
                    data = serotinus, family = binomial, weights=tot_intervals)
AIC(final_mod3)

final_mod4 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_conifer.1500) + visit + (1|subcompartment),
                    data = serotinus, family = binomial, weights=tot_intervals)
AIC(final_mod4)

final_mod5 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_grassland.3000) + visit + (1|subcompartment),
                    data = serotinus, family = binomial, weights=tot_intervals)
AIC(final_mod5)

# same model selected (habitat only) 
# same habitat variables significant (none)

## check model assumptions
testDispersion(sp_act_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = sp_act_hab, plot = F)
plot(simulationOutput) ## no assumptions violated!
## Check for multicolinearity
car::vif(sp_act_hab) ## looking at GIF^(1/(2*Df)) ==> all look good (under 3)

## put all 5 model summaries into one data frame and save
sp_act_hab_sum <- as.data.frame(coef(summary(sp_act_hab))) #selecting full model coefficient averages
sp_act_hab_sum$parameters <- row.names(sp_act_hab_sum)
row.names(sp_act_hab_sum) <- 1:nrow(sp_act_hab_sum)
sp_act_hab_sum$model <- "habitat"
sp_act_hab_sum$AIC <- AIC(sp_act_hab)

final_mod2_sum <- as.data.frame(coef(summary(final_mod2))) #selecting full model coefficient averages
final_mod2_sum$parameters <- row.names(final_mod2_sum)
row.names(final_mod2_sum) <- 1:nrow(final_mod2_sum)
final_mod2_sum$model <- "arable"
final_mod2_sum$AIC <- AIC(final_mod2)

final_mod3_sum <- as.data.frame(coef(summary(final_mod3))) #selecting full model coefficient averages
final_mod3_sum$parameters <- row.names(final_mod3_sum)
row.names(final_mod3_sum) <- 1:nrow(final_mod3_sum)
final_mod3_sum$model <- "broadleaf"
final_mod3_sum$AIC <- AIC(final_mod3)

final_mod4_sum <- as.data.frame(coef(summary(final_mod4))) #selecting full model coefficient averages
final_mod4_sum$parameters <- row.names(final_mod4_sum)
row.names(final_mod4_sum) <- 1:nrow(final_mod4_sum)
final_mod4_sum$model <- "conifer"
final_mod4_sum$AIC <- AIC(final_mod4)

final_mod5_sum <- as.data.frame(coef(summary(final_mod5))) #selecting full model coefficient averages
final_mod5_sum$parameters <- row.names(final_mod5_sum)
row.names(final_mod5_sum) <- 1:nrow(final_mod5_sum)
final_mod5_sum$model <- "grassland"
final_mod5_sum$AIC <- AIC(final_mod5)

sp_act_hab_mods <- rbind(sp_act_hab_sum, final_mod2_sum, final_mod3_sum, final_mod4_sum, final_mod5_sum)
write.csv(sp_act_hab_mods, file="Results/Bats/Eptesicus_serotinus_habitat.csv", row.names=FALSE)

###########################################################################################################################

## Sp9. Myotis daubentonii
daubentonii <- bat_act[bat_act$scientific_name=="Myotis daubentonii",]
tot_intervals <- daubentonii$tot_intervals

round(cor(daubentonii[,c("dbh_average","basal_area","per_broadleaf_canopy","complexity_score","canopy_openess",
                       "deadwood_snags","fallen_deadwood", "dis_to_edge", "prop_grassland.500")]),3)
cor.test(daubentonii$deadwood_snags, daubentonii$per_broadleaf_canopy) ## significant and r>0.7
## include the variable most strongly correlated with response (sp richness)
cor.test(daubentonii$activity_rate, daubentonii$per_broadleaf_canopy)
cor.test(daubentonii$activity_rate, daubentonii$deadwood_snags) 
## remove deadwood snags from analysis

tot_act_arable1 <- glmer(activity_rate ~ prop_arable.500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=daubentonii)
summary(tot_act_arable1)
r.squaredGLMM(tot_act_arable1) # 0.001527346             
tot_act_arable2 <- glmer(activity_rate ~ prop_arable.1500 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=daubentonii)
summary(tot_act_arable2) 
r.squaredGLMM(tot_act_arable2) # 0.009789410             
tot_act_arable3 <- glmer(activity_rate ~ prop_arable.3000 + visit + (1|subcompartment),
                         family=binomial, weights=tot_intervals, data=daubentonii)
summary(tot_act_arable3)    
r.squaredGLMM(tot_act_arable3) # 0.019141925           
#### ARABLE 3000M #### 

tot_act_broad1 <- glmer(activity_rate ~ prop_broadleaf.500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=daubentonii)
summary(tot_act_broad1) 
r.squaredGLMM(tot_act_broad1) # 0.0026634222             
tot_act_broad2 <- glmer(activity_rate ~ prop_broadleaf.1500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=daubentonii)
summary(tot_act_broad2)
r.squaredGLMM(tot_act_broad2) # 0.00027003019              
tot_act_broad3 <- glmer(activity_rate ~ prop_broadleaf.3000 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=daubentonii)
summary(tot_act_broad3) 
r.squaredGLMM(tot_act_broad3) # 0.00048121014             
#### BROADLEAF 500M #### (1500m previously)


tot_act_conifer1 <- glmer(activity_rate ~ prop_conifer.500 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=daubentonii)
summary(tot_act_conifer1)   
r.squaredGLMM(tot_act_conifer1) # 0.022854302               
tot_act_conifer2 <- glmer(activity_rate ~ prop_conifer.1500 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=daubentonii)
summary(tot_act_conifer2)       
r.squaredGLMM(tot_act_conifer2) # 0.047190872              
tot_act_conifer3 <- glmer(activity_rate ~ prop_conifer.3000 + visit + (1|subcompartment),
                          family=binomial, weights=tot_intervals, data=daubentonii)
summary(tot_act_conifer3)         
r.squaredGLMM(tot_act_conifer3) # 0.018337843              
#### CONIFER 1500M #### 

tot_act_grass1 <- glmer(activity_rate ~ prop_grassland.500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=daubentonii)
summary(tot_act_grass1)        
r.squaredGLMM(tot_act_grass1) # 0.032421360                
tot_act_grass2 <- glmer(activity_rate ~ prop_grassland.1500 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=daubentonii)
summary(tot_act_grass2)                   
r.squaredGLMM(tot_act_grass2) # 0.0063150656                
tot_act_grass3 <- glmer(activity_rate ~ prop_grassland.3000 + visit + (1|subcompartment),
                        family=binomial, weights=tot_intervals, data=daubentonii)
summary(tot_act_grass3)          
r.squaredGLMM(tot_act_grass3) # 0.014684258               
#### GRASSLAND 500M #### (3000m previously)

# same buffer sizes selected random vs fixed

# save results
tot_act_land <- data.frame(response=rep("Myotis daubentonii"),
                           var=c(names(model.frame(tot_act_arable3))[2], names(model.frame(tot_act_broad1))[2], 
                                 names(model.frame(tot_act_conifer2))[2], names(model.frame(tot_act_grass1))[2]),
                           estimate=c(summary(tot_act_arable3)$coefficients[2,1], summary(tot_act_broad1)$coefficients[2,1],
                                      summary(tot_act_conifer2)$coefficients[2,1], summary(tot_act_grass1)$coefficients[2,1]), 
                           se=c(summary(tot_act_arable3)$coefficients[2,2], summary(tot_act_broad1)$coefficients[2,2],
                                summary(tot_act_conifer2)$coefficients[2,2], summary(tot_act_grass1)$coefficients[2,2]), 
                           p_value=c(summary(tot_act_arable3)$coefficients[2,4], summary(tot_act_broad1)$coefficients[2,4],
                                     summary(tot_act_conifer2)$coefficients[2,4], summary(tot_act_grass1)$coefficients[2,4]),
                           r2=c(r.squaredGLMM(tot_act_arable3)[1,1], r.squaredGLMM(tot_act_broad1)[1,1], 
                                r.squaredGLMM(tot_act_conifer2)[1,1], r.squaredGLMM(tot_act_grass1)[1,1]))           

write.csv(tot_act_land, file="Results/Bats/Myotis_daubentonii_land_cover.csv", row.names=FALSE)

## Now run models with habitat variables only 
sp_act_hab <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + visit + (1|subcompartment),
                    data = daubentonii, family = binomial, weights=tot_intervals)
summary(sp_act_hab)
AIC(sp_act_hab)

final_mod2 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_arable.3000) + visit + (1|subcompartment),
                    data = daubentonii, family = binomial, weights=tot_intervals)
AIC(final_mod2)

final_mod3 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_broadleaf.500) + visit + (1|subcompartment),
                    data = daubentonii, family = binomial, weights=tot_intervals)
AIC(final_mod3)

final_mod4 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_conifer.1500) + visit + (1|subcompartment),
                    data = daubentonii, family = binomial, weights=tot_intervals)
AIC(final_mod4)
summary(final_mod4)
plot(ggpredict(final_mod4, terms = "prop_conifer.1500"))

pred <- ggpredict(final_mod4, terms = "prop_conifer.1500 [all]")
daubentonii_conifer <- ggplot(pred, aes(x, predicted)) +
        geom_line() +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
        geom_point(data=daubentonii, aes(x=prop_conifer.1500, y=activity_rate))+
        labs(x="Proportion of conifer within 1500m", y="B. barbastellus \nactivity rate")+
        theme_classic()
daubentonii_conifer
ggsave(daubentonii_conifer, file="Graphs/Bats/Daubentonii_conifer_1500m.png")

final_mod5 <- glmer(activity_rate ~ scale(dbh_average) + scale(basal_area) +
                      scale(per_broadleaf_canopy) + scale(complexity_score) + 
                      scale(canopy_openess) + 
                      scale(fallen_deadwood) + scale(dis_to_edge) + 
                      scale(prop_grassland.500) + visit + (1|subcompartment),
                    data = daubentonii, family = binomial, weights=tot_intervals)
AIC(final_mod5)

# same model selected (conifer)
# same significant habitat variables (none)

## check model assumptions
testDispersion(sp_act_hab) ## no underdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = sp_act_hab, plot = F)
plot(simulationOutput) ## no assumptions violated!
## Check for multicolinearity
car::vif(final_mod4) ## looking at GIF^(1/(2*Df)) ==> all look good (under 3)

## put all 5 model summaries into one data frame and save
sp_act_hab_sum <- as.data.frame(coef(summary(sp_act_hab))) #selecting full model coefficient averages
sp_act_hab_sum$parameters <- row.names(sp_act_hab_sum)
row.names(sp_act_hab_sum) <- 1:nrow(sp_act_hab_sum)
sp_act_hab_sum$model <- "habitat"
sp_act_hab_sum$AIC <- AIC(sp_act_hab)

final_mod2_sum <- as.data.frame(coef(summary(final_mod2))) #selecting full model coefficient averages
final_mod2_sum$parameters <- row.names(final_mod2_sum)
row.names(final_mod2_sum) <- 1:nrow(final_mod2_sum)
final_mod2_sum$model <- "arable"
final_mod2_sum$AIC <- AIC(final_mod2)

final_mod3_sum <- as.data.frame(coef(summary(final_mod3))) #selecting full model coefficient averages
final_mod3_sum$parameters <- row.names(final_mod3_sum)
row.names(final_mod3_sum) <- 1:nrow(final_mod3_sum)
final_mod3_sum$model <- "broadleaf"
final_mod3_sum$AIC <- AIC(final_mod3)

final_mod4_sum <- as.data.frame(coef(summary(final_mod4))) #selecting full model coefficient averages
final_mod4_sum$parameters <- row.names(final_mod4_sum)
row.names(final_mod4_sum) <- 1:nrow(final_mod4_sum)
final_mod4_sum$model <- "conifer"
final_mod4_sum$AIC <- AIC(final_mod4)

final_mod5_sum <- as.data.frame(coef(summary(final_mod5))) #selecting full model coefficient averages
final_mod5_sum$parameters <- row.names(final_mod5_sum)
row.names(final_mod5_sum) <- 1:nrow(final_mod5_sum)
final_mod5_sum$model <- "grassland"
final_mod5_sum$AIC <- AIC(final_mod5)

sp_act_hab_mods <- rbind(sp_act_hab_sum, final_mod2_sum, final_mod3_sum, final_mod4_sum, final_mod5_sum)
write.csv(sp_act_hab_mods, file="Results/Bats/Myotis_daubentonii_habitat.csv", row.names=FALSE)


### Put significant habitat variable graphs together

library(ggpubr)

hab_vars <- ggarrange(brandtii_per_broadleaf, pipistrellus_openness, pipistrellus_complex, pipistrellus_snags, pipistrellus_fallen,
                         pygmaeus_dbh, pygmaeus_open, pygmaeus_snags, pygmaeus_fallen, richness_dbh, labels=c("(a)", 
                         "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)"), nrow=5, ncol=2)
hab_vars
ggsave(hab_vars, file="Graphs/Bats/Habitat_variables_all.png", height=11, width=7)

###################### Graphs

## try and make similar plot to Fuentes 2012 figure 4
richness <- read.csv("Results/Bats/Total_richness_land_cover.csv", header=TRUE)
auritus <- read.csv("Results/Bats/Plecotus_auritus_land_cover.csv", header=TRUE)
pygmaeus <- read.csv("Results/Bats/Pipistrellus_pygmaeus_land_cover.csv", header=TRUE)
pipistrellus <- read.csv("Results/Bats/Pipistrellus_pipistrellus_land_cover.csv", header=TRUE)
noctula <- read.csv("Results/Bats/Nyctalus_noctula_land_cover.csv", header=TRUE)
nattereri <- read.csv("Results/Bats/Myotis_nattereri_land_cover.csv", header=TRUE)
brandtii <- read.csv("Results/Bats/Myotis_mystacinus_brandtii_land_cover.csv", header=TRUE)
daubentonii <- read.csv("Results/Bats/Myotis_daubentonii_land_cover.csv", header=TRUE)
serotinus <- read.csv("Results/Bats/Eptesicus_serotinus_land_cover.csv", header=TRUE)
barbastellus <- read.csv("Results/Bats/Barbastella_barbastellus_land_cover.csv", header=TRUE)

land_cover_r2 <- rbind(richness,auritus,pygmaeus,pipistrellus,noctula,nattereri,brandtii,daubentonii,serotinus,barbastellus)
land_cover_r2$Scale <- sub(".*\\.", "", land_cover_r2$Land_cover)
land_cover_r2$Land_cover <- sub("(.*?)[\\.|:].*", "\\1", land_cover_r2$Land_cover)
land_cover_r2$Scale <- as.factor(land_cover_r2$Scale)
land_cover_r2$Scale <- ordered(land_cover_r2$Scale, levels=c(500,1500,3000))

land_cover <- ggplot(land_cover_r2, aes(x=Scale, y=R2, colour=Land_cover, group=Land_cover))+
        geom_line(lwd=1)+
        geom_point(size=3, shape=15)+
        labs(x="Spatial scale (m)", y=bquote("R"^{"2"}))+
        theme_classic()+
        facet_wrap(~ Response)
land_cover
ggsave(land_cover, file="Graphs/Bats/Landscape_bats.png", height=8, width=14)


## Landscape results
(layout_matrix <- matrix(c(1,1,2,2,3,3,4,4,8,5,5,6,6,7,7,8), nrow = 2, byrow = TRUE))
landscape <- gridExtra::grid.arrange(barbastellus_conifer, daubentonii_conifer, noctula_conifer, noctula_grass, pipistrellus_conifer,
                       pygmaeus_broadleaf, pygmaeus_grass, layout_matrix = layout_matrix)
landscape
ggsave(landscape, file="Graphs/Bats/Land_cover_results.png", height=8, width=14)



