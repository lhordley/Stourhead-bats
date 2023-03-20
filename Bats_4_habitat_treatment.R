##########################
#### user: Lisbeth Hordley
#### date: November 2021
#### info: Stourhead bats: calculate differences in habitat between treatments at bat sites

rm(list = ls())
options(scipen=999)

# Packages
library(plyr)
library(dplyr)
library(ggeffects)
library(lme4)
library(glmmTMB)
library(multcomp)
library(lmerTest)
library(data.table)
library(glmmADMB)
library(broom)

## read in data
hab_data <- read.csv("Data/Stourhead_hab_structure_draft.csv", header=TRUE)

## only interested in 7 habitat variables:
# 1. basal area
# 2. average dbh
# 3. canopy openness
# 4. % broadleaf canopy
# 5. complexity score
# 6. deadwood snags
# 7. fallen deadwood

hab_data <- hab_data[,c("plot", "survey_year", "treatment", "subcompartment", "basal_area", "dbh_average",
                        "per_broadleaf_canopy", "fallen_deadwood", "deadwood_snags", "canopy_openess",
                        "complexity_score")]

bat_data <- read.csv("Data/Stourhead_bats_richness_final.csv", header=TRUE)
bat_plots <- unique(bat_data[,c("plot","year")]) ## 32 plots

hab_data <- merge(hab_data, bat_plots, by.x=c("plot", "survey_year"), by.y=c("plot", "year"), all.y=TRUE) ## 155 plots 
hab_data$basal_area <- as.numeric(hab_data$basal_area)
hab_data$canopy_openess <- as.numeric(hab_data$canopy_openess)
summary(hab_data) 

## Change treatment names
# Clear Fell = Preparatory Stage (PS)
# Early Initial transition irregular =  Regeneration Initiation Stage (RIS)
# Irregular = Structural Development Stage (SDS)

hab_data$treatment <- recode(hab_data$treatment, "Clear Fell" = "Stage 1", "Early Transitioning Irregular" = "Stage 2", "Irregular" = "Stage 3")

summary(hab_data)
str(hab_data)

########## 1. BASAL AREA ##########
## poisson model (log linear model had significant homogeneity of variance)
basal_mod <- glmer(basal_area ~ treatment + (1|subcompartment),
                   family="poisson", hab_data)
summary(basal_mod)

sim1 <- DHARMa::simulateResiduals(fittedModel = basal_mod)
DHARMa::testDispersion(sim1)
plot(sim1) ## no assumptions violated

## Then, analyze post hoc comparisons as:
Treat.comp<-glht(basal_mod, mcp(treatment="Tukey"))
summary(Treat.comp)
# save results
res1 <- data.frame(broom::tidy(Treat.comp)[,c(1:2,4:5,7)])
res1$term <- "basal area"

letters <- data.frame(treatment = names(cld(Treat.comp)$mcletters$Letters),
                      L = cld(Treat.comp)$mcletters$Letters)
hdata1 <- merge(hab_data, letters)
maxy_basal <- max(ddply(hdata1, "treatment", summarise, fivenum(basal_area)[5])[,2])*1.1

basal_area_p <- ggplot(hdata1, aes(x=treatment, y = basal_area)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA) + 
  geom_text(aes(x=treatment, y=maxy_basal, label = L), size=7) + 
  ylab("Basal area")+
  xlab("Treatment") +
  theme_classic() +
  theme(text=element_text(size=20))
basal_area_p

########## 2. AVERAGE DBH ##########
## poisson model (log linear model had significant homogeneity of variance)
dbh_mod <- lmer(log(dbh_average) ~ treatment + (1|subcompartment), hab_data)
summary(dbh_mod)

sim1 <- DHARMa::simulateResiduals(fittedModel = dbh_mod)
DHARMa::testDispersion(sim1)
plot(sim1) ## significant test for homogeneity of variance

## Then, analyze post hoc comparisons as:
Treat.comp<-glht(dbh_mod,mcp(treatment='Tukey'))
summary(Treat.comp)
# save results
res2 <- data.frame(broom::tidy(Treat.comp)[,c(1:2,4:5,7)])
res2$term <- "average dbh"
res_final <- rbind(res1, res2)

letters <- data.frame(treatment = names(cld(Treat.comp)$mcletters$Letters),
                      L = cld(Treat.comp)$mcletters$Letters)
hdata1 <- merge(hab_data, letters)
maxy_dbh <- max(ddply(hdata1, "treatment", summarise, fivenum(dbh_average)[5])[,2])*1.1

average_dbh_p <- ggplot(hdata1, aes(x=treatment, y = dbh_average)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA) + 
  geom_text(aes(treatment, y=maxy_dbh, label = L), size=7) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Average DBH")+
  xlab("Treatment") +
  theme_classic() +
  theme(text=element_text(size=20))
average_dbh_p

########## 3. CANOPY OPENNESS ##########
## poisson model (log linear model had significant homogeneity of variance)
canopy_mod <- lmer(log(canopy_openess) ~ treatment + (1|subcompartment), hab_data)
summary(canopy_mod)

sim1 <- DHARMa::simulateResiduals(fittedModel = canopy_mod)
DHARMa::testDispersion(sim1)
plot(sim1) ## all good

## Then, analyze post hoc comparisons as:
Treat.comp<-glht(canopy_mod,mcp(treatment='Tukey'))
summary(Treat.comp)
# save results
res3 <- data.frame(broom::tidy(Treat.comp)[,c(1:2,4:5,7)])
res3$term <- "canopy openness"
res_final <- rbind(res_final, res3)

letters <- data.frame(treatment = names(cld(Treat.comp)$mcletters$Letters),
                      L = cld(Treat.comp)$mcletters$Letters)
hdata1 <- merge(hab_data, letters)
maxy_canopy <- max(ddply(hdata1, "treatment", summarise, fivenum(canopy_openess)[5])[,2])*1.1

canopy_p <- ggplot(hdata1, aes(x=treatment, y = canopy_openess)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA) + 
  geom_text(aes(treatment, y=maxy_canopy, label = L), size=7) + 
  ylab("Canopy openness")+
  xlab("Treatment") +
  theme_classic() +
  theme(text=element_text(size=20))
canopy_p

########## 4. % BROADLEAF CANOPY ##########
## might need zero infalted model
hab_data$treatment = as.factor(hab_data$treatment)
broadleaf_mod <- glmmTMB(per_broadleaf_canopy ~ treatment + (1|subcompartment),
                         data=hab_data, ziformula=~1, family=poisson)
summary(broadleaf_mod)

sim1 <- DHARMa::simulateResiduals(fittedModel = broadleaf_mod)
DHARMa::testDispersion(sim1)
plot(sim1) ## looks good

## Then, analyze post hoc comparisons as:
Treat.comp<-glht(broadleaf_mod,mcp(treatment='Tukey'))
summary(Treat.comp)
# save results
res4<- data.frame(broom::tidy(Treat.comp)[,c(1:2,4:5,7)])
res4$term <- "percentage broadleaf"
res_final <- rbind(res_final, res4)

letters <- data.frame(treatment = names(cld(Treat.comp)$mcletters$Letters),
                      L = cld(Treat.comp)$mcletters$Letters)
hdata1 <- merge(hab_data, letters)
maxy_broad <- max(ddply(hdata1, "treatment", summarise, fivenum(per_broadleaf_canopy)[5])[,2])*1.1

broadleaf_p <- ggplot(hdata1, aes(x=treatment, y = per_broadleaf_canopy)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA) + 
  geom_text(aes(treatment, y=maxy_broad, label = L), size=7) + 
  ylab("Percentage \nbroadleaf canopy")+
  xlab("Treatment") +
  theme_classic()+
  theme(text=element_text(size=20))
broadleaf_p

########## 5. COMPLEXITY SCORE ##########
complex_mod <- lmer(complexity_score ~ treatment + (1|subcompartment), hab_data)
summary(complex_mod) 

sim1 <- DHARMa::simulateResiduals(fittedModel = complex_mod)
DHARMa::testDispersion(sim1)
plot(sim1) ## looks good

## Then, analyze post hoc comparisons as:
Treat.comp<-glht(complex_mod,mcp(treatment='Tukey'))
summary(Treat.comp)
# save results
res5 <- data.frame(broom::tidy(Treat.comp)[,c(1:2,4:5,7)])
res5$term <- "complexity score"
res_final <- rbind(res_final, res5)

letters <- data.frame(treatment = names(cld(Treat.comp)$mcletters$Letters),
                      L = cld(Treat.comp)$mcletters$Letters)
hdata1 <- merge(hab_data, letters)
maxy_complex <- max(ddply(hdata1, "treatment", summarise, fivenum(complexity_score)[5])[,2])*1.1

complex_p <- ggplot(hdata1, aes(x=treatment, y = complexity_score)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA) + 
  geom_text(aes(treatment, y=maxy_complex, label = L), size=7) + 
  ylab("Complexity score")+
  xlab("Treatment") +
  theme_classic() +
  theme(text=element_text(size=20))
complex_p

########## 6. DEADWOOD SNAGS ##########
## poisson model (log linear model had significant homogeneity of variance)
deadwood_mod <- lmer(sqrt(deadwood_snags) ~ treatment + (1|subcompartment), hab_data)
summary(deadwood_mod) ## singular fit error

sim1 <- DHARMa::simulateResiduals(fittedModel = deadwood_mod)
DHARMa::testDispersion(sim1)
plot(sim1) ## deviations and homogeneity of variance
hist(resid(deadwood_mod)) ## good

## Then, analyze post hoc comparisons as:
Treat.comp<-glht(deadwood_mod,mcp(treatment='Tukey'))
summary(Treat.comp)
# save results
res6 <- data.frame(broom::tidy(Treat.comp)[,c(1:2,4:5,7)])
res6$term <- "deadwood snags"
res_final <- rbind(res_final, res6)

letters <- data.frame(treatment = names(cld(Treat.comp)$mcletters$Letters),
                      L = cld(Treat.comp)$mcletters$Letters)
hdata1 <- merge(hab_data, letters)
maxy_dead1 <- max(ddply(hdata1, "treatment", summarise, fivenum(deadwood_snags)[5])[,2])*1.1

deadwood_p1 <- ggplot(hdata1, aes(x=treatment, y = deadwood_snags)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA) + 
  geom_text(aes(treatment, y=maxy_dead1, label = L), size=7) + 
  ylab("Length of deadwood \nsnags (m)")+
  xlab("Treatment") +
  theme_classic() +
  theme(text=element_text(size=20))
deadwood_p1

########## 7. FALLEN DEADWOOD ##########
## poisson model (log linear model had significant homogeneity of variance)
deadwood_mod2 <- lmer(sqrt(fallen_deadwood) ~ treatment + (1|subcompartment), hab_data)
summary(deadwood_mod2) ## singular fit error

sim1 <- DHARMa::simulateResiduals(fittedModel = deadwood_mod2)
DHARMa::testDispersion(sim1)
plot(sim1) ## slight deviation but looks good
hist(resid(deadwood_mod2)) ## good

## Then, analyze post hoc comparisons as:
Treat.comp<-glht(deadwood_mod2,mcp(treatment='Tukey'))
summary(Treat.comp)
# save results
res7 <- data.frame(broom::tidy(Treat.comp)[,c(1:2,4:5,7)])
res7$term <- "fallen deadwood"
res_final <- rbind(res_final, res7)
write.csv(res_final, file="Results/Bats/Habitat_treatment_tukey.csv", row.names=FALSE)

letters <- data.frame(treatment = names(cld(Treat.comp)$mcletters$Letters),
                      L = cld(Treat.comp)$mcletters$Letters)
hdata1 <- merge(hab_data, letters)
maxy_dead2 <- max(ddply(hdata1, "treatment", summarise, fivenum(fallen_deadwood)[5])[,2])*1.1

deadwood_p2 <- ggplot(hdata1, aes(x=treatment, y = fallen_deadwood)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA) + 
  geom_text(aes(treatment, y=maxy_dead2, label = L), size=7) + 
  ylab("Length of fallen \ndeadwood (m)")+
  xlab("Treatment") +
  theme_classic() +
  theme(text=element_text(size=20))
deadwood_p2

###### Put all plots together
# 3 x 3 plots
library(ggpubr)
library(gridExtra)
library(grid)

layout_matrix <- matrix(c(1,1,2,2,3,3,4,4,8,5,5,6,6,7,7,8), nrow = 2, byrow = TRUE)
habitat_plots <- gridExtra::grid.arrange(basal_area_p, average_dbh_p, canopy_p, broadleaf_p,
                                         complex_p, deadwood_p1, deadwood_p2,
                                         layout_matrix = layout_matrix)
habitat_plots

ggsave(habitat_plots, file="Graphs/Habitat/Habitat_treatment_bats.png", height=8, width=14)






