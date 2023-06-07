##########################
#### user: Lisbeth Hordley
#### date: November 2021
#### info: Stourhead bats: run models testing differences between treatments

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
library(lme4)
library(glmmTMB)
library(AICcmodavg)
library(multcomp)
library(lmerTest)

bat_rich <- read.csv("Data/Stourhead_bats_richness_final.csv", header=TRUE)
bat_act <- read.csv("Data/Stourhead_bats_activity_final.csv", header=TRUE)

##################### DIFFERENCES BETWEEN TREATMENTS
############# GLMS ##############
## Random effects: month + subcomparment
## Predictors: treatment
# Covariates: distance to edge

##### Total species richness
# Which distribution?
par(mfrow=c(2,2))
car::qqp(bat_rich$richness, "norm")
car::qqp(bat_rich$richness, "lnorm")
nbinom <- fitdistr(bat_rich$richness, "Negative Binomial")
car::qqp(bat_rich$richness, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
poisson <- fitdistr(bat_rich$richness, "Poisson")
car::qqp(bat_rich$richness, "pois", lambda=poisson$estimate)

#  Try Poisson
richness_mod <- glmer(richness ~ treatment + scale(dis_to_edge) + visit + (1|subcompartment),
                      data = bat_rich,
                      family = poisson(),
                      na.action = "na.fail")
## singular fit errors from both random effects - not enough variation in subcompartment or visit
summary(richness_mod)

## min temp, cloud cover (both positive), moon cycle, and both distances (negative) significant
## find R2 values
library(MuMIn)
r.squaredGLMM(richness_mod)
## delta marginal = 2% - not good
## check model assumptions
library(DHARMa)
testDispersion(richness_mod) ## underdispersion 
simulationOutput <- simulateResiduals(fittedModel = richness_mod, plot = F)
plot(simulationOutput) ## significant deviation - underdispersion

overdisp_fun <- function(model) {
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  } 
  model.df <- sum(sapply(VarCorr(model), vpars)) + length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model, type = "pearson") # computes pearson residuals
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
} 
overdisp_fun(richness_mod) ## ratio=0.6 - no overdispersion but there is underdispersion
## if ratio >2 = excessive overdispersion

## calc mean and SD of richness in each treatment
bat_rich %>%
  group_by(treatment, visit) %>%  
  summarise(mean_richness = mean(richness),
            sd_richness = sd(richness))

## underdispersion most likely caused by lack of independence
## underdispersed  - mean higher than the variance
## underdispersion leads to conservatism in statistical inference (e.g. decreased power, lowered type I errors)
## type 1 error: significant p value but should be non-sig (i.e. reject the null hypothesis, when the null is correct)

## try using a COM-Poisson 
library(glmmTMB)
library(COMPoissonReg)
library (tibble)
bat_rich <-as_tibble (bat_rich)%>%
  mutate (treatment = factor (treatment)) ## change treatment to factor

richness_mod2 <- glmmTMB(richness ~ treatment + scale(dis_to_edge) + visit + (1|subcompartment),
                         data = bat_rich,
                         family = genpois(link="log"),
                         na.action = "na.fail")

summary(richness_mod2)
testDispersion(richness_mod2) ## looks good 
simulationOutput <- simulateResiduals(fittedModel = richness_mod2, plot = F)
plot(simulationOutput) ## no assumptions violated
qqnorm(resid(richness_mod2))
qqline(resid(richness_mod2))
## looks fine to me 

# Pairwise differences in management effects
glht1 <- glht(richness_mod2, mcp(treatment="Tukey"))
summary(glht1) ## none are significant 

library(broom)
CI <- summary(glht1)
tot_rich_treatment <- data.frame(tidy(CI))
write.csv(tot_rich_treatment, file="Results/Bats/Total_richness_treatment_tukey.csv", row.names=FALSE)

## plot result
library(broom)
CI <- confint(glht1)
tidy(CI)
rich_treatment <- ggplot(CI, aes(lhs, estimate, ymin = lwr, ymax = upr)) +
  geom_hline(yintercept=0, linetype="11", colour="grey60") +
  geom_errorbar(width=0.1) + 
  labs(x="Treatment", y="Difference") +
  geom_point() +
  coord_flip() +
  theme_classic()
rich_treatment
## no significant differences
ggsave(rich_treatment, file="Graphs/Bat_richness_treatment.png")

## alternate plot (boxplot with letters indicating significance)
letters <- data.frame(treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(bat_rich, letters)

rich_treatment2 <- ggplot(mdata1, aes(x=treatment, y = richness)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point() +
  geom_text(aes(treatment, y=max(ddply(bat_rich, "treatment", summarise, fivenum(richness)[5])[,2])*1.1, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("Bat species richness")+
  xlab("Treatment") +
  theme_classic()
rich_treatment2

# 


###################### EACH SPECIES ACTIVITY ~ TREATMENT
## Sp1. Nyctalus noctula

# Binomial model
noctula <- bat_act[bat_act$scientific_name=="Nyctalus noctula",]
tot_intervals <- noctula$tot_intervals
noctula_mod <- glmer(activity_rate ~ treatment + scale(dis_to_edge) + visit + (1|subcompartment),
                     data = noctula,
                     family=binomial, weights=tot_intervals, na.action = "na.fail")
## no errors
summary(noctula_mod)
## check model assumptions
library("DHARMa")
check_gamma_model <- simulateResiduals(fittedModel = noctula_mod, n = 500)
plot(check_gamma_model) ## no assumptions violated

# Pairwise differences in management effects
glht1 <- glht(noctula_mod, mcp(treatment="Tukey"))

summary(glht1) 
confint(glht1)

library(broom)
CI <- summary(glht1)
noctula_rich_treatment <- data.frame(tidy(CI))
write.csv(noctula_rich_treatment, file="Results/Bats/Noctula_richness_treatment_tukey.csv", row.names=FALSE)

## plot result
library(broom)
CI <- confint(glht1)
tidy(CI)
letters <- data.frame(treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(noctula, letters)

noctula_treatment <- ggplot(mdata1, aes(x=treatment, y = activity_rate)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point() +
  geom_text(aes(treatment, y=max(ddply(noctula, "treatment", summarise, fivenum(activity_rate)[5])[,2])*1.1, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("N. noctula activity rate")+
  xlab("Treatment") +
  theme_classic()
noctula_treatment

## Sp2. Plecotus auritus

# Binomial model
auritus <- bat_act[bat_act$scientific_name=="Plecotus auritus",]
tot_intervals <- auritus$tot_intervals
auritus_mod <- glmer(activity_rate ~ treatment + scale(dis_to_edge) + visit + (1|subcompartment),
                     data = auritus, family=binomial, weights=tot_intervals, na.action = "na.fail")
## no errors
summary(auritus_mod)
## check model assumptions
library("DHARMa")
check_gamma_model <- simulateResiduals(fittedModel = auritus_mod, n = 500)
plot(check_gamma_model) ## no assumptions violated

# Pairwise differences in management effects
glht1 <- glht(auritus_mod, mcp(treatment="Tukey"))

summary(glht1) 
confint(glht1)

library(broom)
CI <- summary(glht1)
auritus_rich_treatment <- data.frame(tidy(CI))
write.csv(auritus_rich_treatment, file="Results/Bats/Auritus_richness_treatment_tukey.csv", row.names=FALSE)

## plot result
library(broom)
CI <- confint(glht1)
tidy(CI)
letters <- data.frame(treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(auritus, letters)

auritus_treatment <- ggplot(mdata1, aes(x=treatment, y = activity_rate)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point() +
  geom_text(aes(treatment, y=max(ddply(auritus, "treatment", summarise, fivenum(activity_rate)[5])[,2])*1.1, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("P. auritus activity rate")+
  xlab("Treatment") +
  theme_classic()
auritus_treatment

## Sp3. Myotis nattereri

# Binomial model
nattereri <- bat_act[bat_act$scientific_name=="Myotis nattereri",]
tot_intervals <- nattereri$tot_intervals
nattereri_mod <- glmer(activity_rate ~ treatment + scale(dis_to_edge) + visit + (1|subcompartment),
                       data = nattereri, family=binomial, weights=tot_intervals, 
                       na.action = "na.fail")
## singular fit warning
summary(nattereri_mod)
## check model assumptions
library("DHARMa")
check_gamma_model <- simulateResiduals(fittedModel = nattereri_mod, n = 500)
plot(check_gamma_model) ## no assumptions violated

# Pairwise differences in management effects
glht1 <- glht(nattereri_mod, mcp(treatment="Tukey"))

summary(glht1) 
confint(glht1)

library(broom)
CI <- summary(glht1)
nattereri_rich_treatment <- data.frame(tidy(CI))
write.csv(nattereri_rich_treatment, file="Results/Bats/Nattereri_richness_treatment_tukey.csv", row.names=FALSE)

## plot result
library(broom)
CI <- confint(glht1)
tidy(CI)
letters <- data.frame(treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(nattereri, letters)

nattereri_treatment <- ggplot(mdata1, aes(x=treatment, y = activity_rate)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point() +
  geom_text(aes(treatment, y=max(ddply(nattereri, "treatment", summarise, fivenum(activity_rate)[5])[,2])*1.1, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("M. nattereri activity rate")+
  xlab("Treatment") +
  theme_classic()
nattereri_treatment

## Sp4. Myotis mystacinus/brandtii
# Binomial model
brandtii <- bat_act[bat_act$scientific_name=="Myotis mystacinus/brandtii",]
tot_intervals <- brandtii$tot_intervals
brandtii_mod <- glmer(activity_rate ~ treatment + scale(dis_to_edge) + visit + (1|subcompartment),
                      data = brandtii, family=binomial, weights=tot_intervals, 
                      na.action = "na.fail")
## no errors
summary(brandtii_mod)
## check model assumptions
library("DHARMa")
check_gamma_model <- simulateResiduals(fittedModel = brandtii_mod, n = 500)
plot(check_gamma_model) ## significant quantile deviations - not an assumption of binomial models
plotResiduals(check_gamma_model, brandtii$visit, xlab = "visit", main=NULL)
## all predictors non-significant so continue

# Pairwise differences in management effects
glht1 <- glht(brandtii_mod, mcp(treatment="Tukey"))

summary(glht1) 
confint(glht1)

library(broom)
CI <- summary(glht1)
brandtii_rich_treatment <- data.frame(tidy(CI))
write.csv(brandtii_rich_treatment, file="Results/Bats/Brandtii_richness_treatment_tukey.csv", row.names=FALSE)

## plot result
library(broom)
CI <- confint(glht1)
tidy(CI)
letters <- data.frame(treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(brandtii, letters)

brandtii_treatment <- ggplot(mdata1, aes(x=treatment, y = activity_rate)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point() +
  geom_text(aes(treatment, y=max(ddply(brandtii, "treatment", summarise, fivenum(activity_rate)[5])[,2])*1.1, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("M. brandtii activity rate")+
  xlab("Treatment") +
  theme_classic()
brandtii_treatment

## Sp5. Pipistrellus pipistrellus
# Binomial model
pipistrellus <- bat_act[bat_act$scientific_name=="Pipistrellus pipistrellus",]
tot_intervals <- pipistrellus$tot_intervals
pipistrellus_mod <- glmer(activity_rate ~ treatment + scale(dis_to_edge) + visit + (1|subcompartment),
                          data = pipistrellus, family=binomial, weights=tot_intervals, 
                          na.action = "na.fail")
## no errors
summary(pipistrellus_mod)
## check model assumptions
## check model fit
library("DHARMa")
check_gamma_model <- simulateResiduals(fittedModel = pipistrellus_mod, n = 500)
plot(check_gamma_model) ## no assumptions violated
plotResiduals(check_gamma_model, pipistrellus$visit, xlab = "visit", main=NULL)
## quantile deviations with distance to edge but the test is non-significant
## within group deviations from uniformity significant for visit - much lower in visit 2

# Pairwise differences in management effects
glht1 <- glht(pipistrellus_mod, mcp(treatment="Tukey"))

summary(glht1) 
confint(glht1)

library(broom)
CI <- summary(glht1)
pipistrellus_rich_treatment <- data.frame(tidy(CI))
write.csv(pipistrellus_rich_treatment, file="Results/Bats/Pipistrellus_richness_treatment_tukey.csv", row.names=FALSE)

## plot result
library(broom)
CI <- confint(glht1)
tidy(CI)
letters <- data.frame(treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(pipistrellus, letters)

pipistrellus_treatment <- ggplot(mdata1, aes(x=treatment, y = activity_rate)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point() +
  geom_text(aes(treatment, y=max(ddply(pipistrellus, "treatment", summarise, fivenum(activity_rate)[5])[,2])*1.1, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("P. pipistrellus activity rate")+
  xlab("Treatment") +
  theme_classic()
pipistrellus_treatment

## Sp6. Pipistrellus pygmaeus
# Binomial model
pygmaeus <- bat_act[bat_act$scientific_name=="Pipistrellus pygmaeus",]
tot_intervals <- pygmaeus$tot_intervals
pygmaeus_mod <- glmer(activity_rate ~ treatment + scale(dis_to_edge) + visit + (1|subcompartment),
                      data = pygmaeus, family=binomial, weights=tot_intervals, 
                      na.action = "na.fail")
## no errors
summary(pygmaeus_mod)
## check model assumptions
library("DHARMa")
check_gamma_model <- simulateResiduals(fittedModel = pygmaeus_mod, n = 500)
plot(check_gamma_model) ## significant test for homogeneity of variance
plotResiduals(check_gamma_model, pygmaeus$visit, xlab = "visit", main=NULL)
## significant levene test for homogeneity of variance in treatment
## quantile deviations in distance to edge
## looks much better than common pip so carry on

# Pairwise differences in management effects
glht1 <- glht(pygmaeus_mod, mcp(treatment="Tukey"))

summary(glht1) 
confint(glht1)

library(broom)
CI <- summary(glht1)
pygmaeus_rich_treatment <- data.frame(tidy(CI))
write.csv(pygmaeus_rich_treatment, file="Results/Bats/Pygmaeus_richness_treatment_tukey.csv", row.names=FALSE)

## plot result
library(broom)
CI <- confint(glht1)
tidy(CI)
letters <- data.frame(treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(pygmaeus, letters)

pygmaeus_treatment <- ggplot(mdata1, aes(x=treatment, y = activity_rate)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point() +
  geom_text(aes(treatment, y=max(ddply(pygmaeus, "treatment", summarise, fivenum(activity_rate)[5])[,2])*1.1, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("P. pygmaeus activity rate")+
  xlab("Treatment") +
  theme_classic()
pygmaeus_treatment

## Sp7. Barbastella barbastellus
# Binomial model
barbastellus <- bat_act[bat_act$scientific_name=="Barbastella barbastellus",]
tot_intervals <- barbastellus$tot_intervals
barbastellus_mod <- glmer(activity_rate ~ treatment + scale(dis_to_edge) + visit + (1|subcompartment),
                          data = barbastellus, family=binomial, weights=tot_intervals, 
                          na.action = "na.fail")
## singular fit warning
summary(barbastellus_mod)
## check model assumptions
library("DHARMa")
check_gamma_model <- simulateResiduals(fittedModel = barbastellus_mod, n = 500)
plot(check_gamma_model) ## no assumptions violated

# Pairwise differences in management effects
glht1 <- glht(barbastellus_mod, mcp(treatment="Tukey"))

summary(glht1) 
confint(glht1)

library(broom)
CI <- summary(glht1)
barbastellus_rich_treatment <- data.frame(tidy(CI))
write.csv(barbastellus_rich_treatment, file="Results/Bats/Barbastellus_richness_treatment_tukey.csv", row.names=FALSE)

## plot result
library(broom)
CI <- confint(glht1)
tidy(CI)
letters <- data.frame(treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(barbastellus, letters)

barbastellus_treatment <- ggplot(mdata1, aes(x=treatment, y = activity_rate)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point() +
  geom_text(aes(treatment, y=max(ddply(barbastellus, "treatment", summarise, fivenum(activity_rate)[5])[,2])*1.1, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("B. barbastellus activity rate")+
  xlab("Treatment") +
  theme_classic()
barbastellus_treatment

## Sp8. Eptesicus serotinus
# Binomial model
serotinus <- bat_act[bat_act$scientific_name=="Eptesicus serotinus",]
tot_intervals <- serotinus$tot_intervals
serotinus_mod <- glmer(activity_rate ~ treatment + scale(dis_to_edge) + visit + (1|subcompartment),
                       data = serotinus, family=binomial, weights=tot_intervals, 
                       na.action = "na.fail")
## no errors
summary(serotinus_mod)
## check model assumptions
library("DHARMa")
check_gamma_model <- simulateResiduals(fittedModel = serotinus_mod, n = 500)
plot(check_gamma_model) ## no assumptions violated

# Pairwise differences in management effects
glht1 <- glht(serotinus_mod, mcp(treatment="Tukey"))

summary(glht1) 
confint(glht1)

library(broom)
CI <- summary(glht1)
serotinus_rich_treatment <- data.frame(tidy(CI))
write.csv(serotinus_rich_treatment, file="Results/Bats/Serotinus_richness_treatment_tukey.csv", row.names=FALSE)

## plot result
library(broom)
CI <- confint(glht1)
tidy(CI)
letters <- data.frame(treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(serotinus, letters)

serotinus_treatment <- ggplot(mdata1, aes(x=treatment, y = activity_rate)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point() +
  geom_text(aes(treatment, y=max(ddply(serotinus, "treatment", summarise, fivenum(activity_rate)[5])[,2])*1.1, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("E. serotinus activity rate")+
  xlab("Treatment") +
  theme_classic()
serotinus_treatment


## Sp9. Myotis daubentonii
# Binomial model
daubentonii <- bat_act[bat_act$scientific_name=="Myotis daubentonii",]
tot_intervals <- daubentonii$tot_intervals
daubentonii_mod <- glmer(activity_rate ~ treatment + scale(dis_to_edge) + visit + (1|subcompartment),
                         data = daubentonii, family=binomial, weights=tot_intervals, 
                         na.action = "na.fail")
## no errors
summary(daubentonii_mod)
## check model assumptions
library("DHARMa")
check_gamma_model <- simulateResiduals(fittedModel = daubentonii_mod, n = 500)
plot(check_gamma_model) ## no assumptions violated

# Pairwise differences in management effects
glht1 <- glht(daubentonii_mod, mcp(treatment="Tukey"))

summary(glht1) 
confint(glht1)

library(broom)
CI <- summary(glht1)
daubentonii_rich_treatment <- data.frame(tidy(CI))
write.csv(daubentonii_rich_treatment, file="Results/Bats/Daubentonii_richness_treatment_tukey.csv", row.names=FALSE)

## plot result
library(broom)
CI <- confint(glht1)
tidy(CI)
letters <- data.frame(treatment = names(cld(glht1)$mcletters$Letters),
                      L = cld(glht1)$mcletters$Letters)
mdata1 <- merge(daubentonii, letters)

daubentonii_treatment <- ggplot(mdata1, aes(x=treatment, y = activity_rate)) +
  geom_boxplot(outlier.shape = NA, fill="white", outlier.colour=NA, position=position_dodge(width=0.9)) + 
  geom_point() +
  geom_text(aes(treatment, y=max(ddply(daubentonii, "treatment", summarise, fivenum(activity_rate)[5])[,2])*1.1, label = L), position=position_dodge(width=0.9)) + 
  guides(colour=guide_legend(title="Year")) +
  ylab("M. daubentonii activity rate")+
  xlab("Treatment") +
  theme_classic()
daubentonii_treatment

##### Put all treatment plots together
library(ggpubr)

bat_treatment <- ggarrange(rich_treatment2, auritus_treatment, barbastellus_treatment, brandtii_treatment,
                           daubentonii_treatment, nattereri_treatment, noctula_treatment, pipistrellus_treatment,
                           pygmaeus_treatment, serotinus_treatment,
                           labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)"),
                           hjust=0.05, ncol = 3, nrow = 4)
bat_treatment
ggsave(bat_treatment, file="Graphs/Bats/Bat_treatment_plots.png", height=15, width=12)


########## SRE species richness ########## 
# # Which distribution?
# hist(bat_hab_data$sre_rich)
# par(mfrow=c(2,2))
# car::qqp(bat_hab_data$sre_rich, "norm")
# car::qqp(bat_hab_data$sre_rich, "lnorm")
# nbinom <- fitdistr(bat_hab_data$sre_rich, "Negative Binomial")
# car::qqp(bat_hab_data$sre_rich, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
# poisson <- fitdistr(bat_hab_data$sre_rich, "Poisson")
# car::qqp(bat_hab_data$sre_rich, "pois", lambda=poisson$estimate)
# 
# #  Try Poisson
# richness_mod3 <- glmer(sre_rich ~ treatment + Dist_edge_S + (1|visit) + (1|subcompartment),
#                       data = bat_hab_data,
#                       family = poisson(),
#                       na.action = "na.fail")
# ## singular fit errors from both random effects - not enough variation in subcompartment or visit
# summary(richness_mod3)
# ## min temp, cloud cover (both positive), moon cycle, and both distances (negative) significant
# ## find R2 values
# library(MuMIn)
# r.squaredGLMM(richness_mod3)
# ## delta marginal = 2% - not good
# ## check model assumptions
# library(DHARMa)
# testDispersion(richness_mod3) ## underdispersion 
# simulationOutput <- simulateResiduals(fittedModel = richness_mod3, plot = F)
# plot(simulationOutput) ## significant deviation - underdispersion
# 
# ## try using a COM-Poisson 
# library(glmmTMB)
# library(COMPoissonReg)
# library (tibble)
# bat_hab_data <-as_tibble (bat_hab_data)%>%
#   mutate (treatment = factor (treatment)) ## change treatment to factor
# richness_mod4 <- glmmTMB(mre_rich ~ treatment + Dist_edge_S + (1|visit) + (1|subcompartment),
#                          data = bat_hab_data,
#                          family = "compois",
#                          na.action = "na.fail")
# 
# summary(richness_mod4) ## NaN values - model doesn't work for some reason..
# # testDispersion(richness_mod4) ## looks good 
# # simulationOutput <- simulateResiduals(fittedModel = richness_mod4, plot = F)
# # plot(simulationOutput) ## no assumptions violated
# # qqnorm(resid(richness_mod4))
# # qqline(resid(richness_mod4))
# # ## looks fine to me 
# # 
# # # Pairwise differences in management effects
# # glht1 <- glht(richness_mod4, mcp(treatment="Tukey"))
# # summary(glht1) ## non are significant 
# # 
# # ## plot result
# # library(broom)
# # CI <- confint(glht1)
# # tidy(CI)
# # rich_treatment <- ggplot(CI, aes(lhs, estimate, ymin = lwr, ymax = upr)) +
# #   geom_hline(yintercept=0, linetype="11", colour="grey60") +
# #   geom_errorbar(width=0.1) + 
# #   labs(x="Treatment", y="Difference") +
# #   geom_point() +
# #   coord_flip() +
# #   theme_classic()
# # rich_treatment
# # ## no significant differences


  