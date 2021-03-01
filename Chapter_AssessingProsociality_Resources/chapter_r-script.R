#
# R-Script to accompany the chapter 'Assessing Prosociality'.
# Last changes 17.12.2020 by RH 
#

# Clear workspace and load packages.
rm(list = ls())
set.seed(7)
require(tidyverse)
require(RcppRoll)
require(lme4)
require(wesanderson)
require(foreach)
require(doParallel)

# Build the population.
population.size = 100000
prob.sharing.cont = 0.5
prob.sharing.exp = 0.6
intercept = log(prob.sharing.cont/(1-prob.sharing.cont))
slope = (log(prob.sharing.exp/(1-prob.sharing.exp)) - log(prob.sharing.cont/(1-prob.sharing.cont)))/(1-0)

population.data <- tibble(ID = c(1:population.size) %>% as.factor(),
                          Gender = sample(c("male", "female"), population.size, replace = T) %>%
                            as.factor(),
                          Age = sample(seq(3,5.9,0.1), size=population.size, replace=T),
                          Condition = sample(c("experimental", "control"), population.size, replace = T) %>%
                            as.factor()) %>% 
  mutate(Age.group = floor(Age)) %>%
  mutate(Age.scale = as.vector(scale(Age, scale=FALSE))) %>%
  mutate(Linear.predictor = intercept + slope*(Condition=="experimental") + 0.22*Age.scale) %>%
  mutate(Shared = rbinom(population.size,10, prob=exp(Linear.predictor)/(1+exp(Linear.predictor))))

# Function to fit models.
fit.model <- function(s.data, version){
  
  # Base models needed for the different versions below.
  sample.model.reduced <- lm(Shared ~ Gender, data = s.data)    
  sample.model.reduced3 <- lm(Shared ~ Condition + Gender, data = s.data)    
  
  if(version==1){
    sample.model <- lm(Shared ~ Condition*as.factor(Age.group) + Gender, data = s.data)  
    sample.model.reduced2 <- lm(Shared ~ Condition + as.factor(Age.group)+ Gender, data = s.data)    
    model.comp = anova(sample.model.reduced, sample.model)
    model.comp2 = anova(sample.model.reduced2, sample.model)
    model.comp3 = anova(sample.model.reduced3, sample.model)
    #
    sample.coeff <- as_tibble(coefficients(summary(sample.model)))
    return(c(model.comp$"Pr(>F)"[2], model.comp2$"Pr(>F)"[2], model.comp3$"Pr(>F)"[2], mean(sample.coeff[3:4,]$`Std. Error`)))
  }else if(version==2){
    sample.model <- lm(Shared ~ Condition*Age.group + Gender, data = s.data)  
    sample.model.reduced2 <- lm(Shared ~ Condition + Age.group+ Gender, data = s.data)    
    model.comp = anova(sample.model.reduced, sample.model)
    model.comp2 = anova(sample.model.reduced2, sample.model)
    model.comp3 = anova(sample.model.reduced3, sample.model)
    #
    sample.coeff <- as_tibble(coefficients(summary(sample.model)))
    return(c(model.comp$"Pr(>F)"[2], model.comp2$"Pr(>F)"[2], model.comp3$"Pr(>F)"[2], sample.coeff[3,]$`Std. Error`))
  }else if(version==3){
    sample.model <- lm(Shared ~ Condition*Age + Gender, data = s.data)  
    sample.model.reduced2 <- lm(Shared ~ Condition + Age + Gender, data = s.data)    
    model.comp = anova(sample.model.reduced, sample.model)
    model.comp2 = anova(sample.model.reduced2, sample.model)
    model.comp3 = anova(sample.model.reduced3, sample.model)
    #
    sample.coeff <- as_tibble(coefficients(summary(sample.model)))
    return(c(model.comp$"Pr(>F)"[2], model.comp2$"Pr(>F)"[2], model.comp3$"Pr(>F)"[2], sample.coeff[3,]$`Std. Error`))
  }else if(version==4){
    contr=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))
    sample.model=glmer(cbind(Shared, (10-Shared)) ~ Condition*Age + Gender + (1|ID), data=s.data, family=binomial, control=contr)
    sample.model.reduced <- glmer(cbind(Shared, (10-Shared)) ~ Gender + (1|ID), data=s.data, family=binomial, control=contr)
    sample.model.reduced2 <- glmer(cbind(Shared, (10-Shared)) ~ Condition + Age + Gender + (1|ID), data=s.data, family=binomial, control=contr)
    sample.model.reduced3 <- glmer(cbind(Shared, (10-Shared)) ~ Condition + Gender + (1|ID), data=s.data, family=binomial, control=contr)
    model.comp = anova(sample.model.reduced, sample.model, test="Chisqu")
    model.comp2 = anova(sample.model.reduced2, sample.model, test="Chisqu")
    model.comp3 = anova(sample.model.reduced3, sample.model, test="Chisqu")
    #
    sample.coeff <- as_tibble(coefficients(summary(sample.model)))
    return(c(model.comp$"Pr(>Chisq)"[2], model.comp2$"Pr(>Chisq)"[2], model.comp3$"Pr(>Chisq)"[2], sample.coeff[3,]$`Std. Error`))
  }
}

# Main simulation.
to.simulate = 0

if(to.simulate == 1){
  numCores <- detectCores()
  registerDoParallel(numCores)
  #
  nr.sim <- 1000
  sample.size <- seq(20,100,by = 5)
  res.simulations = as_tibble(expand_grid(SampleSize = sample.size, SimulationNr = 1:nr.sim)) %>%
  mutate(PVer1 = as.numeric(NA), PVer1.Inter = as.numeric(NA), PVer1.Age = as.numeric(NA), PVer1.SE = as.numeric(NA)) %>%
  mutate(PVer2 = as.numeric(NA), PVer2.Inter = as.numeric(NA), PVer2.Age = as.numeric(NA), PVer2.SE = as.numeric(NA)) %>%
  mutate(PVer3 = as.numeric(NA), PVer3.Inter = as.numeric(NA), PVer3.Age = as.numeric(NA),PVer3.SE = as.numeric(NA)) %>%
  mutate(PVer4 = as.numeric(NA), PVer4.Inter = as.numeric(NA), PVer4.Age = as.numeric(NA), PVer4.SE = as.numeric(NA))

# for(a in 1:length(sample.size)){
foreach (a = 1:length(sample.size)) %dopar% {
    
  #for(b in 1: nr.sim){
  foreach (b = 1:nr.sim) %dopar% {
      
    # Draw sample and ensure that we have equal numbers for Gender, Age (group), & Condition.
    repeat{
      sample.data <- slice_sample(population.data, n = sample.size[a], replace=F)
      
      if(sum(sample.data$Gender=="female")>2 && 
         sum(sample.data$Gender=="male")>2 && 
         sum(sample.data$Age.group=="3" & sample.data$Condition=="experimental")>2 &&
         sum(sample.data$Age.group=="3" & sample.data$Condition=="control")>2 &&
         sum(sample.data$Age.group=="4" & sample.data$Condition=="experimental")>2 &&
         sum(sample.data$Age.group=="4" & sample.data$Condition=="control")>2 &&
         sum(sample.data$Age.group=="5" & sample.data$Condition=="experimental")>2 &&
         sum(sample.data$Age.group=="5" & sample.data$Condition=="control")>2){break}
    }
    
    # Run models versions 1 to 4.
    res.simulations[res.simulations$SampleSize==sample.size[a] & res.simulations$SimulationNr == b,3:6] <- 
      t(fit.model(sample.data, 1))
    #
    res.simulations[res.simulations$SampleSize==sample.size[a] & res.simulations$SimulationNr == b,7:10] <- t(fit.model(sample.data, 2))
    #
    res.simulations[res.simulations$SampleSize==sample.size[a] & res.simulations$SimulationNr == b,11:14] <- t(fit.model(sample.data, 3))
    #
    res.simulations[res.simulations$SampleSize==sample.size[a] & res.simulations$SimulationNr == b,15:18] <- t(suppressMessages(fit.model(sample.data, 4)))
    
  }
}

stopImplicitCluster()
save.image(file=paste("Simulation_Results_",Sys.Date(),"_",floor(abs(round(rnorm(1),2))*100),".RData",collapse="",sep = ""))

}else if(to.simulate == 0){
  load("Simulation_Results_2020-12-08_5.RData")
}


# Plot power of the omnibus test.
res.summarized <- group_by(res.simulations, SampleSize) %>%
  summarise(PowerVer1 = sum(PVer1< 0.05)/nr.sim, PowerVer2 = sum(PVer2< 0.05)/nr.sim, PowerVer3 = sum(PVer3< 0.05)/nr.sim, PowerVer4 = sum(PVer4< 0.05)/nr.sim)
#
res.summarized <- res.summarized %>% 
  mutate(PowerVer1 = roll_mean(PowerVer1, n=5, align="left", , fill=0)) %>%
  mutate(PowerVer2 = roll_mean(PowerVer2, n=5, align="left", , fill=0)) %>%
  mutate(PowerVer3 = roll_mean(PowerVer3, n=5, align="left", , fill=0)) %>%
  mutate(PowerVer4 = roll_mean(PowerVer4, n=5, align="left", , fill=0)) %>%
  pivot_longer(c(`PowerVer1`, `PowerVer2`, `PowerVer3`, `PowerVer4`), names_to = "Version", values_to = "Power")
#
res.summarized$Version <- str_replace_all(res.summarized$Version, "PowerVer1", "ANOVA\n(Age group as factor)\n")
res.summarized$Version <- str_replace_all(res.summarized$Version, "PowerVer2", "Linear Regression\n(Age group as number)\n")
res.summarized$Version <- str_replace_all(res.summarized$Version, "PowerVer3", "Linear Regression\n(Age as number)\n")
res.summarized$Version <- str_replace_all(res.summarized$Version, "PowerVer4", "binomial GLM\n(Age as number)\n")

res.summarized <- res.summarized %>% 
  mutate(across(where(is_character),as_factor))

ggplot(data = res.summarized) + 
  geom_hline(yintercept = 0.8) +
  ylim(0.25, 1) +
  xlim(20,100) +
  theme_bw() + 
  geom_line(mapping = aes(x = SampleSize, y = Power,color = Version), size=1) +
  scale_color_manual(values = wes_palette("Darjeeling1", n=4)) +
  labs(color = "Type of Fiited Model", x = "Sample Size", y = "Statistical Power") +
  geom_text(x=91, y=0.78, label="1 - ", color="grey30") +
  geom_text(x=94, y=0.78, label=expression(beta), color="grey30") +
  geom_text(x=98, y=0.78, label="= .8", color="grey30")


# Plot power of the standard error.
res.summarized <- group_by(res.simulations, SampleSize) %>%
  summarise(PowerVer1.SE = mean(PVer1.SE), PowerVer2.SE = mean(PVer2.SE), PowerVer3.SE = mean(PVer3.SE), PowerVer4.SE = mean(PVer4.SE))

res.summarized <- res.summarized %>% 
  pivot_longer(c(`PowerVer1.SE`, `PowerVer2.SE`, `PowerVer3.SE`, `PowerVer4.SE`), names_to = "Version", values_to = "Power")
#

res.summarized$Version <- str_replace_all(res.summarized$Version, "PowerVer1.SE", "ANOVA\n(Age group as factor)\n")
res.summarized$Version <- str_replace_all(res.summarized$Version, "PowerVer2.SE", "Linear Regression\n(Age group as number)\n")
res.summarized$Version <- str_replace_all(res.summarized$Version, "PowerVer3.SE", "Linear Regression\n(Age as number)\n")
res.summarized$Version <- str_replace_all(res.summarized$Version, "PowerVer4.SE", "binomial GLM\n(Age as number)\n")

res.summarized <- res.summarized %>% 
  mutate(across(where(is_character),as_factor))

ggplot(data = res.summarized) + 
  ylim(0, 1.3) +
  theme_bw() + 
  geom_line(mapping = aes(x = SampleSize, y = Power,color = Version), size=1.75) +
  scale_color_manual(values = wes_palette("Darjeeling1", n=4)) +
  labs(color = "Type of Fiited Model", x = "Sample Size", y = "Standard Error\n (estimating age)")
