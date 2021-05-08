################################################################################
#
# Power analysis based on simulations.
# Linear Regression
#
# Here a short skit...
#
# Dependent measure: Change in eye muscle activity, Gaussian distribution
# Predictor 1 (Condition): Factor, 2 levels
# We model the difference between two focal levels ('IJA' and 'NOJA'). The respective 
# betas (reflecting the mean difference in conditions) is modeled to range
# between 0 (no effect) and (0.29 - (-0.32))*2 (twice the largest effect between 
# conditions reported in Addabbo et al., 2019, p. 4).
#
# Predictor 2 (Age): Continuous
# This effect is modeled to be 0.

# Predictor 3 (Gender): Factor, 2 levels
# This effect is modeled to be 0. 
#
# The effect parameters are set below with these commands:
# coefs=c("(Intercept)"= 0, "ConditionIJA" = betas, "Age.scale"= 0 , "Gendermale"= 0)
# beta.effect.pool = seq(from=0, to = (0.29 - (-0.32))*2, by = 0.01)
# 
# We simulated 200 datasets  with same means, sds and ANOVA p-values as Addabbo et al. (2019)
# and calculated the mean residual variation based on these 200 simulated datasets. 
# This value will be used here: res.var = 0.21
# and used here: 
# sample.fit.model(population.data, 1, now.sample, now.beta, population.size, 
# res.var)
#
# We use the same value to estimate standardized Cohen's d:
# mutate(Cohen.D = Beta.Mean/0.06, Beta.CIhigh = Beta.CIhigh/res.var,  # warum /0.06?
# Beta.CIlow = Beta.CIlow/res.var)
#
# Last changes 08.05. by DS 
#
################################################################################

# Clear work space and load packages.
rm(list = ls())
RNGkind("L'Ecuyer-CMRG")
set.seed(7)
require(tidyverse)
require(parallel)

# Detect number of cores:
nr.cores = detectCores()

# Build the population.
population.size = 100000
population.data <- tibble(ID = c(1:population.size) %>% as.factor(),
                          Gender = sample(c("male", "female"), population.size, replace = T) %>%
                            as.factor(),
                          Age.group = sample(c(4,10), size=population.size, replace=T),
                          Condition = sample(c("IJA", "NOJA"), population.size, replace = T) %>%
                            as.factor()) %>% 
  mutate(Age.days = Age.group+rnorm(population.size, 0, 0.1)) %>%
  mutate(Age.scale = as.vector(scale(Age.days, scale=FALSE)))
equal.N.levels <- length(c("IJA", "NOJA", "4", "10"))

# Function to fit models.
sample.fit.model <- function(pop.data, version, s.size, betas, pop.size, var.est){

  # Add dependent measure.  
  s.data = pop.data
  matrix.dummy.coding = model.matrix(object =~ Condition + Age.scale + Gender, data = s.data)
  coefs=c("(Intercept)"= 0, "Conditionhelping" = betas, "Conditionsocial" = betas, "Age.scale"= 0 , "Gendermale"= 0)
  s.data <- s.data %>% 
            add_column(Change = as.vector(matrix.dummy.coding%*%coefs)+rnorm(n= pop.size, sd= var.est, mean=0))
  
  # Draw sample.
  sample.data.4IJA <- slice_sample(s.data[c(which(s.data$Age.group=="4" & s.data$Condition=="IJA")),], n = s.size/equal.N.levels, replace=F) # Unsere Stichprobe wird sein: Identische VPAnzahl pro Altersgruppe, jedes Subjekt hat einen Wert bei IJA und einen Wert bei NOJA, das wird hier nicht abgebildet, oder? 
  sample.data.10IJA <- slice_sample(s.data[c(which(s.data$Age.group=="10" & s.data$Condition=="IJA")),], n = s.size/equal.N.levels, replace=F)
  sample.data.4NOJA <- slice_sample(s.data[c(which(s.data$Age.group=="4" & s.data$Condition=="NOJA")),], n = s.size/equal.N.levels, replace=F)
  sample.data.10NOJA <- slice_sample(s.data[c(which(s.data$Age.group=="10" & s.data$Condition=="NOJA")),], n = s.size/equal.N.levels, replace=F)
  sample.data <- rbind.data.frame(sample.data.4IJA, sample.data.10IJA, sample.data.4NOJA, sample.data.10NOJA)
  
  # Base model.
  sample.model.reduced <- lm(Change ~ Gender, data = sample.data)    
  
  if(version==1){
    sample.model <- lm(Change ~ Condition*Age.group + Gender, data = sample.data)  
    model.comp = anova(sample.model.reduced, sample.model)
    return(model.comp$"Pr(>F)"[2])
  }
  
}    

# Main simulation.
to.simulate = 1 # Set this to actually run simulations. If set to 0, the image will be loaded.
nr.sim = 10 # Number of simulation to be run. Set to 10 'to test the code'. Anything above 10k will take a while... 
            # Note that the same value is used (1) for the number of times a sample is drawn and a model is fitted and (2)
            # for number of 1-beta values that are calculated (see below). 
alpha = 0.05 # The 'significance' level.
res.var = 0.21
sample.pool = seq(from=20, to = 200, by =12)
beta.effect.pool = seq(from=0, to = (0.29 - (-0.32))*2, by = 0.01)
a <- 0
res.simulations = as_tibble(expand_grid(SimulationNr = 1:nr.sim, SampleSize = as.numeric(NA), Beta = as.numeric(NA), Power = as.numeric(NA)))

system.time(  
repeat{
  
  a <- a + 1
  print(a)
  now.beta <- sample(beta.effect.pool, size=1)
  now.sample <- sample(sample.pool, size=1)
  
  f <- function(i) {
    sample.fit.model(population.data, 1, now.sample, now.beta, population.size, res.var)
  }  
  
  # Run simulations.
  # system.time(p.res <- lapply(1:10, f))
  # Running the following with half the number of available cores:
  res.simulations[a,2:4] <- t(c(now.sample, now.beta, sum(unlist(mclapply(1:nr.sim, f, mc.cores = nr.cores/2))<alpha)/nr.sim))
  if(a==nrow(res.simulations)){break}
}
)

res.aggregate <- res.simulations %>%
  filter(Power>=0.8) %>%
  mutate(Beta = abs(Beta)) %>%
  group_by(SampleSize) %>% 
  summarise(Beta.Mean = mean(Beta, na.rm = TRUE),
            Beta.Sd = sd(Beta, na.rm = TRUE),
            Beta.N = n(),
            Beta.Low = min(Beta, na.rm = TRUE),
            Beta.High = max(Beta, na.rm = TRUE),
            Beta.CIlow = quantile(Beta, probs = c(0.025,0.975))[1],
            Beta.CIhigh = quantile(Beta, probs = c(0.025,0.975))[2],
            ) %>% 
  mutate(Cohen.D = Beta.Mean/res.var, Beta.CIhigh = Beta.CIhigh/res.var, Beta.CIlow = Beta.CIlow/res.var)
  
ggplot(data = res.aggregate) + 
  ylim(0, 2) +
  xlim(20,200) +
  theme_bw() + 
  geom_point(mapping = aes(x = SampleSize, y = Beta.CIlow), size=3)+
  ylab("Cohen's D (lower 95% CI; bootstrapped)")+
  ggtitle("The expected (95%) lower effect sizes that achieve power >= 0.8")
  #geom_errorbar(aes(x = SampleSize, ymin=Beta.CIlow, ymax=Beta.CIhigh), colour="black", width=.1)
  

# To get a sense of the distribution:
ggplot(filter(res.simulations, Power>=0.8), aes(x=as.factor(SampleSize), y=Beta)) + geom_boxplot()

save.image(file=paste("Simulation_Results_LR_EMG_",Sys.Date(),"_",floor(abs(round(rnorm(1),2))*100),".RData",collapse="",sep = ""))

# Reference
# Kaiser, J., Crespo-Llado, M. M., Turati, C., & Geangu, E. (2017). 
# The development of spontaneous facial responses to others’ emotions in 
# infancy: An EMG study. Scientific reports, 7(1), 1-10.

# ToDo
# Add option to load Image in case to.simulate is set to 0.
#
