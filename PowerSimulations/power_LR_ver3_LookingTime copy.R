################################################################################
#
# Power analysis based on simulations.
# Linear Regression
#
# LOOKING TIME EXAMPLE
#
# This is not working properly yet. the questions is how to simulate the dependent
# measure properly.
#
# Here a short skit...
#
# Dependent measure: Proportion of looking time
#
# Predictor 1 (Condition): Factor, 2 levels
# We base our effect sizes on Experiment 2 in Thiele et al. (submitted). 
# Their study included a social and a control condition. At age 0, the intercept
# was 0.34.
#
# Predictor 2 (Age): Continuous
# This effect is modeled to be beta = 0.02 (SE = 0.0043).

# Predictor 3 (Gender): Factor, 2 levels
# This effect is modeled to be beta = -0.02 (SE = 0.016) . 
#
# We are modeling the residual variation to be 0.041 based on the data provided
# by Thiele et al. (submitted)
#
# The standard deviation to calculate Cohen's D is estimate to be 0.083 based 
# on the data provided by Thiele et al. (submitted)
#
# Last changes 10.04. by RH 
#
################################################################################

# Clear work space and load packages.
rm(list = ls())
RNGkind("L'Ecuyer-CMRG")
set.seed(7)
require(tidyverse)
require(broom)
require(parallel)

# Detect number of cores:
nr.cores = detectCores()

# Build the population.
population.size = 100000
population.data <- tibble(ID = c(1:population.size) %>% as.factor(),
                          Gender = sample(c("male", "female"), population.size, replace = T) %>%
                            as.factor(),
                          Age.group = sample(c(10,14,24), size=population.size, replace=T),
                          Condition = sample(c("social", "helping"), population.size, replace = T) %>%
                            as.factor())

# Function to fit models.
sample.fit.model <- function(pop.data, version, s.size, beta.age, beta.gender, pop.size, var.est){
  # browser()
  # Add dependent measure.  
  s.data = pop.data
  matrix.dummy.coding = model.matrix(object =~ Condition + Age.group + Gender, data = s.data)
  coefs=c("(Intercept)"= 0.34, "Conditionsocial" = 0, "Age.group"= beta.age , "Gendermale"= beta.gender)
  s.data <- s.data %>% 
    add_column(Change = as.vector(matrix.dummy.coding%*%coefs)+rnorm(n= pop.size, sd= var.est, mean=0))
  s.data <- s.data %>% 
    add_column(Change100 = rbinom(pop.size, 100, prob=exp(s.data$Change)/(1+exp(s.data$Change)))) %>%
    mutate(Change100 = Change100/100)
  
  # This distribution looks bad
  #hist(s.data$Change)
  # This one is better:
  # hist(s.data$Change100)
  
  # Draw sample.
  sample.data <- slice_sample(s.data, n = s.size, replace=F)
  
  # Base model.
  sample.model.reduced <- lm(Change100 ~ Gender, data = sample.data)    
  
  if(version==1){
    sample.model <- lm(Change100 ~ Condition*Age.group + Gender, data = sample.data)  
    sample.model = tidy(summary(sample.model)) 
    return(as.numeric(sample.model %>% filter(term=="Age.group") %>% select(p.value)))
    # model.comp = anova(sample.model.reduced, sample.model)
    # return(model.comp$"Pr(>F)"[2])
  }
  
}    

# Main simulation.
to.simulate = 1 # Set this to actually run simulations. If set to 0, the image will be loaded.
nr.sim = 10 # Number of simulation to be run. Set to 10 'to test the code'. Anything above 10k will take a while... 
            # Note that the same value is used (1) for the number of times a sample is drawn and a model is fitted and (2)
            # for number of 1-beta values that are calculated (see below). 
alpha = 0.05 # The 'significance' level.
res.var = 0.041
dm.var = 0.083
sample.pool = seq(from=20, to = 200, by =10)
beta.effect.pool.age = seq(from=0, to = 0.021*2, by = 0.001)
beta.effect.pool.gender = -seq(from=0, to = (0.02)*2, by = 0.01)

a <- 0
res.simulations = as_tibble(expand_grid(SimulationNr = 1:nr.sim, SampleSize = as.numeric(NA), Beta = as.numeric(NA), Power = as.numeric(NA)))

system.time(  
repeat{
  
  a <- a + 1
  now.beta.age <- sample(beta.effect.pool.age, size=1)
  now.beta.gender <- sample(beta.effect.pool.gender, size=1)
  now.sample <- sample(sample.pool, size=1)
  
  run.for.power <- function(i) {
    sample.fit.model(population.data, 1, now.sample, now.beta.age, now.beta.gender, population.size, res.var)
  }  
  
  # Run simulations.
  # system.time(p.res <- lapply(1:10, f))
  # Running the following with half the number of available cores:
  res.simulations[a,2:4] <- t(c(now.sample, now.beta.age, sum(unlist(mclapply(1:nr.sim, run.for.power, mc.cores = nr.cores/2))<alpha)/nr.sim))
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

save.image(file=paste("Simulation_Results_LR_LookingTime_",Sys.Date(),"_",floor(abs(round(rnorm(1),2))*100),".RData",collapse="",sep = ""))

# Reference
# Kaiser, J., Crespo-Llado, M. M., Turati, C., & Geangu, E. (2017). 
# The development of spontaneous facial responses to others’ emotions in 
# infancy: An EMG study. Scientific reports, 7(1), 1-10.
