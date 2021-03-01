################################################################################
#
# Power analysis based on simulations.
# Linear Regression
#
# Here a short skit.
#
# Dependent measure: Change in pupil dilation, Gaussian distribution
# Predictor 1 (Condition): Factor, 3 levels
# Predictor 2 (Age): Continuous
# Predictor 3 (Gender): Factor, 2 levels
#
# Last changes 27.02. by RH 
#
################################################################################

# Clear work space and load packages.
rm(list = ls())
set.seed(7)
require(tidyverse)
require(parallel)
#require(foreach)
#require(doParallel)

# Build the population.
population.size = 100000
population.data <- tibble(ID = c(1:population.size) %>% as.factor(),
                          Gender = sample(c("male", "female"), population.size, replace = T) %>%
                            as.factor(),
                          Age.group = sample(c(10,14,24), size=population.size, replace=T),
                          Condition = sample(c("social", "helping", "control"), population.size, replace = T) %>%
                            as.factor()) %>% 
  mutate(Age.days = Age.group+rnorm(population.size, 0, 0.1)) %>%
  mutate(Age.scale = as.vector(scale(Age.days, scale=FALSE)))
  
# Function to fit models.
sample.fit.model <- function(pop.data, version, s.size, betas, pop.size, var.est){

  # Add dependent measure.  
  s.data = pop.data
  matrix.dummy.coding = model.matrix(object =~ Condition + Age.scale + Gender, data = s.data)
  coefs=c("(Intercept)"= 0, "Conditionhelping" = betas, "Conditionsocial" = betas, "Age.scale"= 0 , "Gendermale"= 0)
  s.data <- s.data %>% 
    add_column(Change = as.vector(matrix.dummy.coding%*%coefs)+rnorm(n= pop.size, sd= var.est, mean=0))
  
  # Draw sample.
  sample.data <- slice_sample(s.data, n = s.size, replace=F)
  
  # Base model.
  sample.model.reduced <- lm(Change ~ Gender, data = sample.data)    
  
  if(version==1){
    sample.model <- lm(Change ~ Condition*Age.group + Gender, data = sample.data)  
    model.comp = anova(sample.model.reduced, sample.model)
    return(model.comp$"Pr(>F)"[2])
  }
  
}    

# Main simulation.
to.simulate = 1
nr.sim = 1000
alpha = 0.05
res.var = 0.04 # Based on Hepach & Hermann, 2020. This is the residual variation!
sample.pool = seq(from=20, to = 200, by =10)
beta.effect.pool = seq(from=0, to = 0.2, by = 0.01)
a <- 0
res.simulations = as_tibble(expand_grid(SimulationNr = 1:nr.sim, SampleSize = as.numeric(NA), Beta = as.numeric(NA), Power = as.numeric(NA)))

system.time(  
repeat{
  
  a <- a + 1
  now.beta <- sample(beta.effect.pool, size=1)
  now.sample <- sample(sample.pool, size=1)
  
  f <- function(i) {
    sample.fit.model(population.data, 1, now.sample, now.beta, population.size, res.var)
  }  
  
  # Run simulations.
  # system.time(p.res <- lapply(1:10, f))
  res.simulations[a,2:4] <- t(c(now.sample, now.beta, sum(unlist(mclapply(1:nr.sim, f))<alpha)/nr.sim))
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
  mutate(Cohen.D = Beta.Mean/0.06, Beta.CIhigh = Beta.CIhigh/0.06, Beta.CIlow = Beta.CIlow/0.06)
  
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

save.image(file=paste("Simulation_Results_LR_",Sys.Date(),"_",floor(abs(round(rnorm(1),2))*100),".RData",collapse="",sep = ""))