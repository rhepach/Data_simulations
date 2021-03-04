################################################################################
#
# R-Script to accompany the chapter 'Assessing Prosociality'.
# 
# Last changes 03.03.2020 by RH 
#
#
################################################################################

# Clear workspace, load packages, and source scripts.
rm(list = ls())
require(tidyverse)
require(parallel)
require(RcppRoll)
require(wesanderson)
require(lme4)
source("chapter_fncts.r")

# Set random number generator.
RNGkind("L'Ecuyer-CMRG")
set.seed(7)

# Detect number of cores.
current.cores <- detectCores()

# Build the population.
# Note that this could also be done later in the processs through drawing from 
# an 'infinite' population.
population.size = 100000
prob.sharing.cont = 0.5
prob.sharing.exp = 0.6
intercept = log(prob.sharing.cont/(1-prob.sharing.cont))
slope = (log(prob.sharing.exp/(1-prob.sharing.exp)) - log(prob.sharing.cont/(1-prob.sharing.cont)))/(1-0)
slope.age = 0.22

population.data <- tibble(ID = c(1:population.size) %>% as.factor(),
                          Gender = sample(c("male", "female"), population.size, replace = T) %>%
                            as.factor(),
                          Age = sample(seq(3,5.9,0.1), size=population.size, replace=T),
                          Condition = sample(c("experimental", "control"), population.size, replace = T) %>%
                            as.factor()) %>% 
  mutate(Age.group = floor(Age)) %>%
  mutate(Age.scale = as.vector(scale(Age, scale=FALSE))) %>%
  mutate(Linear.predictor = intercept + slope*(Condition=="experimental") + slope.age*Age.scale) %>%
  mutate(Shared = rbinom(population.size,10, prob=exp(Linear.predictor)/(1+exp(Linear.predictor))))

# These are the slopes to later calculate the absolute deviation values.
slope.age.fact = mean(coefficients(lm(Shared ~ Condition+as.factor(Age.group) + Gender, data = population.data))[3:4])
slope.age.group = as.vector(coefficients(lm(Shared ~ Condition+Age.group + Gender, data = population.data))[3])
slope.age.cont = as.vector(coefficients(lm(Shared ~ Condition+Age + Gender, data = population.data))[3])

# Main simulation.
to.simulate = 0 # Setting this to 1 will run the simulations. Set to 0 to load last workspace.
nr.sims = 1000
alpha = 0.05
sample.pool = seq(20,120,by = 5)
res.simulations = as_tibble(expand_grid(SampleSize = sample.pool)) %>%
  mutate(PVer1 = as.numeric(NA), PVer2 = as.numeric(NA)) %>%
  mutate(PVer3 = as.numeric(NA), PVer4 = as.numeric(NA)) %>%
  mutate(BetaVer1 = as.numeric(NA), BetaVer2 = as.numeric(NA)) %>%
  mutate(BetaVer3 = as.numeric(NA), BetaVer4 = as.numeric(NA))

if(to.simulate == 1){
system.time(
for(a in 1:length(sample.pool)){
	
	# Run simulations.
	res.sims <- mclapply(1:nr.sims, function(i) {
	
		# Step 1. Draw sample.
		current.sample <- sample.balance(population.data, sample.pool[a])

		# Step 2. Fit model2.
		model1.res <- fit.model(current.sample, 1)
		model2.res <- fit.model(current.sample, 2)
		model3.res <- fit.model(current.sample, 3)
		model4.res <- suppressMessages(fit.model(current.sample, 4))
    return(c(model1.res, model2.res, model3.res,model4.res))
  	}, mc.cores = 6)
	
	# Calculate summaries from the simulations.
	res.matrix = matrix(unlist(res.sims), nrow= 8)
	
	# Subtract simulated estimate from population estimate.
	res.matrix[2,] <- abs(res.matrix[2,]-slope.age.fact)
	res.matrix[4,] <- abs(res.matrix[4,]-slope.age.group)
	res.matrix[6,] <- abs(res.matrix[6,]-slope.age.cont)
	res.matrix[8,] <- abs(res.matrix[8,]-0.22)

	res.power = apply(res.matrix[seq(1,7,b=2),],1,function(x) sum(x<alpha)/nr.sims)
	res.estimate = apply(res.matrix[seq(2,8,b=2),],1,function(x) mean(x))
	res.simulations[a,] <- t(c(sample.pool[a], res.power, res.estimate))
	rm("res.sims", "res.matrix", "res.power", "res.estimate")
})}else if(to.simulate == 0){load("Simulation_Results_2021-03-03_67.RData")}

# Plot Power.
results.power <- res.simulations %>%
	select(SampleSize, PVer1, PVer2, PVer3, PVer4 ) %>%  
	mutate(PVer1 = roll_mean(PVer1, n=5, align="left", , fill=1)) %>%
  	mutate(PVer2 = roll_mean(PVer2, n=5, align="left", , fill=1)) %>%
  	mutate(PVer3 = roll_mean(PVer3, n=5, align="left", , fill=1)) %>%
  	mutate(PVer4 = roll_mean(PVer4, n=5, align="left", , fill=1)) %>%
	pivot_longer(c('PVer1', 'PVer2', 'PVer3', 'PVer4'), names_to = "Version", values_to = "Power")

results.power$Version <- str_replace_all(results.power $Version, "PVer1", "ANOVA\n(Age group as factor)\n")
results.power$Version <- str_replace_all(results.power $Version, "PVer2", "Linear Regression\n(Age group as number)\n")
results.power$Version <- str_replace_all(results.power $Version, "PVer3", "Linear Regression\n(Age as number)\n")
results.power$Version <- str_replace_all(results.power $Version, "PVer4", "binomial GLM\n(Age as number)\n")

results.power <- results.power %>% 
  mutate(across(where(is_character),as_factor))

quartz(width=8,height=4)
ggplot(data = results.power) + 
  geom_hline(yintercept = 0.8) +
  ylim(0.25, 1) +
  xlim(20,100) +
  theme_bw() + 
  geom_line(mapping = aes(x = SampleSize, y = Power,color = Version), size=1) +
  scale_color_manual(values = wes_palette("Darjeeling1", n=4)) +
  labs(color = "Type of Fiited Model", x = "Sample Size", y = "Statistical Power") +
  geom_text(x=91, y=0.76, label="1 - ", color="grey30") +
  geom_text(x=94, y=0.76, label=expression(beta), color="grey30") +
  geom_text(x=98, y=0.76, label="= .8", color="grey30")

# Plot estimation accuracy.
results.esti <- res.simulations %>%
	select(SampleSize, BetaVer1, BetaVer2, BetaVer3, BetaVer4 ) %>%  
	mutate(BetaVer1 = roll_mean(BetaVer1, n=5, align="left", , fill=0)) %>%
  	mutate(BetaVer2 = roll_mean(BetaVer2, n=5, align="left", , fill=0)) %>%
  	mutate(BetaVer3 = roll_mean(BetaVer3, n=5, align="left", , fill=0)) %>%
  	mutate(BetaVer4 = roll_mean(BetaVer4, n=5, align="left", , fill=0)) %>%
	pivot_longer(c('BetaVer1', 'BetaVer2', 'BetaVer3', 'BetaVer4'), names_to = "Version", values_to = "Estimate")

results.esti $Version <- str_replace_all(results.esti $Version, "BetaVer1", "ANOVA\n(Age group as factor)\n")
results.esti $Version <- str_replace_all(results.esti $Version, "BetaVer2", "Linear Regression\n(Age group as number)\n")
results.esti $Version <- str_replace_all(results.esti $Version, "BetaVer3", "Linear Regression\n(Age as number)\n")
results.esti $Version <- str_replace_all(results.esti $Version, "BetaVer4", "binomial GLM\n(Age as number)\n")

results.esti <- results.esti %>% 
  mutate(across(where(is_character),as_factor))

quartz(width=8,height=4)
ggplot(data = results.esti) + 
  ylim(0, 0.55) +
  xlim(20, 100) +
  theme_bw() + 
  geom_line(mapping = aes(x = SampleSize, y = Estimate,color = Version), size=1.5) +
  scale_color_manual(values = wes_palette("Darjeeling1", n=4)) +
  labs(color = "Type of Fitted Model", x = "Sample Size", y = "Deviation \n (estimating 'Age')")

save.image(file=paste("Simulation_Results_",Sys.Date(),"_",floor(abs(round(rnorm(1),2))*100),".RData",collapse="",sep = ""))