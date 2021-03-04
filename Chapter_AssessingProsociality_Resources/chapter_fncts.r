# Function to draw sample with balanced numbers.
sample.balance <- function(pop.data, sam.size){

  		repeat{
      		sample.data <- slice_sample(pop.data, n = sam.size, replace=F) 
      			if(sum(sample.data$Gender=="female")>2 && 
         		sum(sample.data$Gender=="male")>2 && 
         		sum(sample.data$Age.group=="3" & sample.data$Condition=="experimental")>2 &&
         		sum(sample.data$Age.group=="3" & sample.data$Condition=="control")>2 &&
         		sum(sample.data$Age.group=="4" & sample.data$Condition=="experimental")>2 &&
         		sum(sample.data$Age.group=="4" & sample.data$Condition=="control")>2 &&
         		sum(sample.data$Age.group=="5" & sample.data$Condition=="experimental")>2 &&
         		sum(sample.data$Age.group=="5" & sample.data$Condition=="control")>2){break}
    	}  
return(sample.data)
}

# Function to fit models.
fit.model <- function(s.data, version){
  
  # Base models needed for the different versions below.
  sample.model.reduced <- lm(Shared ~ Gender, data = s.data)    
  sample.model.reduced3 <- lm(Shared ~ Condition + Gender, data = s.data)    
  
  if(version==1){
    sample.model <- lm(Shared ~ Condition+as.factor(Age.group) + Gender, data = s.data)  
    sample.model.reduced2 <- lm(Shared ~ Condition + as.factor(Age.group)+ Gender, data = s.data)    
    model.comp = anova(sample.model.reduced, sample.model)
    model.comp2 = anova(sample.model.reduced2, sample.model)
    model.comp3 = anova(sample.model.reduced3, sample.model)
    #
    sample.coeff <- as_tibble(coefficients(summary(sample.model)))
    #return(c(model.comp$"Pr(>F)"[2], model.comp2$"Pr(>F)"[2], model.comp3$"Pr(>F)"[2], mean(sample.coeff[3:4,]$`Estimate`)))
    return(c(model.comp$"Pr(>F)"[2], mean(sample.coeff[3:4,]$`Estimate`)))
  }else if(version==2){
    sample.model <- lm(Shared ~ Condition+Age.group + Gender, data = s.data)  
    sample.model.reduced2 <- lm(Shared ~ Condition + Age.group+ Gender, data = s.data)    
    model.comp = anova(sample.model.reduced, sample.model)
    model.comp2 = anova(sample.model.reduced2, sample.model)
    model.comp3 = anova(sample.model.reduced3, sample.model)
    #
    sample.coeff <- as_tibble(coefficients(summary(sample.model)))
    #return(c(model.comp$"Pr(>F)"[2], model.comp2$"Pr(>F)"[2], model.comp3$"Pr(>F)"[2], sample.coeff[3,]$`Estimate`))
    return(c(model.comp$"Pr(>F)"[2], sample.coeff[3,]$`Estimate`))
  }else if(version==3){
    sample.model <- lm(Shared ~ Condition+Age + Gender, data = s.data)  
    sample.model.reduced2 <- lm(Shared ~ Condition + Age + Gender, data = s.data)    
    model.comp = anova(sample.model.reduced, sample.model)
    model.comp2 = anova(sample.model.reduced2, sample.model)
    model.comp3 = anova(sample.model.reduced3, sample.model)
    #
    sample.coeff <- as_tibble(coefficients(summary(sample.model)))
    #return(c(model.comp$"Pr(>F)"[2], model.comp2$"Pr(>F)"[2], model.comp3$"Pr(>F)"[2], sample.coeff[3,]$`Estimate`))
	return(c(model.comp$"Pr(>F)"[2], sample.coeff[3,]$`Estimate`))
  }else if(version==4){
    contr=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))
    sample.model=glmer(cbind(Shared, (10-Shared)) ~ Condition+Age + Gender + (1|ID), data=s.data, family=binomial, control=contr)
    sample.model.reduced <- glmer(cbind(Shared, (10-Shared)) ~ Gender + (1|ID), data=s.data, family=binomial, control=contr)
    sample.model.reduced2 <- glmer(cbind(Shared, (10-Shared)) ~ Condition + Age + Gender + (1|ID), data=s.data, family=binomial, control=contr)
    sample.model.reduced3 <- glmer(cbind(Shared, (10-Shared)) ~ Condition + Gender + (1|ID), data=s.data, family=binomial, control=contr)
    model.comp = anova(sample.model.reduced, sample.model, test="Chisqu")
    model.comp2 = anova(sample.model.reduced2, sample.model, test="Chisqu")
    model.comp3 = anova(sample.model.reduced3, sample.model, test="Chisqu")
    #
    sample.coeff <- as_tibble(coefficients(summary(sample.model)))
    #return(c(model.comp$"Pr(>Chisq)"[2], model.comp2$"Pr(>Chisq)"[2], model.comp3$"Pr(>Chisq)"[2], sample.coeff[3,]$`Estimate`))
    return(c(model.comp$"Pr(>Chisq)"[2], sample.coeff[3,]$`Estimate`))
  }
}


