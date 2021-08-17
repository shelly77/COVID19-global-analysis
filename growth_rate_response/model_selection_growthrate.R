#######################################################
# Supporting R script for: The relative contribution  #
# of environmental, demographic and socioeconomic     #
# factors to global variation of COVID-19 transmission#                                  
#          Compiled by Yihan Cao                      #
#          University of Oslo, Aug, 2021              #
#######################################################

#This R script is used for model selection when the growth
#rate of confirmed cases is taken as the response variable.

################################################################
# complie and load the template

library(TMB)
compile("global_change_growthrate.cpp")
dyn.load(dynlib("global_change_growthrate"))

#################################################################
#load data (run global_change_growthrate.R first)

# Initialize parameters
parameters <- list(
  beta0 = 1,
  beta_temp = 0,
  beta_temp_sq = 0,
  beta_temp_lag1 = 0,
  beta_temp_lag2 = 0,
  
  beta_windsp = 0,
  beta_windsp_lag1 = 0,
  beta_windsp_lag2 = 0,
  
  beta_humidity = 0,
  beta_humidity_lag1 = 0,
  beta_humidity_lag2 = 0,
  
  beta_uv = 0,
  beta_uv_lag1 = 0,
  beta_uv_lag2 = 0,
  

  beta_population = 0,
  beta_dens = 0,
  beta_mage = 0,
  
  beta_gdp = 0,
  beta_newtests = 0,
  
  beta_days = 0,
  beta_days_sq = 0, 
  
  beta_mobility = 0,
  beta_mobility_lag1 = 0,
  beta_mobility_lag2 = 0,
  
  beta_debtrelief = 0,
  beta_debtrelief_lag1 = 0,
  beta_debtrelief_lag2 = 0,
  
  beta_healthinvestment = 0,
  beta_healthinvestment_lag1 = 0,
  beta_healthinvestment_lag2 = 0,
  
  beta_contatracing = 0,
  beta_contatracing_lag1 = 0,
  beta_contatracing_lag2 = 0,
  
  logalpha0 = 0.2,
  logalpha1 = 0.2,
  logalpha2 = 0.2,
  
  logbeta_spatioLevel = 0.2,
  mu_spatioLevel = rep(0,nlevels(data$spatioLevel)),
  logsd_randm_slop_days = 0.2,
  mu_spatioLevel_days = rep(0,nlevels(data$spatioLevel)),
  logsd_randm_slop_sqdays = 0.2,
  mu_spatioLevel_sqdays = rep(0,nlevels(data$spatioLevel)),
  logsd_randm_slop_mobility = 0.2,
  mu_spatioLevel_mobility= rep(0,nlevels(data$spatioLevel)),
  logsd_randm_slop_mobility_lag1 = 0.2,
  mu_spatioLevel_mobility_lag1 = rep(0,nlevels(data$spatioLevel)),
  logsd_randm_slop_mobility_lag2 = 0.2,
  mu_spatioLevel_mobility_lag2 = rep(0,nlevels(data$spatioLevel)),
  logsd_randm_slop_temp = 0.2,
  mu_spatioLevel_temp = rep(0,nlevels(data$spatioLevel)),
  logsd_randm_slop_sqtemp = 0.2,
  mu_spatioLevel_sqtemp = rep(0,nlevels(data$spatioLevel))
)


# turn off most of parameters in map
#######################################################
library("numDeriv")
source("functions.R")

map0 = list(
  beta_temp = factor(NA),
  beta_temp_sq = factor(NA),
  beta_temp_lag1 = factor(NA),
  beta_temp_lag2 = factor(NA),
  
  beta_windsp = factor(NA),
  beta_windsp_lag1 = factor(NA),
  beta_windsp_lag2 = factor(NA),
  
  beta_humidity = factor(NA),
  beta_humidity_lag1 = factor(NA),
  beta_humidity_lag2 = factor(NA),
  
  beta_uv = factor(NA),
  beta_uv_lag1 = factor(NA),
  beta_uv_lag2 = factor(NA),
  
  beta_population = factor(NA),
  beta_dens = factor(NA),
  beta_mage = factor(NA),
  
  beta_gdp = factor(NA),
  beta_newtests = factor(NA),
  
  beta_days = factor(NA),
  beta_days_sq = factor(NA),
  
  beta_mobility = factor(NA),
  beta_mobility_lag1 = factor(NA),
  beta_mobility_lag2 = factor(NA),
  
  beta_debtrelief = factor(NA),
  beta_debtrelief_lag1 = factor(NA),
  beta_debtrelief_lag2 = factor(NA),
  
  beta_healthinvestment = factor(NA),
  beta_healthinvestment_lag1 = factor(NA),
  beta_healthinvestment_lag2 = factor(NA),
  
  beta_contatracing = factor(NA),
  beta_contatracing_lag1 = factor(NA),
  beta_contatracing_lag2 = factor(NA),
  
  logalpha1 = factor(NA),
  logalpha2 = factor(NA),
  
  logbeta_spatioLevel = factor(NA),
  mu_spatioLevel = factor(rep(NA,nlevels(data$spatioLevel))),
  
  logsd_randm_slop_days  = factor(NA),
  mu_spatioLevel_days = factor(rep(NA,nlevels(data$spatioLevel))),
  
  logsd_randm_slop_sqdays  = factor(NA),
  mu_spatioLevel_sqdays = factor(rep(NA,nlevels(data$spatioLevel))),
  
  logsd_randm_slop_mobility  = factor(NA),
  mu_spatioLevel_mobility = factor(rep(NA,nlevels(data$spatioLevel))),
  
  logsd_randm_slop_mobility_lag1 = factor(NA),
  mu_spatioLevel_mobility_lag1  = factor(rep(NA,nlevels(data$spatioLevel))),
  
  logsd_randm_slop_mobility_lag2 = factor(NA),
  mu_spatioLevel_mobility_lag2  = factor(rep(NA,nlevels(data$spatioLevel))),
  
  logsd_randm_slop_temp  = factor(NA),
  mu_spatioLevel_temp = factor(rep(NA,nlevels(data$spatioLevel))),
  
  logsd_randm_slop_sqtemp  = factor(NA),
  mu_spatioLevel_sqtemp = factor(rep(NA,nlevels(data$spatioLevel)))
)


#################################################################
#define a null model
model0 <- list(
  parameters = parameters,   
  map=map0,
  args.MakeADFun = list(
    DLL = "global_change_growthrate",
    data = tmb_data_growrate,
    random = c("mu_spatioLevel","mu_spatioLevel_days","mu_spatioLevel_sqdays",
               "mu_spatioLevel_mobility","mu_spatioLevel_mobility_lag1","mu_spatioLevel_mobility_lag2",
               "mu_spatioLevel_temp","mu_spatioLevel_sqtemp"),
    silent = FALSE)
)
class(model0) <- "tmbmodel"

#model selection procedure
#################################################################
modeltest <- update(model0,optimizer="optim",list(beta0 = NULL,                   
                                                  beta_temp_lag1  = NULL,         
                                                  beta_windsp_lag2  = NULL,         
                                                  beta_humidity_lag2 = NULL,      
                                                  beta_uv   = NULL,               
                                                  beta_population   = NULL,       
                                                  beta_mage      = NULL,        
                                                  beta_newtests   = NULL,        
                                                  beta_mobility_lag2   = NULL,    
                                                  beta_contatracing   = NULL,    
                                                  beta_days  = NULL, 
                                                
                                                  mu_spatioLevel=  NULL,
                                                  logbeta_spatioLevel =  NULL,
                                                  mu_spatioLevel_mobility_lag2 = NULL,
                                                  logsd_randm_slop_mobility_lag2 = NULL,
                                                  mu_spatioLevel_days = NULL,
                                                  logsd_randm_slop_days = NULL,
                                                  mu_spatioLevel_sqdays = NULL,
                                                  logsd_randm_slop_sqdays = NULL),
                  
                  cv=0.1,seed=123)
modeltest
modeltest$sdreport
#aic= 2237.98 
#p = 16

rep_value <- summary(modeltest$sdreport,"report")


#plot the random effect on mobility, but the random effects are insignificant.
#check on the estimates and reported values of the best model
rep_fixed <- summary(modeltest$sdreport,"fixed",p.value = TRUE)
rep_random <- summary(modeltest$sdreport,"random")
rep_value <- summary(modeltest$sdreport,"report")

random_slop_mobility <- as.data.frame(rep_random[rownames(rep_random)=="mu_spatioLevel_mobility_lag2",])
rownames(random_slop_mobility) <- unique(data$spatioLevelfactor)
colnames(random_slop_mobility) <- c("Estimate", "SE")
random_slop_mobility$Country <- rownames(random_slop_mobility)

beta_mobility_fixed <- rep_fixed["beta_mobility_lag2","Estimate"]
random_slops <- random_slop_mobility %>%
                mutate(Estimate = Estimate + beta_mobility_fixed) %>%  
                arrange(Estimate)  %>%  
                mutate(Country = as.character(Country))  

level_order <-random_slops$Country #this vector is used for reorder the countries according to the effect sizes

#ggplot 
plot_random_mobility <- ggplot(random_slops, aes(x=factor(Country,level = level_order), Estimate)) +
  geom_pointrange(aes(ymin = Estimate - 1.96*SE, ymax = Estimate + SE),color="darkred") + 
  geom_hline(yintercept=0,color = "red")  + coord_flip() +
  theme_bw(base_size=10)+
  labs(y = "Effect of population mobility on transmission growth rate change", x = "Country")+
  theme(
    axis.title.x = element_text(color = "darkred", size = 10, face = "bold"),
    axis.title.y = element_text(color = "black", size = 10, face = "bold")
  )

