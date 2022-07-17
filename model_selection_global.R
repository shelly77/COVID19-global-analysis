#######################################################
# Supporting R script for: The relative contribution  #
# of environmental, demographic and socioeconomic     #
# factors to variation global of COVID-19 transmission#                                  
#          Compiled by Yihan Cao                      #
#          University of Oslo, Aug, 2021              #
#######################################################

#This R script contains model selection procedure.


################################################################
# complie and load the template
library(TMB)
compile("state_space_global_autocor.cpp")
dyn.load(dynlib("state_space_global_autocor"))

#################################################################
# load data (run data_pre_global.R first)

# save(data,file="data.RData")

# Initialize parameters
parameters <- list(

# fixed effect parameters
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
  
  beta_pweekcases = 0,
  
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
  
  # residuals related parameters
  logalpha0 = 0.2,
  logalpha1 = 0.2,
  logalpha2 = 0.2,
  
 # random effects related parameters
  logsd_spatioLevel = 0.2,
  mu_spatioLevel = rep(0,nlevels(tmb_scale_pweek$spatioLevel)),
  logsd_randm_slop_days = 0.2,
  mu_spatioLevel_days = rep(0,nlevels(tmb_scale_pweek$spatioLevel)),
  logsd_randm_slop_sqdays = 0.2,
  mu_spatioLevel_sqdays = rep(0,nlevels(tmb_scale_pweek$spatioLevel)),
  logsd_randm_slop_mobility = 0.2,
  mu_spatioLevel_mobility= rep(0,nlevels(tmb_scale_pweek$spatioLevel)),
  logsd_randm_slop_mobility_lag1 = 0.2,
  mu_spatioLevel_mobility_lag1 = rep(0,nlevels(tmb_scale_pweek$spatioLevel)),
  logsd_randm_slop_mobility_lag2 = 0.2,
  mu_spatioLevel_mobility_lag2 = rep(0,nlevels(tmb_scale_pweek$spatioLevel)),
  logsd_randm_slop_temp = 0.2,
  mu_spatioLevel_temp = rep(0,nlevels(tmb_scale_pweek$spatioLevel)),
  logsd_randm_slop_sqtemp = 0.2,
  mu_spatioLevel_sqtemp = rep(0,nlevels(tmb_scale_pweek$spatioLevel))
)


# turn off most of parameters in map
# so that we start with simple models
#######################################################
library("numDeriv")
source("functions_maic.R")


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
  
  beta_pweekcases = factor(NA),
  
  beta_population = factor(NA),
  beta_dens = factor(NA),
  beta_mage = factor(NA),
  # 
  beta_gdp = factor(NA),
  beta_newtests = factor(NA),
  # 
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
  
  logsd_spatioLevel = factor(NA),
  mu_spatioLevel = factor(rep(NA,nlevels(tmb_scale_pweek$spatioLevel))),
  
  logsd_randm_slop_days  = factor(NA),
  mu_spatioLevel_days = factor(rep(NA,nlevels(tmb_scale_pweek$spatioLevel))),
  
  logsd_randm_slop_sqdays  = factor(NA),
  mu_spatioLevel_sqdays = factor(rep(NA,nlevels(tmb_scale_pweek$spatioLevel))),
  
  logsd_randm_slop_mobility  = factor(NA),
  mu_spatioLevel_mobility = factor(rep(NA,nlevels(tmb_scale_pweek$spatioLevel))),
  
  logsd_randm_slop_mobility_lag1 = factor(NA),
  mu_spatioLevel_mobility_lag1  = factor(rep(NA,nlevels(tmb_scale_pweek$spatioLevel))),
  
  logsd_randm_slop_mobility_lag2 = factor(NA),
  mu_spatioLevel_mobility_lag2  = factor(rep(NA,nlevels(tmb_scale_pweek$spatioLevel))),
  
  logsd_randm_slop_temp  = factor(NA),
  mu_spatioLevel_temp = factor(rep(NA,nlevels(tmb_scale_pweek$spatioLevel))),
  
  logsd_randm_slop_sqtemp  = factor(NA),
  mu_spatioLevel_sqtemp = factor(rep(NA,nlevels(tmb_scale_pweek$spatioLevel)))
)


#################################################################

model0 <- list(
  parameters = parameters,   
  map=map0,
  args.MakeADFun = list(
    DLL = "state_space_global_autocor",
    data = tmb_scale_pweek,
    random = c("mu_spatioLevel","mu_spatioLevel_days","mu_spatioLevel_sqdays",
               "mu_spatioLevel_mobility","mu_spatioLevel_mobility_lag1","mu_spatioLevel_mobility_lag2",
               "mu_spatioLevel_temp","mu_spatioLevel_sqtemp"),
    silent = FALSE)
)
class(model0) <- "tmbmodel"

#model selection procedure
#################################################################


model1<- update(model0,optimizer="optim",cv=0.01,seed=123)
model1
model1$sdreport
#maic=  4292.753 
#p = 2


model1a <- update(model1,optimizer="nlminb",list(logalpha1 = NULL),
                  cv=0.01,seed=123)
model1a
model1a$sdreport
#maic= 3940.568 
#p = 3


model2 <- update(model1a,optimizer="nlminb",list(beta_temp = NULL,
                                                 beta_windsp = NULL,
                                                 beta_humidity = NULL,
                                                 beta_uv = NULL),
                 cv=0.1,seed=123)
model2
model2$sdreport
#maic= 3909.219 
#p = 7

model2a <- update(model1a,optimizer="nlminb",list(beta_temp_lag1 = NULL,
                                                  beta_windsp_lag1 = NULL,
                                                  beta_humidity_lag1 = NULL,
                                                  beta_uv_lag1 = NULL),
                  cv=0.1,seed=123)
model2a
model2a$sdreport
#maic=  3915.625 
#p = 7

model2b <- update(model1a,optimizer="nlminb",list(beta_temp_lag2 = NULL,
                                                  beta_windsp_lag2 = NULL,
                                                  beta_humidity_lag2 = NULL,
                                                  beta_uv_lag2 = NULL),
                  cv=0.1,seed=123)
model2b
model2b$sdreport
#maic= 3914.904 
#p = 7


model3 <- update(model2,optimizer="nlminb",list(beta_population = NULL),
                 
                 cv=0.02,seed=123)
model3
model3$sdreport
#maic= 3813.76 
#p = 8

model4 <- update(model3,optimizer="nlminb",list(beta_dens = NULL),
                 
                 cv=0.02,seed=123)
model4
model4$sdreport
#maic=3799.339 
#p = 9

model5 <- update(model4,optimizer="nlminb",list(beta_mage = NULL),
                 
                 cv=0.02,seed=123)
model5
model5$sdreport
#maic= 3797.695 
#p = 10

model6 <- update(model5,optimizer="nlminb",list(beta_gdp= NULL),
                 
                 cv=0.02,seed=123)
model6
model6$sdreport
#maic= 3755.954 
#p = 11

model7 <- update(model6,optimizer="nlminb",list(beta_newtests= NULL),
                 
                 cv=0.02,seed=123)
model7
model7$sdreport
#maic=   3617.314 
#p = 12


model8 <- update(model7,optimizer="nlminb",list(beta_mobility = NULL),
                 
                 cv=0.02,seed=123)
model8
model8$sdreport
#maic= 3511.552 
#p = 13

model8a <- update(model7,optimizer="nlminb",list(beta_mobility_lag1 = NULL),
                  
                  cv=0.02,seed=123)
model8a
model8a$sdreport
#maic=  3516.514 
#p = 13

model8b <- update(model7,optimizer="nlminb",list(beta_mobility_lag2 = NULL),
                  
                  cv=0.02,seed=123)
model8b
model8b$sdreport
#maic=   3534.921 
#p = 13

model9 <- update(model8,optimizer="nlminb",list(beta_debtrelief = NULL),
                 
                 cv=0.02,seed=123)
model9
model9$sdreport
#maic=  3511.943 
#p = 14

model9a <- update(model8,optimizer="nlminb",list(beta_debtrelief_lag1 = NULL),
                  
                  cv=0.02,seed=123)
model9a
model9a$sdreport
#maic= 3510.171 
#p = 14

model9b <- update(model8,optimizer="nlminb",list(beta_debtrelief_lag2 = NULL),
                  
                  cv=0.02,seed=123)
model9b
model9b$sdreport
#maic=  3507.503 
#p = 14


model10 <- update(model9b,optimizer="nlminb",list(beta_healthinvestment = NULL),
                  
                  cv=0.02,seed=123)
model10
model10$sdreport
#maic=  3502.161 
#p = 15

model10a <- update(model9b,optimizer="nlminb",list(beta_healthinvestment_lag1 = NULL),
                   
                   cv=0.02,seed=123)
model10a
model10a$sdreport
#maic=  3504.631 
#p = 15

model10b <- update(model9b,optimizer="nlminb",list(beta_healthinvestment_lag2 = NULL),
                   
                   cv=0.02,seed=123)
model10b
model10b$sdreport
#maic=  3508.198 
#p = 15

model11 <- update(model10,optimizer="nlminb",list(beta_pweekcases = NULL),
                  
                  cv=0.02,seed=123)
model11
model11$sdreport
#maic= 2305.321 
#p = 16


model12 <- update(model11,optimizer="nlminb",list(beta_contatracing = NULL),
                  
                  cv=0.02,seed=123)
model12
model12$sdreport
#maic=   2307.271 
#p = 17

model12a <- update(model11,optimizer="nlminb",list(beta_contatracing_lag1 = NULL),
                   
                   cv=0.02,seed=123)
model12a
model12a$sdreport
#maic= 2307.321 
#p = 17

model12b <- update(model11,optimizer="nlminb",list(beta_contatracing_lag2 = NULL),
                   
                   cv=0.02,seed=123)
model12b
model12b$sdreport
#maic= 2307.291 
#p = 17


model13 <- update(model11,optimizer="nlminb",list(beta_days = NULL),
                  
                  cv=0.02,seed=123)
model13
model13$sdreport
#maic= 2202.081 
#p = 17

model14<- update(model13,optimizer="nlminb",list(beta_days_sq = NULL),
                 
                 cv=0.02,seed=123)
model14
model14$sdreport
#maic= 2186.936 
#p = 18


model15 <- update(model14,optimizer="nlminb",list(beta_temp_sq = NULL),
                  
                  cv=0.02,seed=123)
model15
model15$sdreport
#maic= 2181.498 
#p = 19

model16 <- update(model15,optimizer="nlminb",list(logalpha1 = factor(NA)),
                  
                  cv=0.02,seed=123)
model16
model16$sdreport
#maic=  1287 
#p = 18

#######Include random effects
model17 <- update(model16,optimizer="nlminb",list(  mu_spatioLevel = NULL,
                                                    logsd_spatioLevel = NULL),
                  
                  cv=0.1,seed=123)
model17
model17$sdreport
#aic=  1003.696  
#p = 19

#keep either model 18 or model19, not both, otherwise too complicated random effects
model18 <- update(model17,optimizer="nlminb",list(mu_spatioLevel_temp = NULL,
                                                  logsd_randm_slop_temp= NULL,
                                                  mu_spatioLevel_sqtemp = NULL,
                                                  logsd_randm_slop_sqtemp= NULL),
                  cv=0.01,seed=123)
model18
model18$sdreport
#maic=  936.7905 
#p = 21

model19 <- update(model17,optimizer="optim",list(mu_spatioLevel_days = NULL,
                                                 logsd_randm_slop_days = NULL,
                                                 mu_spatioLevel_sqdays = NULL,
                                                 logsd_randm_slop_sqdays = NULL),
                  
                  
                  cv=0.1,seed=123)
model19
model19$sdreport
#aic=   395.2815 
#p = 21


model20 <- update(model19,optimizer="optim",list(beta_mobility = NULL ,
                                                 beta_mobility_lag1 = factor(NA),
                                                 mu_spatioLevel_mobility = NULL,
                                                 logsd_randm_slop_mobility = NULL),
                  
                  
                  cv=0.1,seed=123)
model20
model20$sdreport
#maic= 338.8311 
#p = 22


model20a <- update(model19,optimizer="nlminb",list(beta_mobility = factor(NA),
                                                   beta_mobility_lag1 = NULL,
                                                   mu_spatioLevel_mobility_lag1 = NULL,
                                                   logsd_randm_slop_mobility_lag1 = NULL),
                   
                   cv=0.1,seed=123)
model20a
model20a$sdreport
#maic= 397.0454 
#p = 22

model20b <- update(model19,optimizer="nlminb",list(beta_mobility = factor(NA),
                                                   beta_mobility_lag2 = NULL,
                                                   mu_spatioLevel_mobility_lag2 = NULL,
                                                   logsd_randm_slop_mobility_lag2 = NULL),
                   
                   cv=0.05,seed=123)
model20b
model20b$sdreport
#maic= 461.8485 
#p = 22


###############################################
####So far, model20 is the best
##play around it,check the neighbour models

model21 <- update(model20,optimizer="nlminb",list(beta_temp  =  factor(NA),
                                                  beta_windsp   =  factor(NA),
                                                  beta_humidity   =  factor(NA),
                                                  beta_uv   =  factor(NA),
                                                  beta_temp_lag1 = NULL,
                                                  beta_windsp_lag1 = NULL,
                                                  beta_humidity_lag1 = NULL,
                                                  beta_uv_lag1 = NULL),
                  cv=0.01,seed=123)
model21
model21$sdreport
##maic= 343.8342 
#p = 22

model22 <- update(model20,optimizer="nlminb",list(beta_temp  =  factor(NA),
                                                  beta_windsp   =  factor(NA),
                                                  beta_humidity   =  factor(NA),
                                                  beta_uv   =  factor(NA),
                                                  beta_temp_lag2 = NULL,
                                                  beta_windsp_lag2 = NULL,
                                                  beta_humidity_lag2 = NULL,
                                                  beta_uv_lag2 = NULL),
                  
                  cv=0.1,seed=123)
model22
model22$sdreport
##maic=  342.7529 
#p = 21

model23 <- update(model20,optimizer="nlminb",list(beta_temp_sq = factor(NA)),
                  
                  cv=0.1,seed=123)
model23
model23$sdreport
#m#aic=  337.4667 
#p = 21


model24 <- update(model23,optimizer="nlminb",list(beta_gdp =  factor(NA)),
                  cv=0.01,seed=123)
model24
model24$sdreport
#maic=  336.679 
#p = 20

model25 <- update(model24,optimizer="nlminb",list(beta_dens =  factor(NA)),
                  cv=0.01,seed=123)
model25
model25$sdreport
#aic= 335.0525 
#p = 19

model26 <- update(model25,optimizer="nlminb",list(beta_debtrelief_lag2 =  factor(NA)),
                  cv=0.01,seed=123)
model26
model26$sdreport
#maic=  333.2368 
#p = 18


model27 <- update(model26,optimizer="nlminb",list(beta_healthinvestment  =  factor(NA)),
                  cv=0.01,seed=123)
model27
model27$sdreport
#maic= 331.4557 
#p = 17


model28 <- update(model27,optimizer="nlminb",list(beta_days_sq =  factor(NA)),
                  cv=0.01,seed=123)
model28
model28$sdreport
#maic=  329.4946 
#p = 16

model29 <- update(model28,optimizer="nlminb",list(beta_windsp = factor(NA)),
                  cv=0.01,seed=123)
model29
model29$sdreport
#maic= 327.4956 
#p = 15

model30 <- update(model29,optimizer="nlminb",list(beta_humidity = factor(NA)),
                  cv=0.01,seed=123)
model30
model30$sdreport
#maic=  325.6383 
#p = 14

model31 <- update(model30,optimizer="nlminb",list(beta_mage = factor(NA)),
                  cv=0.01,seed=123)
model31
model31$sdreport
#maic= 326.3435 
#p = 13

model32 <- update(model30,optimizer="nlminb",list(beta_temp = factor(NA)),
                  cv=0.01,seed=123)
model32
model32$sdreport
#maic= 327.9912 
#p = 13

model33 <- update(model30,optimizer="nlminb",list(  beta_mobility_lag1 = NULL,
                                                    mu_spatioLevel_mobility_lag1  = NULL,
                                                    logsd_randm_slop_mobility_lag1 = NULL,
                                                    beta_mobility = factor(NA),
                                                    mu_spatioLevel_mobility= factor(rep(NA,nlevels(tmb_scale_pweek$spatioLevel))),
                                                    logsd_randm_slop_mobility =factor(NA)),
                  cv=0.2,seed=123)
model33
model33$sdreport
#maic=  385.8462 
#p = 14

model34 <- update(model30,optimizer="nlminb",list(  beta_contatracing = NULL),
                  
                  cv=0.2,seed=123)
model34
model34$sdreport
#maic= 321.124 
#p = 15


model35 <- update(model34,optimizer="nlminb",list(beta_debtrelief = NULL),
                  
                  cv=0.2,seed=123)
model35
model35$sdreport
#maic=  320.5725 
#p = 16


#####################################################################

#play around model35 to confirm that it is the best model
#320.5725 

modelN1 <- update(model35,optimizer="nlminb",list( beta_temp= factor(NA),
                                                   beta_uv = factor(NA),
                                                   beta_temp_lag1 = NULL,
                                                   beta_uv_lag1  = NULL),
                  
                  cv=0.1,seed=123)
modelN1
modelN1$sdreport
##aic=  325.6553 
#p = 16
320.5725 - 325.6553 

modelN2 <- update(model35,optimizer="nlminb",list( beta_temp = factor(NA),
                                                   beta_uv = factor(NA),
                                                   beta_temp_lag2 = NULL,
                                                   beta_uv_lag2 = NULL),
                  
                  cv=0.01,seed=123)
modelN2
modelN2$sdreport
#maic= 326.433 
#p = 16
320.5725 - 326.433 

modelN3 <- update(model35,optimizer="nlminb",list(beta_humidity = NULL),
                  
                  cv=0.01,seed=123)
modelN3
modelN3$sdreport
#maic=   322.3718 
#p = 17
320.5725  -  322.3718 

modelN4 <- update(model35,optimizer="nlminb",list(beta_windsp = NULL),
                  
                  cv=0.01,seed=123)
modelN4
modelN4$sdreport
#maic=  322.5568 
#p = 17
320.5725  -322.5568 

modelN5 <- update(model35,optimizer="nlminb",list(beta_population  = factor(NA)),
                  
                  cv=0.01,seed=123)
modelN5
modelN5$sdreport
#maic= 324.8179 
#p = 15
320.5725  -324.8179 

modelN6 <- update(model35,optimizer="nlminb",list(beta_mage = factor(NA)),
                  
                  cv=0.01,seed=123)
modelN6
modelN6$sdreport
#maic=   321.3573 
#p = 15
320.5725  - 321.3573 

modelN7 <- update(model35,optimizer="nlminb",list(beta_gdp = NULL),
                  
                  cv=0.01,seed=123)
modelN7
modelN7$sdreport
#maic= 320.7418 
#p = 17
320.5725  -320.7418 

modelN8 <- update(model35,optimizer="nlminb",list(beta_newtests = factor(NA)),
                  
                  cv=0.01,seed=123)
modelN8
modelN8$sdreport
#maic=  323.3097 
#p = 15
320.5725  - 323.3097 

modelN9 <- update(model35,optimizer="nlminb",list( beta_contatracing = factor(NA)),
                  
                  cv=0.01,seed=123)
modelN9
modelN9$sdreport
#maic=  324.9148 
#p = 15
320.5725  - 324.9148 


modelN10 <- update(model35,optimizer="nlminb",list( beta_contatracing_lag1 = NULL,
                                                    beta_contatracing = factor(NA)),
                   cv=0.01,seed=123)
modelN10
modelN10$sdreport
#maic=   324.202 
#p = 16
320.5725  -  324.202 


modelN11 <- update(model35,optimizer="nlminb",list( beta_contatracing_lag2 = NULL,
                                                    beta_contatracing = factor(NA)),
                   
                   cv=0.01,seed=123)
modelN11
modelN11$sdreport
#maic= 324.9806 
#p = 16
320.5725  - 324.9806 

modelN12 <- update(model35,optimizer="nlminb",list( beta_debtrelief = factor(NA)),
                   
                   cv=0.01,seed=123)
modelN12
modelN12$sdreport
#maic= 321.124 
#p = 15
320.5725  - 321.124 

modelN13 <- update(model35,optimizer="nlminb",list(beta_debtrelief_lag1 = NULL,
                                                   beta_debtrelief = factor(NA)),
                   
                   
                   cv=0.01,seed=123)
modelN13
modelN13$sdreport
#maic= 323.0687 
#p = 16
320.5725  -323.0687 


modelN14 <- update(model35,optimizer="nlminb",list(beta_debtrelief_lag2 = NULL,
                                                   beta_debtrelief = factor(NA)),
                   
                   
                   cv=0.01,seed=123)
modelN14
modelN14$sdreport
#maic=   323.0413 
#p = 16
320.5725  -  323.0413 


modelN15 <- update(model35,optimizer="nlminb",list(beta_healthinvestment = NULL),
                   
                   cv=0.01,seed=123)
modelN15
modelN15$sdreport
#maic=   322.3048 
#p = 16
320.5725  -  322.3048 


modelN16 <- update(model35,optimizer="nlminb",list(beta_healthinvestment_lag1 = NULL),
                   
                   cv=0.01,seed=123)
modelN16
modelN16$sdreport
#maic=   321.3223 
#p = 17
320.5725  - 321.3223 

modelN17 <- update(model35,optimizer="nlminb",list(beta_healthinvestment_lag2 = NULL),
                   
                   cv=0.01,seed=123)
modelN17
modelN17$sdreport
#maic=  322.501 
#p = 17
320.5725  - 322.501 

modelN18 <- update(model35,optimizer="nlminb",list(beta_dens = NULL),
                   
                   cv=0.01,seed=123)
modelN18
modelN18$sdreport
#maic=  322.1646 
#p = 17
320.5725  -322.1646 

modelN19 <- update(model35,optimizer="nlminb",list(beta_temp_sq = NULL),
                   
                   cv=0.01,seed=123)
modelN19
modelN19$sdreport
#maic= 321.9138 
#p = 17
320.5725 - 321.9138 

modelN20 <- update(modelonestep,optimizer="nlminb",list(beta_days_sq = NULL),
                   
                   cv=0.01,seed=123)
modelN20
modelN20$sdreport



#################################################
modelbest <- model35
#maic= 320.5725  
#p = 16
model35$sdreport


#one step to the best model
modelonestep <- update(model0,optimizer="nlminb",list(beta0 = NULL,                   
                                                      beta_temp  = NULL,         
                                                      beta_uv  = NULL,               
                                                      beta_population   = NULL,       
                                                      beta_mage      = NULL,   
                                                      beta_newtests   = NULL,  
                                                      beta_debtrelief = NULL,
                                                      beta_mobility   = NULL,    
                                                      beta_contatracing   = NULL,    
                                                      beta_days  = NULL, 
                                                      beta_pweekcases = NULL,
                                                      
                                                      mu_spatioLevel=  NULL,
                                                      logsd_spatioLevel =  NULL,
                                                      mu_spatioLevel_mobility = NULL,
                                                      logsd_randm_slop_mobility = NULL,
                                                      mu_spatioLevel_days = NULL,
                                                      logsd_randm_slop_days = NULL,
                                                      mu_spatioLevel_sqdays = NULL,
                                                      logsd_randm_slop_sqdays = NULL),
                       
                       cv=0.01,seed=123)
modelonestep
modelonestep$sdreport
#aic = 320.5725  
#p = 16










