#######################################################
# Supporting R script for: The relative contribution  #
# of environmental, demographic and socioeconomic     #
# factors to global variation of COVID-19 transmission#                                  
#          Compiled by Yihan Cao                      #
#          University of Oslo, Jan, 2021              #
#######################################################

#This R script is used for model selection.

#The only difference between this file and "modelSelecton/global2.cpp"
#is that this file contains the number of cases in previous week as a 
#covariate in the model

################################################################
# complie and load the template
#library(TMB)
#compile("NegBinoModelCountry.cpp")
#dyn.load(dynlib("NegBinoModelCountry"))

library(TMB)
#compile("stateSpaceCountry.cpp")
#dyn.load(dynlib("stateSpaceCountry"))
compile("stateSpaceGlobalAutoCor3.cpp")
dyn.load(dynlib("stateSpaceGlobalAutocor3"))

#################################################################
# # load data (run dataPreGlobal.R)

#save(data,file="data.RData")

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
  # 
  beta_uv = 0,
  beta_uv_lag1 = 0,
  beta_uv_lag2 = 0,
  
  beta_pweekcases = 0,
  # 
  beta_population = 0,
  beta_dens = 0,
  beta_mage = 0,
  # 
  beta_gdp = 0,
  beta_newtests = 0,
  # 
  beta_days = 0,
  beta_days_sq = 0,
  # 
  
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
  
  #rhoRan = rep(-0.5,3),
  logsd_spatioLevel = 0.2,
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
  mu_spatioLevel_sqtemp = rep(0,nlevels(tmb_scale_pweek$spatioLevel))
)


# turn off most of parameters in map
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
  # 
  beta_uv = factor(NA),
  beta_uv_lag1 = factor(NA),
  beta_uv_lag2 = factor(NA),
  
  beta_pweekcases = factor(NA),
  
  # 
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
  # 
  logalpha1 = factor(NA),
  logalpha2 = factor(NA),
  
  # rhoRan = factor(rep(NA,3)),
  
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
    DLL = "stateSpaceGlobalAutocor3",
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
#maic= 4186.039 
#p = 2


model1a <- update(model1,optimizer="nlminb",list(logalpha1 = NULL),
                  cv=0.01,seed=123)
model1a
model1a$sdreport
#maic= 3830.073 
#p = 3


#####model1b <- update(model1a,optimizer="nlminb",list(logalpha2 = NULL),
  #                cv=0.01,seed=123)
#####model1b
#######model1b$sdreport
#aic= 3833.986 
#p = 4


model2 <- update(model1a,optimizer="nlminb",list(beta_temp = NULL,
                                                 beta_windsp = NULL,
                                                 beta_humidity = NULL,
                                                 beta_uv = NULL),
                 cv=0.1,seed=123)
model2
model2$sdreport
#maic= 3798.949 
#p = 7

model2a <- update(model1a,optimizer="nlminb",list(beta_temp_lag1 = NULL,
                                                  beta_windsp_lag1 = NULL,
                                                  beta_humidity_lag1 = NULL,
                                                  beta_uv_lag1 = NULL),
                  cv=0.1,seed=123)
model2a
model2a$sdreport
#maic= 3805.419 
#p = 7

model2b <- update(model1a,optimizer="nlminb",list(beta_temp_lag2 = NULL,
                                                  beta_windsp_lag2 = NULL,
                                                  beta_humidity_lag2 = NULL,
                                                  beta_uv_lag2 = NULL),
                  cv=0.1,seed=123)
model2b
model2b$sdreport
#maic= 3807.628 
#p = 7


model3 <- update(model2,optimizer="nlminb",list(beta_population = NULL),
                 
                 cv=0.02,seed=123)
model3
model3$sdreport
#maic= 3704.614 
#p = 8

model4 <- update(model3,optimizer="nlminb",list(beta_dens = NULL),
                 
                 cv=0.02,seed=123)
model4
model4$sdreport
#maic= 3691.407 
#p = 9

model5 <- update(model4,optimizer="nlminb",list(beta_mage = NULL),
                 
                 cv=0.02,seed=123)
model5
model5$sdreport
#maic= 3687.167
#p = 10

model6 <- update(model5,optimizer="nlminb",list(beta_gdp= NULL),
                 
                 cv=0.02,seed=123)
model6
model6$sdreport
#maic= 3641.956
#p = 11

model7 <- update(model6,optimizer="nlminb",list(beta_newtests= NULL),
                 
                 cv=0.02,seed=123)
model7
model7$sdreport
#maic=  3501.307 
#p = 12


model8 <- update(model7,optimizer="nlminb",list(beta_mobility = NULL),
                 
                 cv=0.02,seed=123)
model8
model8$sdreport
#maic= 3422.894 
#p = 13

model8a <- update(model7,optimizer="nlminb",list(beta_mobility_lag1 = NULL),
                  
                  cv=0.02,seed=123)
model8a
model8a$sdreport
#maic=  3424.405 
#p = 13

model8b <- update(model7,optimizer="nlminb",list(beta_mobility_lag2 = NULL),
                  
                  cv=0.02,seed=123)
model8b
model8b$sdreport
#maic=  3436.291 
#p = 13

model9 <- update(model8,optimizer="nlminb",list(beta_debtrelief = NULL),
                 
                 cv=0.02,seed=123)
model9
model9$sdreport
#maic=  3424.85 
#p = 14

model9a <- update(model8,optimizer="nlminb",list(beta_debtrelief_lag1 = NULL),
                  
                  cv=0.02,seed=123)
model9a
model9a$sdreport
#maic=   3423.322 
#p = 14

model9b <- update(model8,optimizer="nlminb",list(beta_debtrelief_lag2 = NULL),
                  
                  cv=0.02,seed=123)
model9b
model9b$sdreport
#maic=   3418.993 
#p = 14


model10 <- update(model9b,optimizer="nlminb",list(beta_healthinvestment = NULL),
                  
                  cv=0.02,seed=123)
model10
model10$sdreport
#maic= 3413.474
#p = 15

model10a <- update(model9b,optimizer="nlminb",list(beta_healthinvestment_lag1 = NULL),
                   
                   cv=0.02,seed=123)
model10a
model10a$sdreport
#maic=  3415.675 
#p = 15

model10b <- update(model9b,optimizer="nlminb",list(beta_healthinvestment_lag2 = NULL),
                   
                   cv=0.02,seed=123)
model10b
model10b$sdreport
#maic=  3419.534 
#p = 15

model11 <- update(model10,optimizer="nlminb",list(beta_pweekcases = NULL),
                   
                   cv=0.02,seed=123)
model11
model11$sdreport
#maic= 2234.986 
#p = 16


model12 <- update(model11,optimizer="nlminb",list(beta_contatracing = NULL),
                  
                  cv=0.02,seed=123)
model12
model12$sdreport
#maic=   2236.959 
#p = 17

model12a <- update(model11,optimizer="nlminb",list(beta_contatracing_lag1 = NULL),
                   
                   cv=0.02,seed=123)
model12a
model12a$sdreport
#maic= 2236.986 
#p = 17

model12b <- update(model11,optimizer="nlminb",list(beta_contatracing_lag2 = NULL),
                   
                   cv=0.02,seed=123)
model12b
model12b$sdreport
#maic= 2236.982 
#p = 17


model13 <- update(model11,optimizer="nlminb",list(beta_days = NULL),
                  
                  cv=0.02,seed=123)
model13
model13$sdreport
#maic= 2137.811 
#p = 17

model14<- update(model13,optimizer="nlminb",list(beta_days_sq = NULL),
                 
                 cv=0.02,seed=123)
model14
model14$sdreport
#maic= 2119.161 
#p = 18


model15 <- update(model14,optimizer="nlminb",list(beta_temp_sq = NULL),
                  
                  cv=0.02,seed=123)
model15
model15$sdreport
#maic=  2112.472
#p = 19

model16 <- update(model15,optimizer="nlminb",list(logalpha1 = factor(NA)),
                  
                  cv=0.02,seed=123)
model16
model16$sdreport
#maic=   1215.427 
#p = 18

#######Include random effects
model17 <- update(model16,optimizer="nlminb",list(  mu_spatioLevel = NULL,
                                                    logsd_spatioLevel = NULL),
                  
                  cv=0.1,seed=123)
model17
model17$sdreport
#aic=   923.1482 
#p = 19

#keep either model 18 or model19, not both, otherwise too complicated random effects
model18 <- update(model17,optimizer="nlminb",list(mu_spatioLevel_temp = NULL,
                                                  logsd_randm_slop_temp= NULL,
                                                  mu_spatioLevel_sqtemp = NULL,
                                                  logsd_randm_slop_sqtemp= NULL),
                  cv=0.01,seed=123)
model18
model18$sdreport
#maic=  856.6193 
#p = 21

model19 <- update(model17,optimizer="optim",list(mu_spatioLevel_days = NULL,
                                                 logsd_randm_slop_days = NULL,
                                                 mu_spatioLevel_sqdays = NULL,
                                                 logsd_randm_slop_sqdays = NULL),
                  
                  
                  cv=0.1,seed=123)
model19
model19$sdreport
#aic=  266.4989 
#p = 21


model20 <- update(model19,optimizer="optim",list(beta_mobility = NULL ,
                                                 beta_mobility_lag1 = factor(NA),
                                                 mu_spatioLevel_mobility = NULL,
                                                 logsd_randm_slop_mobility = NULL),
                  
                  
                  cv=0.1,seed=123)
model20
model20$sdreport
#maic=  230.2352 
#p = 22


model20a <- update(model19,optimizer="nlminb",list(beta_mobility = factor(NA),
                                                  beta_mobility_lag1 = NULL,
                                                  mu_spatioLevel_mobility_lag1 = NULL,
                                                  logsd_randm_slop_mobility_lag1 = NULL),
                  
                  cv=0.1,seed=123)
model20a
model20a$sdreport
#maic= 253.1135 
#p = 22

model20b <- update(model19,optimizer="nlminb",list(beta_mobility = factor(NA),
                                                  beta_mobility_lag2 = NULL,
                                                  mu_spatioLevel_mobility_lag2 = NULL,
                                                  logsd_randm_slop_mobility_lag2 = NULL),
                  
                  cv=0.05,seed=123)
model20b
model20b$sdreport
#maic=  300.2673 
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
##maic= 236.2979 
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
##maic=  235.0132 
#p = 21

model23 <- update(model20,optimizer="nlminb",list(beta_temp_sq = factor(NA)),
              
                   cv=0.1,seed=123)
model23
model23$sdreport
#m#aic=   228.462 
#p = 21


model24 <- update(model23,optimizer="nlminb",list(beta_gdp =  factor(NA)),
                  cv=0.01,seed=123)
model24
model24$sdreport
#maic=  227.5183 
#p = 20

model25 <- update(model24,optimizer="nlminb",list(beta_dens =  factor(NA)),
                  cv=0.01,seed=123)
model25
model25$sdreport
#aic=  225.9339 
#p = 19

model26 <- update(model25,optimizer="nlminb",list(beta_debtrelief_lag2 =  factor(NA)),
                  cv=0.01,seed=123)
model26
model26$sdreport
#maic= 224.2857
#p = 18


model27 <- update(model26,optimizer="nlminb",list(beta_healthinvestment  =  factor(NA)),
                  cv=0.01,seed=123)
model27
model27$sdreport
#maic= 222.4958 
#p = 17


model28 <- update(model27,optimizer="nlminb",list(beta_days_sq =  factor(NA)),
                  cv=0.01,seed=123)
model28
model28$sdreport
#maic= 220.6224 
#p = 16

model29 <- update(model28,optimizer="nlminb",list(beta_windsp = factor(NA)),
                  cv=0.01,seed=123)
model29
model29$sdreport
#maic= 218.7137 
#p = 15

model30 <- update(model29,optimizer="nlminb",list(beta_humidity = factor(NA)),
                  cv=0.01,seed=123)
model30
model30$sdreport
#maic= 216.7766 
#p = 14

model31 <- update(model30,optimizer="nlminb",list(beta_mage = factor(NA)),
                  cv=0.01,seed=123)
model31
model31$sdreport
#maic= 217.3561 
#p = 13

model32 <- update(model30,optimizer="nlminb",list(beta_temp = factor(NA)),
                  cv=0.01,seed=123)
model32
model32$sdreport
#maic=  218.2942 
#p = 13

model33 <- update(model30,optimizer="nlminb",list(  beta_mobility_lag1 = NULL,
                                                    mu_spatioLevel_mobility_lag1  = NULL,
                                                    logsd_randm_slop_mobility_lag1 = NULL,
                                                    beta_mobility = factor(NA),
                                                    mu_spatioLevel_mobility= factor(rep(NA,nlevels(data$spatioLevel))),
                                                    logsd_randm_slop_mobility =factor(NA)),
                  cv=0.2,seed=123)
model33
model33$sdreport
#maic=  242.5004 
#p = 14

model34 <- update(model30,optimizer="nlminb",list(  beta_contatracing = NULL),
                                            
                  cv=0.2,seed=123)
model34
model34$sdreport
#maic=  211.0469 
#p = 15


model35 <- update(model34,optimizer="nlminb",list(beta_debtrelief = NULL),
                  
                  cv=0.2,seed=123)
model35
model35$sdreport
#maic=  208.6521
#p = 16


#####################################################################

#play around model35 to confirm that it is the best model
#216.7766 

modelN1 <- update(model35,optimizer="nlminb",list( beta_temp= factor(NA),
                                                   beta_uv = factor(NA),
                                                   beta_temp_lag1 = NULL,
                                                   beta_uv_lag1  = NULL),
                  
                  cv=0.1,seed=123)
modelN1
modelN1$sdreport
##aic=  215.0133 
#p = 16
208.6521 - 215.0133 

modelN2 <- update(model35,optimizer="nlminb",list( beta_temp = factor(NA),
                                                   beta_uv = factor(NA),
                                                   beta_temp_lag2 = NULL,
                                                   beta_uv_lag2 = NULL),
                  
                  cv=0.01,seed=123)
modelN2
modelN2$sdreport
#maic=214.7715 
#p = 16
208.6521 -214.7715 

modelN3 <- update(model35,optimizer="nlminb",list(beta_humidity = NULL),
                  
                  cv=0.01,seed=123)
modelN3
modelN3$sdreport
#maic=  210.6138 
#p = 17
208.6521 - 210.6138 

modelN4 <- update(model35,optimizer="nlminb",list(beta_windsp = NULL),
                  
                  cv=0.01,seed=123)
modelN4
modelN4$sdreport
#maic=  210.5996
#p = 17
208.6521 -210.5996

modelN5 <- update(model35,optimizer="nlminb",list(beta_population  = factor(NA)),
                  
                  cv=0.01,seed=123)
modelN5
modelN5$sdreport
#maic= 212.8352 
#p = 15
208.6521 -212.8352 

modelN6 <- update(model35,optimizer="nlminb",list(beta_mage = factor(NA)),
                  
                  cv=0.01,seed=123)
modelN6
modelN6$sdreport
#maic=  209.365 
#p = 15
208.6521 -209.365

modelN7 <- update(model35,optimizer="nlminb",list(beta_gdp = NULL),
                  
                  cv=0.01,seed=123)
modelN7
modelN7$sdreport
#maic= 208.8787
#p = 17
208.6521 -208.8787

modelN8 <- update(model35,optimizer="nlminb",list(beta_newtests = factor(NA)),
                  
                  cv=0.01,seed=123)
modelN8
modelN8$sdreport
#maic= 211.5586 
#p = 15
208.6521 -211.5586 

modelN9 <- update(model35,optimizer="nlminb",list( beta_contatracing = factor(NA)),
                  
                  cv=0.01,seed=123)
modelN9
modelN9$sdreport
#maic= 214.0228 
#p = 15
208.6521 -214.0228 


modelN10 <- update(model35,optimizer="nlminb",list( beta_contatracing_lag1 = NULL,
                                                   beta_contatracing = factor(NA)),
                  cv=0.01,seed=123)
modelN10
modelN10$sdreport
#maic=  212.6076 
#p = 16
208.6521 - 212.6076 


modelN11 <- update(model35,optimizer="nlminb",list( beta_contatracing_lag2 = NULL,
                                                    beta_contatracing = factor(NA)),
                                                    
                   cv=0.01,seed=123)
modelN11
modelN11$sdreport
#maic= 213.3388 
#p = 16
#deltaaic = 3.046
208.6521 -213.3388 

modelN12 <- update(model35,optimizer="nlminb",list( beta_debtrelief = factor(NA)),
                  
                  cv=0.01,seed=123)
modelN12
modelN12$sdreport
#maic= 211.0469 
#p = 15
208.6521 -211.0469 

modelN13 <- update(model35,optimizer="nlminb",list(beta_debtrelief_lag1 = NULL,
                                                   beta_debtrelief = factor(NA)),
                                                
                   
                   cv=0.01,seed=123)
modelN13
modelN13$sdreport
#maic=   212.8971 
#p = 16
#deltaaic = 1.11
208.6521 -212.8971 


modelN14 <- update(model35,optimizer="nlminb",list(beta_debtrelief_lag2 = NULL,
                                                   beta_debtrelief = factor(NA)),
                                                
                   
                   cv=0.01,seed=123)
modelN14
modelN14$sdreport
#maic=   212.8444 
#p = 16
#deltaaic = 1.11
208.6521 - 212.8444 


modelN15 <- update(model35,optimizer="nlminb",list(beta_healthinvestment = NULL),
                   
                   cv=0.01,seed=123)
modelN15
modelN15$sdreport
#maic=   210.3786 
#p = 16
#deltaaic = 0.025
208.6521 - 210.3786 


modelN16 <- update(model35,optimizer="nlminb",list(beta_healthinvestment_lag1 = NULL),
                   
                   cv=0.01,seed=123)
modelN16
modelN16$sdreport
#aic=   209.4288 
#p = 17
208.6521 - 209.4288 

modelN17 <- update(model35,optimizer="nlminb",list(beta_healthinvestment_lag2 = NULL),
                   
                   cv=0.01,seed=123)
modelN17
modelN17$sdreport
#aic= 210.6049 
#p = 17
208.6521 -210.6049 

modelN18 <- update(model35,optimizer="nlminb",list(beta_dens = NULL),
                   
                   cv=0.01,seed=123)
modelN18
modelN18$sdreport
#aic= 210.2764 
#p = 17
208.6521 -210.2764 



#################################################
modelbest <- model35
#aic= 208.6521 
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
#aic = 208.6521 
#p = 16










