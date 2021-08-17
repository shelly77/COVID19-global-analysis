
### in this script, we use the one data source that integrates severall
# data sources. see https://rs-delve.github.io/data_software/global-dataset.html

#load the data file
globalraw  <- read.csv("raw_data/combined_dataset_latest.csv", header=T)
head(globalraw)
#globalraw[globalraw$country_name=="Russian Federation",]
###############################################
#A customized function to calculate mean-normalized data
normalizefun <- function(x) { 
  x <- as.matrix(x)
  minAttr=apply(x, 2, min)
  meanAttr=apply(x, 2, mean)
  maxAttr=apply(x, 2, max)
  x <- sweep(x, 2, meanAttr, FUN="-") 
  x <- sweep(x, 2,  maxAttr-minAttr, "/") 
  attr(x, 'normalized:min') = minAttr
  attr(x, 'normalized:max') = maxAttr
  return (x)
} 


#apply the function to our mobility variables
mymatrix <- as.matrix(na.omit(globalraw[,47:52]),ncol=6)
#linear transformation such as normalization and standardization
#doesn't change pearson correlation
normalizedData <- normalizefun(mymatrix)

#calculate the correlation matrix of the mobility variables 
#and check the P values
library(Hmisc)
res <-rcorr(normalizedData,type="pearson")
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
flattenCorrMatrix(res$r, res$P)
# the mobility variables are higly correlated
# may choose one of them and fit our model

###############################################
#Similarly, we check the correlation between following variables:
#income_support
#debt_relief
#fiscal measures
#testing_policy
#contact_tracing
#healthcare_investment
#masks
#test_new_smoothed
#ignore public_information, which is included as part of the variable "goverment_response_stringency_index"
#ignore also international_support which I think does not affect domestic transmission
mymatrix <- as.matrix(na.omit(globalraw[,13:19]),ncol=7)
normalizedData <- mynormalize(mymatrix)
res <-rcorr(normalizedData,type="pearson")
flattenCorrMatrix(res$r, res$P)
#after excluding statistically significant correlated variables,
#we decide to retain these variables:
#debt_relief
#fiscal measures
#contact_tracing
#healthcare_investment

###############################################
#Similarly, we check the correlation between demographic variables
demodf <- with(globalraw,cbind(stats_population,stats_population_density,stats_median_age,stats_population_urban,stats_population_school_age))
demomat <- as.matrix(na.omit(demodf),ncol=5)
normalizedData <- mynormalize(demomat)
res <-rcorr(normalizedData,type="spearman")
flattenCorrMatrix(res$r, res$P)
#all of these variables are significantly correlated
#we decide to retain these variables in the model and 
#conduct model selection to retain one or two of them in the best model
#population density 
#population
#median_age

###############################################
#After a preliminary analysis above,
#we decide to select a subset of variables from the raw dataset.
dateFrom <- as.integer(gsub("-", "","2020-03-01"))
dateTo <-as.integer(gsub("-", "","2020-12-31"))
library("dplyr")
library("aweek")
library("zoo")
globalsub <- globalraw %>%
                mutate(date = as.integer(gsub("-", "",DATE))) %>%
                mutate(cases_new_smoothed = zoo::rollmean(cases_new, k = 7, fill = NA))%>%
                filter(dateFrom <= date & date <= dateTo) %>%
                mutate(weekno = date2week(DATE, numeric = TRUE))  %>%
                 select(country_name,DATE,weekno,
                      cases_new_smoothed, 
                      cases_days_since_first,
  
                       stats_population, 
                       stats_population_density,
                       stats_median_age,
                       
                       stats_gdp_per_capita,
                       npi_stringency_index,
                       tests_new_smoothed,
                       mobility_retail_recreation,
                       npi_debt_relief,
                       npi_contact_tracing,
                       npi_healthcare_investment,
                      
                       weather_precipitation_mean,
                       weather_humidity_mean,
                       weather_sw_radiation_mean, 
                       weather_temperature_mean,
                       weather_wind_speed_mean)



#remove rows that contain NA and then check if the final dataframe contains NA values
colSums(is.na(globalraw))
globalsub <- na.omit(globalsub) 
apply(globalsub, 2, function(x) any(is.na(x)))


#calculate weekly mean of the variables
globalweek <-    with(globalsub,
                      aggregate(cbind(cases_new_smoothed, 
                                cases_days_since_first,
                                
                                stats_population, 
                                stats_population_density,
                                stats_median_age,
                                
                                stats_gdp_per_capita,
                                npi_stringency_index,
                                tests_new_smoothed,
                                mobility_retail_recreation,
                                npi_debt_relief,
                                npi_contact_tracing,
                                npi_healthcare_investment,
                                
                                weather_precipitation_mean,
                                weather_humidity_mean,
                                weather_sw_radiation_mean, 
                                weather_temperature_mean,
                                weather_wind_speed_mean), 
                                by=list(country_name,as.factor(weekno)), mean,na.rm=TRUE))%>% 
                                rename(country = Group.1, weekno = Group.2,
                                       stringency_index = npi_stringency_index,
                                       debt_relief = npi_debt_relief,
                                       contact_tracing = npi_contact_tracing,
                                       healthcare_investment = npi_healthcare_investment) %>% 
                                arrange(country,weekno) 
                                
table(globalweek$country)


##############################################################
##create data on the effect-lagged variables
##############################################################   
#produce dataframes for non lagged, one week lagged and two weeks lagged variables

### create indeices for the effect lagged variables

#index for no lagged effect variables
index <- rep(0,nrow(globalweek))
for(i in 2: (nrow(globalweek)-1)){
  if(!identical(globalweek$country[i],globalweek$country[i-1])){
    index[i] = 1
    index[i+1]= 1}
}
index[1]= index[2]=1

#index for one week lagged effect 
index_lag1 <- rep(0,nrow(globalweek))
for(i in 1: (nrow(globalweek)-1)){
  if(!identical(globalweek$country[i],globalweek$country[i+1])){
    index_lag1[i] = 1
    index_lag1[i+1]= 1}
}
index_lag1[1]=1
index_lag1[nrow(globalweek)]=1

#index for two weeks lagged effect 
index_lag2 <-  rep(0,nrow(globalweek))
for(i in 2: (nrow(globalweek)-1)){
  if(!identical(globalweek$country[i],globalweek$country[i+1])){
    index_lag2[i-1] = 1
    index_lag2[i] = 1
  }
}
index_lag2[nrow(globalweek)-1]=1
index_lag2[nrow(globalweek)]=1

globalweek_nolag <- globalweek %>% filter(index == 0) 
globalweek_lag1 <- globalweek[index_lag1==0,]
colnames(globalweek_lag1) <- paste(colnames(globalweek_nolag), "lag1", sep = "_")
globalweek_lag2 <- globalweek[index_lag2==0,]
colnames(globalweek_lag2) <- paste(colnames(globalweek_nolag), "lag2", sep = "_")

weekDfnew <-  cbind(globalweek_nolag,globalweek_lag1,globalweek_lag2) %>%
                     arrange(country,weekno) %>%
                     select(-c("country_lag1","country_lag2","weekno_lag1","weekno_lag2"))

## We exclude countries with too few data points (<5)
## and countries with poor data qualities by visual inspection
# for Japan we exclude data before march due to the outbreak 
# in cruise Diamand
countryIndx <- as.data.frame(table(weekDfnew$country)>=5)
weekDfnew <- weekDfnew %>%
             filter(country %in% rownames(countryIndx)[countryIndx$`table(weekDfnew$country) >= 5`]) %>%
             filter(!(country %in% c("Fiji","Oman"))) %>%
             filter(!(country =="Japan" & as.numeric(weekno) <= 12))

levels(weekDfnew$country)[21] = "Bolivia"
levels(weekDfnew$country)[46] = "Ivory Coast"
levels(weekDfnew$country)[91] = "Korea"
levels(weekDfnew$country)[140] = "Russian"
levels(weekDfnew$country)[182] = "Vietnam"

#creat a new vector index to indicte if the error is normal distributed or the variance of which should 
#be caculated with autoregressive conditional heteroskedasticity
indexFit <- rep(0,nrow(weekDfnew))
indexFit[1] = 1
for(i in 2: (nrow(weekDfnew)-1)){
  if(identical(weekDfnew$country[i],weekDfnew$country[i-1])){
    indexFit[i] = 0
  }else{
    indexFit[i]= 1
  }
}
indexFit  <- c(0,indexFit[1:(nrow(weekDfnew)-1)]) + indexFit 

globaldata <- cbind(weekDfnew,indexFit)

range(globnorm$weather_precipitation_mean)



###### normalize the predictive variables
#mean normalization
mean_norm <- function(x){
  if(max(x) != min(x)) return((x- mean(x)) /(max(x)-min(x)))
  if(max(x) == min(x)) return(rep(0,length(x)))
}

globnorm <-  cbind(globaldata[,1:3],normalizefun(as.matrix(globaldata[,4:(ncol(globaldata)-1)])),globaldata[,ncol(globaldata)]) %>%
             mutate(logcases_new_smoothed = log(cases_new_smoothed)) %>%
             rename("index" = "globaldata[, ncol(globaldata)]")

globnorm[globnorm == "-Inf"] <- 0


### creat TMB formated data list
attach(globnorm)
tmb_newdata =
  list(
   lognewcasessmooth = exp(logcases_new_smoothed),
    
    temp =  weather_temperature_mean,
    temp_sq =  weather_temperature_mean^2,
    temp_lag1 =  weather_temperature_mean_lag1,
    temp_lag2 = weather_temperature_mean_lag2,

    windsp = weather_wind_speed_mean,
    windsp_lag1 = weather_wind_speed_mean_lag1,
    windsp_lag2 = weather_wind_speed_mean_lag2,

    humidity = weather_humidity_mean,
    humidity_lag1 = weather_humidity_mean_lag1,
    humidity_lag2 = weather_humidity_mean_lag2,
    
    uv = weather_sw_radiation_mean,
    uv_lag1 = weather_sw_radiation_mean_lag1,
    uv_lag2 =  weather_sw_radiation_mean_lag2,

    population = stats_population,
    dens = stats_population_density,
    mage = stats_median_age,

    gdp = stats_gdp_per_capita,
  
    newtests = tests_new_smoothed,
    newtestssq = tests_new_smoothed^2,
    # 
    strindex = stringency_index,
    strindex_lag1 = stringency_index_lag1,
    strindex_lag2 = stringency_index_lag2,
    # 
    mobility = mobility_retail_recreation,
    mobility_lag1 = mobility_retail_recreation_lag1,
    mobility_lag2 = mobility_retail_recreation_lag2,

    debtrelief = debt_relief,
    debtrelief_lag1 = debt_relief_lag1,
    debtrelief_lag2 = debt_relief_lag2,

    healthinvestment = healthcare_investment,
    healthinvestment_lag1 = healthcare_investment_lag1,
    healthinvestment_lag2 = healthcare_investment_lag2,

    contatracing = contact_tracing,
    contatracing_lag1 = contact_tracing_lag1,
    contatracing_lag2 = contact_tracing_lag2,
    
    cases_days_since_first =  cases_days_since_first,
    cases_days_since_first_sq = cases_days_since_first^2,
    
    spatioLevelfactor = factor(globnorm$country),
    index = globnorm$index
  )
detach(globnorm)


#########################################
#instead of mean-normalizing the covariates,
#we can also standardize the covariates
library("dplyr")
globscale <-  cbind(globaldata[,1:3],scale(as.matrix(globaldata[,4:(ncol(globaldata)-1)])),globaldata[,ncol(globaldata)]) %>%
               mutate(logcases_new_smoothed = log10(cases_new_smoothed)) %>%
               rename("index" = "globaldata[, ncol(globaldata)]")

globscale[globscale == "-Inf"] <- 0
  
attach(globscale)
tmb_scaledata =
  list(
    lognewcasessmooth = logcases_new_smoothed,
    #newcasessmooth = cases_new_smoothed,
    
    temp =  weather_temperature_mean,
    temp_sq =  weather_temperature_mean^2,
    temp_lag1 =  weather_temperature_mean_lag1,
    temp_lag2 = weather_temperature_mean_lag2,
    
    windsp = weather_wind_speed_mean,
    windsp_lag1 = weather_wind_speed_mean_lag1,
    windsp_lag2 = weather_wind_speed_mean_lag2,
    
    humidity = weather_humidity_mean,
    humidity_lag1 = weather_humidity_mean_lag1,
    humidity_lag2 = weather_humidity_mean_lag2,
    
    uv = weather_sw_radiation_mean,
    uv_lag1 = weather_sw_radiation_mean_lag1,
    uv_lag2 =  weather_sw_radiation_mean_lag2,
    
    population = stats_population,
    dens = stats_population_density,
    mage = stats_median_age,
    
    gdp = stats_gdp_per_capita,
    
    newtests = tests_new_smoothed,
    newtestssq = tests_new_smoothed^2,
    # 
    strindex = stringency_index,
    strindex_lag1 = stringency_index_lag1,
    strindex_lag2 = stringency_index_lag2,
    # 
    mobility = mobility_retail_recreation,
    mobility_lag1 = mobility_retail_recreation_lag1,
    mobility_lag2 = mobility_retail_recreation_lag2,
    
    debtrelief = debt_relief,
    debtrelief_lag1 = debt_relief_lag1,
    debtrelief_lag2 = debt_relief_lag2,
    
    healthinvestment = healthcare_investment,
    healthinvestment_lag1 = healthcare_investment_lag1,
    healthinvestment_lag2 = healthcare_investment_lag2,
    
    contatracing = contact_tracing,
    contatracing_lag1 = contact_tracing_lag1,
    contatracing_lag2 = contact_tracing_lag2,
    
    cases_days_since_first =  cases_days_since_first,
    cases_days_since_first_sq = cases_days_since_first^2,
    
    spatioLevelfactor = factor(globscale$country),
    index = globscale$index
  )
detach(globscale)   


#########################################
#add the previous week number of cases as an covarites in the model
library("dplyr")
globscale <-  cbind(globaldata[,1:3],scale(as.matrix(globaldata[,4:(ncol(globaldata)-1)])),globaldata[,ncol(globaldata)]) %>%
              mutate(logcases_new_smoothed = log10(cases_new_smoothed)) %>%
              rename("index" = "globaldata[, ncol(globaldata)]")

globscale[globscale == "-Inf"] <- 0

pweekcases <- as.vector(scale(globscale$logcases_new_smoothed[1:c(nrow(globscale)-1)]))

globscale_pweek <- cbind(globscale[2:nrow(globscale),],pweekcases)

attach(globscale_pweek)
tmb_scale_pweek =
  list(
    lognewcasessmooth = logcases_new_smoothed,
    pweekcases = pweekcases,
    
   # cases_new_smoothed = cases_new_smoothed,
    
    temp =  weather_temperature_mean,
    temp_sq =  weather_temperature_mean^2,
    temp_lag1 =  weather_temperature_mean_lag1,
    temp_lag2 = weather_temperature_mean_lag2,
    
    windsp = weather_wind_speed_mean,
    windsp_lag1 = weather_wind_speed_mean_lag1,
    windsp_lag2 = weather_wind_speed_mean_lag2,
    
    humidity = weather_humidity_mean,
    humidity_lag1 = weather_humidity_mean_lag1,
    humidity_lag2 = weather_humidity_mean_lag2,
    
    uv = weather_sw_radiation_mean,
    uv_lag1 = weather_sw_radiation_mean_lag1,
    uv_lag2 =  weather_sw_radiation_mean_lag2,
    
    population = stats_population,
    dens = stats_population_density,
    mage = stats_median_age,
    
    gdp = stats_gdp_per_capita,
    
    newtests = tests_new_smoothed,
    newtestssq = tests_new_smoothed^2,
    # 
    strindex = stringency_index,
    strindex_lag1 = stringency_index_lag1,
    strindex_lag2 = stringency_index_lag2,
    # 
    mobility = mobility_retail_recreation,
    mobility_lag1 = mobility_retail_recreation_lag1,
    mobility_lag2 = mobility_retail_recreation_lag2,
    
    debtrelief = debt_relief,
    debtrelief_lag1 = debt_relief_lag1,
    debtrelief_lag2 = debt_relief_lag2,
    
    healthinvestment = healthcare_investment,
    healthinvestment_lag1 = healthcare_investment_lag1,
    healthinvestment_lag2 = healthcare_investment_lag2,
    
    contatracing = contact_tracing,
    contatracing_lag1 = contact_tracing_lag1,
    contatracing_lag2 = contact_tracing_lag2,
    
    cases_days_since_first =  cases_days_since_first,
    cases_days_since_first_sq = cases_days_since_first^2,
    
    spatioLevelfactor = factor(globscale_pweek$country),
    index = globscale_pweek$index
  )
detach(globscale_pweek)   

df_tmb_pweek <- data.frame(matrix(unlist(tmb_scale_pweek), nrow=length(tmb_scale_pweek$lognewcasessmooth),
                            byrow=FALSE,dimnames = list(NULL,names(tmb_scale_pweek))),stringsAsFactors=FALSE)  %>%
                 mutate(spatioLevelfactor = as.factor(spatioLevelfactor)) 

#now we run the selected model with glmmTMB, by adding the previous week's cases, the model fit improves a lot
#So in the manual model seleciton procedure, we decide to include it as an covariate.
library("glmmTMB")
mbest_glmm_pweek <- glmmTMB(cases_new_smoothed ~  temp +  uv + population  + mage + pweekcases +
                        newtests + mobility + debtrelief  + contatracing + cases_days_since_first +
                        (1 | spatioLevelfactor) + (0 + mobility| spatioLevelfactor) +
                        (0 + cases_days_since_first | spatioLevelfactor) + 
                        (0 + cases_days_since_first_sq | spatioLevelfactor),
                        family=nbinom2, data=df_tmb_pweek)
check_overdispersion(mbest_glmm_pweek)
summary(mbest_glmm_pweek)
#check the residuals


hist(df_tmb_pweek$cases_new_smoothed)

mbest_lmer_pweek <- lmer(lognewcasessmooth ~  temp +  uv + population  + mage + pweekcases +
                     newtests + mobility + debtrelief  + contatracing + cases_days_since_first +
                     (1 | spatioLevelfactor) + (0 + mobility| spatioLevelfactor) +
                     (0 + cases_days_since_first | spatioLevelfactor) + 
                     (0 + cases_days_since_first_sq | spatioLevelfactor),
                     data=df_tmb_pweek)

summary(mbest_lmer_pweek)
AIC(mbest_lmer_pweek) 263.402
         
