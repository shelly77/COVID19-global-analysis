
### in this script, we use the one data source that integrates severall
# data sources. see https://rs-delve.github.io/data_software/global-dataset.html

# this script contains the code dealing with raw data and
# make a data list compatible with TMB, in which change in growth rate 
# is taken as the response variable.

#load the data file
globalraw  <- read.csv("rawData/combined_dataset_latest.csv", header=T)

###############################################
#After a preliminary analysis in datapreGlobal.R,
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


#calculate daily growth rate change for each country
n <- length(unique(globalsub$country_name))
grorate_change <- vector()
for(i in 1:n){
  country <- unique(globalsub$country_name)[i]
  count_one <- globalsub[globalsub$country_name==country,]
  size <- nrow(count_one)
  change_one_count <- 
    (count_one$cases_new_smoothed[2:size]-count_one$cases_new_smoothed[1:(size-1)])/count_one$cases_new_smoothed[1:(size-1)]
  grorate_change <- c(grorate_change,c(NA,change_one_count))
  }

# add calculated growth rate change to the dataset
globalsub <- globalsub  %>%
             mutate(grorate_change = grorate_change )

#remove rows that contain NA and then check if the final dataframe contains NA values
globalsub <- na.omit(globalsub) 
apply(globalsub, 2, function(x) any(is.na(x)))


#calculate weekly mean of the variables
globalweek <-    with(globalsub,
                      aggregate(cbind(grorate_change,
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
                                      weather_wind_speed_mean), 
                                by=list(country_name,as.factor(weekno)), mean,na.rm=TRUE))%>% 
                rename(country = Group.1, weekno = Group.2,
                stringency_index = npi_stringency_index,
                debt_relief = npi_debt_relief,
                contact_tracing = npi_contact_tracing,
                healthcare_investment = npi_healthcare_investment) %>% 
                arrange(country,weekno) 


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


###### normalize the predictive variables
#mean normalization
mean_norm <- function(x){
  if(max(x) != min(x)) return((x- mean(x)) /(max(x)-min(x)))
  if(max(x) == min(x)) return(rep(0,length(x)))
}

globnorm <-  cbind(globaldata[,1:3],normalizefun(as.matrix(globaldata[,4:(ncol(globaldata)-1)])),globaldata[,ncol(globaldata)]) %>%
             rename("index" = "globaldata[, ncol(globaldata)]")

globnorm[globnorm == "-Inf"] <- 0
globnorm[globnorm == "Inf"] <- 0
globnorm[globnorm == "NA"] <- 0


### creat TMB formated data list
attach(globnorm)
tmb_data_growrate =
  list(
    grorate_change =globnorm$grorate_change,
    
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




