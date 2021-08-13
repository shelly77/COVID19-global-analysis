#This file contains self-defined functions used to produce tables and figures in the paper.
#These functions are also applied to the best model and data to get the tables and figures in the paper.

library("stringr")

#load all the model objects fitted in this study and that saved as a big Rdata file
load("allModelsLastVersion.RData")

#save.image(file = "my_work_space.RData")
#####################################################################################
#To run the following codes, you have to run the code for model selection, get the
#candidate models fitted, and save the best model (modelT) object as "BestModel.Rdata".

#The estimates of fixed effects and random effects are accessible in "modelT"
#but the reported values (poi) are only stored in object "bestModeloneStep"
BestModel <- modelonestep

#load TMB format data
#load("data.RData")

#check on the estimates and reported values of the best model
rep_fixed <- summary(BestModel$sdreport,"fixed",p.value = TRUE)
rep_random <- summary(BestModel$sdreport,"random")
rep_value <- summary(BestModel$sdreport,"report")


#Calculate the correlation between random effects
mu_spatioLevel <- as.data.frame(rep_random[rownames(rep_random)=="mu_spatioLevel",])[1]
mu_spatioLevel_mobility <- as.data.frame(rep_random[rownames(rep_random)=="mu_spatioLevel_mobility",])[1]
mu_spatioLevel_days <- as.data.frame(rep_random[rownames(rep_random)== "mu_spatioLevel_days",])[1]
mu_spatioLevel_sqdays  <- as.data.frame(rep_random[rownames(rep_random)== "mu_spatioLevel_sqdays",])[1]
random_est <- as.data.frame((cbind(mu_spatioLevel,mu_spatioLevel_mobility,mu_spatioLevel_days,mu_spatioLevel_sqdays)))
library(Hmisc)
res2 <- rcorr(as.matrix(random_est[,1:4]))
flattenCorrMatrix(res2$r, res2$P)

#check the reported values
estimates_value <- summary(BestModel$sdreport,"report")
error_est <- subset(estimates_value, rownames(estimates_value )=="error")
state_est <- subset(estimates_value, rownames(estimates_value )=="state")

#plot the acf function for the residuals
acf(ts(errorAll[,1]))


#########Table of estimates from the selected model##################

###################################################
################## Table 1 ########################
###################################################

#Function to generate table of parameter estimates 
model <- modelbest

parametertableP <- function(partable,model,file=NULL) {
  if (!is.null(file)) 
    sink(file)
  cat("\\scalebox{1.2}{\n")
  cat("\\begin{tabular}{c|l |l} \n")
  cat(" \textbf{Parameter} &\textbf{Estimate($\pm$s.d.)}  &\textbf{Description} \\\\ \n")
  cat("\\midrule\n")
  
  catpm <- function(est) {
    if (is.na(est[1]))
      cat(" & $0$")
    else
      cat(" & $",format(est[1],nsmall=nsmall,scientific = nsmall),
          "\\pm",format(est[2],nsmall=nsmall,scientific = nsmall),"$")
  }
  
  estimates <- summary(model$sdreport,c("fixed","report"))
  for (parname in rownames(partable)) {      # loop through all possible parameters
    if (parname=="midrule") {
      cat("\\midrule\n")
    } else {
      parlength <- length(model$parameters[['beta0']])
      if (is.null(model$map[[parname]]))          # if not included in map argument 
        map <- 1:parlength                      
      else                                        # otherwise
        map <- model$map[[parname]]              # look at the map argument
      if (sum(!is.na(map))>0) {                   # if at least one level is estimated include in table
         if (substr(parname,1,3)=="log") {         # estimates through ADREPORT
              name <- substring(parname,4)          # estimates appear 
          subset <- rownames(estimates)==name
          est <- estimates[subset,,drop=FALSE]
         }
        else { # the parameter appears in fixed part once for each estimated level 
          name <- parname
          subset <- rownames(estimates)== parname
          if (is.null(model$map[[parname]])) est <- estimates[subset,,drop=FALSE]
          else est <- estimates[subset,,drop=FALSE][model$map[[parname]],,drop=FALSE]
        }
        est[est[,"Std. Error"]==0,] <- NA # Hack to make things work also for parameters that defaults to non-zero values
        digits <- max(0,2-ceiling(log10(min(est[,"Std. Error"], na.rm=TRUE))))
        est <- round(est,digits)
        nsmall <- max(0,digits)
        cat( "$\\hat",partable[parname,"symbol"],"$ ",sep="")
        catpm(est)
        cat("\\\\ \n")
      }
    }
  }
  
  cat("\\bottomrule\n")
  cat("\\end{tabular}}\n")
  ##cat("\}\n")
  cat("\\begin{tablenotes}\n")
  cat("\\small\n")
  cat("\\begin{itemize}\n")
  for (parname in rownames(partable)) {  
    cat("\\item[$-$]", "$\\hat",partable[parname,"symbol"],"$"," : ", partable[parname,"meaning"],"\n",sep="")}
  cat("\\end{itemize}\n")
  cat("\\end{tablenotes}\n")
  
  if (!is.null(file)) 
    sink()
}

partable <- read.table("parameterDescrip.txt",header=TRUE,sep=";",as.is=TRUE,strip.white=TRUE)
parametertableP(partable,BestModel,file=NULL) 

####################################################
################## Figure 1 ########################
####################################################
#plot the map
library(tidyverse)
library(ggplot2)
library(readr)
library(maps)
library(viridis)

world_cities <- read.csv("worldcities.csv", header=T) %>%
                filter(capital == "primary") %>%
                filter(country %in% globaldata$country)  %>%
                select(country, lat,lng)  %>%
                rename(label = country, x = lat, y = lng)

#setdiff(unique(globaldata$country),world_cities$country)

#Add the countries that not in "world_cities" to our data frame
countries_add <- data.frame(country = c("Bosnia and Herzegovina", "Ivory Coast", "Korea", "Myanmar","Russian","Trinidad and Tobago"), 
                            city =  c("Sarajevo", "Yamoussoukro", "Seoul", "Naypyitaw", "Moscow","Trinidad and Tobago"),  
                            lat=  c(43.9159, 7.5400, 37.5665, 19.7633,55.7558,10.6918 ),  
                            lng =  c(17.6791, -5.5471, 126.9780, 96.0785, 37.6173,-61.2225)) 

# Combine the dataframes contaning location data
all_countries = rbind(world_cities,countries_add)

map_cases <- as.data.frame(with(globaldata,
                              aggregate(cases_new_smoothed, 
                                        by=list(country), sum,na.rm=TRUE))) %>%
                 rename(country = Group.1, sumcases = x)

# Merge the dataframes contaning location data and number of cases data
map_data_all <- merge(all_countries,map_cases) %>%
                distinct(country,.keep_all = T) %>%
                mutate(log(Cases) = log(sumcases))
head(map_data_all)

## get the world map
world_location <- map_data("world") %>%
                 dplyr::rename(country = region)

world_data <- merge(world_location,map_cases, by = "country")
nrow(world_data)


log10(sumcases) <- -sort(log(map_data_all$sumcases))

ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat,fill=id, group = group), fill="grey") +
  geom_point(data=map_data_all, aes(x=lng, y=lat,size=log10(sumcases), color=log10(sumcases) ,stroke=F)) +
  scale_colour_gradient(high="dark red",low="yellow")

 # for (i in 1:nrow(map_data_all)){annotate("text",x =map_data_all$lat, y =map_data_all$lng, label = map_data_all$country,size= 3 )} +
  #geom_text(aes(label = world_cities), size = 4) +
  #scale_size_continuous(name="Cases", range=c(1,7),breaks=mybreaks, labels = c("1-999", "1,000-9,999", "10,000+")) +
  # scale_alpha_continuous(name="Cases", trans="log", range=c(0.1, 0.9),breaks=mybreaks) +
  #scale_color_viridis_c(option="inferno",name="Cases", trans="log", breaks=mybreaks, labels = c("1-999",  "1,000-9,999", "10,000+")) +
  theme_void()+
  #guides( colour = guide_legend(),fill = FALSE)+

 # labs(caption = "The number of cumulatative cases used in the model.") +
  theme(
    legend.position = "right",
    text = element_text(color = "black"),
    plot.background = element_rect(fill = "#ffffff", color = NA), 
    panel.background = element_rect(fill = "#ffffff", color = NA), 
    legend.background = element_rect(fill = "#ffffff", color = NA)
  )


####################################################
################## Figure 2 ########################
####################################################
###############################################
# We plot the correlation matrix for the retained variables

# Mark the insignificant coefficients according to the specified p-value significance level
head(globscale_pweek)

globscale_rename <- globscale_pweek %>%
                    rename(previous_week_cases = pweekcases)

head(globscale_rename)
datamat <- scale(as.matrix(globscale_rename[,c(4:8,10:19,56)])) 

require(Hmisc)
require(corrplot)
corvar <- rcorr(datamat,type = "pearson")
M <- corvar$r
p_mat <- corvar$P
col <- colorRampPalette(c("darkorange", "white", "steelblue"))
corrplot(M, method = "color", col = col(200),  
         type = "upper", order = "hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "darkblue", tl.srt = 45, #Text label color and rotation
         # Combine with significance level
         p.mat = p_mat, sig.level = 0.02,  
         # hide correlation coefficient on the principal diagonal
         diag = FALSE 
)

####################################################
################## Figure 3 ########################
####################################################
#plot the pie chart for all the effects

library(plotrix)
par(mar=c(0.1,1.2,2.2,0.1), oma=c(-2, -2, 1, 1),mgp=c(3,0.4,0))
# browser data without "ymax" and "ymin"
browsers <-
  structure(
    list(
      browser = structure(
        c(1L, 2L, 2L, 2L, 3L ),
        #c(1L,2L, 3L,4L,5L),
        .Label = c( "Fixed effects" ,paste("Random","effects", sep = "\n"), "Unexplained"),
        class = "factor"
      ),
      version = structure(
        c(1L,2L, 3L, 4L,5L),
       # c(1L, 1L, 1L, 2L, 2L, 3L,
       #   3L, 3L, 4L, 4L,4L,5L),
        .Label = c(
          "Fixed effects",
          "Random intercepts",
          "Random mobility effects",
          "Random temporal effects",
          "Unexplained"
        ),
        class = "factor"
      ),
      share = c(21.6, 42,2.6, 35, 1.2)
    ),
    .Names = c("parent", "node", "size")
    ,
    row.names = c(NA,-12L),
    class = "data.frame"
  )

# aggregate data for the browser pie chart
browser_data <-
  aggregate(browsers$size,
            by = list(browser = browsers$parent),
            FUN = sum)

# order version data by browser so it will line up with browser pie chart
version_data <- browsers[order(browsers$node), ]

browser_colors <- c("#7ae774","#747ae7", "chocolate1")

# adjust these as desired (currently colors all versions the same as browser)
version_colors <-
  c(
    "#7ae774",
    "#747ae7",
    "#747ae7",
    "#747ae7",
    "chocolate1"
  )

# format labels to display version and % market share
version_labels <- paste(version_data$node, paste(version_data$size, "%"), sep = "\n")

# coordinates for the center of the chart
center_x <- 0.5
center_y <- 0.5

plot.new()

# draw version pie chart first
version_chart <-
  floating.pie(
    xpos = center_x,
    ypos = center_y,
    x = version_data$size,
    radius = 0.35,
    border = "white",
    col = version_colors
  )

# add labels for version pie chart
pie.labels(
  x = center_x,
  y = center_y,
  angles = version_chart,
  labels = version_labels,
  radius = 0.38,
  bg = NULL,
  cex = 0.55,
  font = 2,
  col = "gray40"
)

# overlay browser pie chart
browser_chart <-
  floating.pie(
    xpos = center_x,
    ypos = center_y,
    x = browser_data$x,
    radius = 0.25,
    border = "white",
    col = browser_colors
  )

paste(browser_data$browser,browser_data$x, sep = " ")

# add labels pie chart
pie.labels(
  x = center_x,
  y = center_y,
  angles = browser_chart,
  labels = paste(browser_data$browser,paste(browser_data$x,"%"), sep = "\n"),
  radius = 0.125,
  bg = NULL,
  cex = 0.65,
  font = 2,
  col = "black"
)

title(main="Variance decomposition for all the effects",line=-2)


####################################################
################## Figure 3 ########################
####################################################
#plot the pie chart for the fixed effects
library(plotrix)
par(mar=c(0.1,1.2,2.2,2), oma=c(-2, -2, 1, 1),mgp=c(3,0.4,0))
# browser data without "ymax" and "ymin"
browsers <-
  structure(
    list(
      browser = structure(
        c(1L,1L,2L,2L,3L,3L,3L,3L,4L,4L),
        .Label = c("Climate", "Demography", "Disease control", "Time"),
        class = "factor"
      ),
      version = structure(
        c(1L,2L, 3L, 4L,5L,
          6L, 7L, 8L,9L,10L),
        # c(1L, 1L, 2L, 2L, 3L,
        #   3L, 3L, 4L, 4L,4L,5L),
        .Label = c(
          "UV radiation",
          "Temperature",
          "Population",
          "Medium age",
          "New tests",
          "Population mobility",
          paste("Contact","tracing", sep = "\n"),
          paste("Debt","relief", sep = "\n"),
          "Days since first case",
          paste("No. of cases","in previous week", sep = "\n")
        ),
        class = "factor"
      ),
      share = c( 0.87,0.2, 0.4, 0.19, 0.34, 3.4, 0.5, 0.3, 2.4,13 )
    ),
    .Names = c("parent", "node", "size")
    ,
    row.names = c(NA,-12L),
    class = "data.frame"
  )

# aggregate data for the browser pie chart
browser_data <-
  aggregate(browsers$size,
            by = list(browser = browsers$parent),
            FUN = sum)

# order version data by browser so it will line up with browser pie chart
version_data <- browsers[order(browsers$node), ]

browser_colors <- c("#b3e773", '#e7e273', '#7ce9e3',"aquamarine3")

# adjust these as desired (currently colors all versions the same as browser)
version_colors <-
  c(
    "#b3e773",
    "#b3e773",
    '#e7e273',
    '#e7e273',
    '#7ce9e3',
    '#7ce9e3',
    '#7ce9e3',
    '#7ce9e3',
    "aquamarine3",
    "aquamarine3"
  )

# format labels to display version and % market share
version_labels <- paste(version_data$node, paste(version_data$size, "%"), sep = "\n")

# coordinates for the center of the chart
center_x <- 0.5
center_y <- 0.5

plot.new()

# draw version pie chart first
version_chart <-
  floating.pie(
    xpos = center_x,
    ypos = center_y,
    x = version_data$size,
    radius = 0.35,
    border = "white",
    col = version_colors
  )

# add labels for version pie chart
pie.labels(
  x = center_x,
  y = center_y,
  angles = version_chart,
  labels = version_labels,
  radius = 0.38,
  bg = NULL,
  cex = 0.55,
  font = 2,
  col = "gray40"
)

# overlay browser pie chart
browser_chart <-
  floating.pie(
    xpos = center_x,
    ypos = center_y,
    x = browser_data$x,
    radius = 0.25,
    border = "white",
    col = browser_colors
  )

paste(browser_data$browser,browser_data$x, sep = " ")

# add labels for browser pie chart
pie.labels(
  x = center_x,
  y = center_y,
  angles = browser_chart,
  labels = paste(browser_data$browser,paste(browser_data$x,"%"), sep = "\n"),
  radius = 0.125,
  bg = NULL,
  cex = 0.65,
  font = 2,
  col = "black"
)

title(main="Variance decomposition for fixed effects",line=-5)

####################################################
################## Variance decomposition ##########
####################################################
#Calculate the variance explained by each covariate
#We use third measure in Xu, 2003

BestModel <- modelbest
#maic=  208.6521 
#p = 16

#Function to calculate the R^2 (equation 20 in Xu, 2003)
R_squared <- function(obj_reduc,obj_full){
            N <- length(tmb_scale_pweek$lognewcasessmooth)
            1-exp(-2/N*(obj_reduc-obj_full))
}

obj_full <- BestModel$opt$objective 

#refit the model in which one of the covariates are excluded from the best model
#We regard the best model as the full model here
require(TMB)
mod_temp <-  update(BestModel,optimizer="optim",list(beta_temp = factor(NA)),
                                                 
                    cv=0.01,seed=123)
mod_temp
mod_temp$sdreport
R_squared(mod_temp$opt$value,obj_full)
#maic= 209.9055 
#p = 15
#0.002


mod_uv <-  update(BestModel,optimizer="optim",list( beta_uv = factor(NA)),
                        
                        cv=0.01,seed=123)
mod_uv
mod_uv$sdreport
R_squared(mod_uv$opt$value,obj_full)
#maic=  219.4483 
#p = 15
# 0.0087

mod_population <-  update(BestModel,optimizer="nlminb",list(beta_population = factor(NA)),
                                                    
                    cv=0.01,seed=123)
mod_population 
mod_population$sdreport
R_squared(mod_population$opt$objective,obj_full)
#maic=  212.8352 
#p = 15
#0.004

mod_mage <-  update(BestModel,optimizer="nlminb",list(beta_mage = factor(NA)),
                          
                          cv=0.01,seed=123)
mod_mage 
mod_mage$sdreport
R_squared(mod_mage$opt$objective,obj_full)
#maic=  209.365 
#p = 15
#  0.0019


mod_newtests<-  update(BestModel,optimizer="nlminb",list(beta_newtests = factor(NA)),
                                                   
                    cv=0.01,seed=123)
mod_newtests
mod_newtests$sdreport
R_squared(mod_newtests$opt$objective,obj_full)
#maic=   211.5586  
#p = 15
#0.0034


mod_contatracing <-  update(BestModel,optimizer="nlminb",list(beta_contatracing = factor(NA)),
                            
                            cv=0.01,seed=123)
mod_contatracing 
mod_contatracing$sdreport
R_squared(mod_contatracing$opt$objective,obj_full)
#maic=   214.0228 
#p = 15
# 0.005


mod_debtrelief <-  update(BestModel,optimizer="nlminb",list(beta_debtrelief = factor(NA)),
                            
                            cv=0.01,seed=123)
mod_debtrelief
mod_debtrelief$sdreport
R_squared(mod_debtrelief$opt$objective,obj_full)
#maic=  211.0469 
#p = 15
# 0.003

mod_debtrelief <-  update(BestModel,optimizer="nlminb",list(beta_debtrelief = factor(NA)),
                          
                          cv=0.01,seed=123)
mod_debtrelief
mod_debtrelief$sdreport
R_squared(mod_debtrelief$opt$objective,obj_full)
#maic=  211.0469 
#p = 15
# 0.003

mod_pweekcases <-  update(BestModel,optimizer="nlminb",list(beta_pweekcases = factor(NA)),
                         cv=0.01,seed=123)
mod_pweekcases
mod_pweekcases$sdreport
R_squared(mod_pweekcases$opt$objective,obj_full)
#maic=  413.1346 
#p = 15
# 0.13

mod_mobility_fixed<-  update(BestModel,optimizer="nlminb",list(beta_mobility = factor(NA)),
                                                                          cv=0.01,seed=123)
mod_mobility_fixed
mod_mobility_fixed$sdreport
R_squared(mod_mobility_fixed$opt$objective,obj_full)
#maic=   256.983 
#p = 15
# 0.034

mod_days_fixed<-  update(BestModel,optimizer="nlminb",list(beta_days = factor(NA)),
                             cv=0.01,seed=123)
mod_days_fixed
mod_days_fixed$sdreport
R_squared(mod_days_fixed$opt$objective,obj_full)
#maic= 241.3037 
#p = 15
# 0.024

#remove random effects
mod_random_intercept <-  update(BestModel,optimizer="nlminb",list(logsd_spatioLevel  = factor(NA),
                                                                  mu_spatioLevel = factor(rep(NA,nlevels(tmb_scale_pweek$spatioLevel)))),
                                
                                
                                cv=0.01,seed=123)
mod_random_intercept
mod_random_intercept$sdreport
R_squared(mod_random_intercept$opt$objective,obj_full)
#maic=  1004.773 
#p = 15
#  0.42

mod_mobility_random <-  update(BestModel,optimizer="nlminb",list(
                                                         mu_spatioLevel_mobility = factor(rep(NA,nlevels(tmb_scale_pweek$spatioLevel))),
                                                         logsd_randm_slop_mobility = factor(NA)),
                       
                       cv=0.01,seed=123)
mod_mobility_random
mod_mobility_random$sdreport
R_squared(mod_mobility_random$opt$objective,obj_full)
#maic=  244.3958 
#p = 15
#  0.026


mod_days_random <-  update(BestModel,optimizer="nlminb",list(
                                                      logsd_randm_slop_days  = factor(NA),
                                                      mu_spatioLevel_days = factor(rep(NA,nlevels(tmb_scale_pweek$spatioLevel))),
                                                      logsd_randm_slop_sqdays  = factor(NA),
                                                      mu_spatioLevel_sqdays = factor(rep(NA,nlevels(tmb_scale_pweek$spatioLevel)))),
                       
                       cv=0.01,seed=123)
mod_days_random
mod_days_random$sdreport
R_squared(mod_days_random$opt$objective,obj_full)
#maic=   826.6525 
#p = 14
# 0.35


####################################################
################## Figure 5a ########################
####################################################
#plot the estimated random  slope of moility for each country
#Check the names of reported random effects
library("ggplot2")
# unique(rownames(rep_random))

random_slop_mobility <- as.data.frame(rep_value[rownames(rep_value)=="beta_mobility+mu_spatioLevel_mobility",])
random_slop_mobility$Estimate[order(random_slop_mobility$Estimate)]
# Convert the random effects data into a dataframe that can be fed to ggplot
random_slop_mobility <- as.data.frame(rep_random[rownames(rep_random)=="mu_spatioLevel_mobility",])
rownames(random_slop_mobility) <- unique(tmb_scale_pweek$spatioLevelfactor)
colnames(random_slop_mobility) <- c("Estimate", "SE")
random_slop_mobility$Country <- rownames(random_slop_mobility)
# 
#beta_mobility_fixed <- rep_fixed["beta_mobility","Estimate"]
random_slops <- random_slop_mobility %>%
                #mutate(Estimate = Estimate + beta_mobility_fixed) %>%
                arrange(Estimate)  %>%
                mutate(Country = as.character(Country))
#  
level_order <-random_slops$Country #this vector is used for reorder the countries according to the effect sizes
# 
#ggplot
plot_random_mobility <- ggplot(random_slops, aes(x=factor(Country,level = level_order), Estimate)) +
         geom_pointrange(aes(ymin = Estimate - 1.96*SE, ymax = Estimate + SE),color="darkblue") +
         geom_hline(yintercept=0,color = "red")  + coord_flip() +
         theme_bw(base_size=10)+
         labs(y = "Random effects of population mobility", x = NULL)+
         theme(
         axis.title.x = element_text(color = "black", size = 12, face = "bold"),
         axis.title.y = element_text(color = "black", size = 12, face = "bold")
         )

#plot the estimated random effect plus fixed effect reported by the model directly
# random_slop_mobility <- as.data.frame(rep_value[rownames(rep_value)=="beta_humidity+mu_spatioLevel_mobility",])
# rownames(random_slop_mobility) <- unique(tmb_scale_pweek$spatioLevelfactor)
# colnames(random_slop_mobility) <- c("Estimate", "SE")
# random_slop_mobility$Country <- rownames(random_slop_mobility)
# random_slops <- random_slop_mobility %>%
#   
# library("dplyr")
# random_slops <-  arrange(as.data.frame(random_slop_mobility), Estimate) %>%  
#                 mutate(Country = as.character(Country))  
#      
# level_order <-random_slops$Country 
# 
# plot_random_mobility <- ggplot(random_slops, aes(x=factor(Country,level = level_order), Estimate)) +
#   geom_pointrange(aes(ymin = Estimate - 1.96*SE, ymax = Estimate + SE),color="darkred") + 
#   geom_hline(yintercept=0,color = "red")  + coord_flip() +
#   theme_bw(base_size=10)+
#   labs(y = "Effect of population mobility", x = "Country")+
#   theme(
#     axis.title.x = element_text(color = "black", size = 10, face = "bold"),
#     axis.title.y = element_text(color = "black", size = 10, face = "bold")
#   )

####################################################
################## Figure 5b ########################
####################################################
#### plot observed number of cases and mobility on the same plot

par(mfrow=c(9,10))
par(mar=c(1,1,1.5,1), oma=c(1, 0.1,0.1, 0.2),
    mgp=c(4,0.3,0))
for (i in 1:length(unique(globaldata$country))){

  obsData <- globaldata[which(as.integer(tmb_scale_pweek$spatioLevelfactor)==i),]
  log_cases <- log(obsData$cases_new_smoothed)
  log_cases[log_cases == "-Inf"] <- 0
  title <-as.character(unique(globaldata$country)[i])
  #plot(as.numeric(obsData$weekno),log(obsData$cases_new_smoothed),pch=20,col = "darkblue", cex.axis=0.75)
  #if (i %in% c(c(0:10)*8+1)){mtext(side=2, line=1.3, "Log(cases)", font=1,cex=0.6)}
 # title(title, line = 0.2, cex.main=0.8)
 # lines(obsData$weekno,mean_norm(obsData$stringency_index),col="red")
  # polygon(c(rev(obsData$weekno), obsData$weekno), c(rev(estData[,1]+1.96*estData[,2]), estData[,1]-1.96*estData[,2]), col = rgb(0.7,0.7,0.7,0.7) , border = NA)
  # lines(obsData$weekno,estData[,1]+1.96*estData[,2],col="red",lty=2)
  #lines(obsData$weekno,estData[,1]-1.96*estData[,2],col="red",lty=2)
  plot.new()
  plot.window(xlim=range(as.numeric(obsData$weekno)), ylim=range(log_cases))
  points(as.numeric(obsData$weekno), log_cases, col="darkblue", type="l")
  axis(1,cex.axis=0.5)
  axis(2, col.axis="darkblue",cex.axis=0.5)
  box()
  title(title, line = 0.2, cex.main=0.8)
  
  if (i %in% c(81:87)){mtext(side=1, line=1.1, "Week", font=1,cex=0.6)}
  
  plot.window(xlim=range(as.numeric(obsData$weekno)), ylim=c(-1,1))
  points(as.numeric(obsData$weekno), mean_norm(obsData$mobility_retail_recreation), col="darkred", type="l")
  axis(4, col.axis="darkred",cex.axis=0.5)
}


####################################################
################## Figure 6 ########################
####################################################
#### plot estimated log(cases) and observed number of cases on the same plot

estAll <-as.data.frame(cbind(subset(estimates_value, rownames(estimates_value)=="state"),tmb_scale_pweek$spatioLevelfactor))
#unique(data$spatioLevelfactor)
par(mfrow=c(11,8))
par(mar=c(0.1,1.2,2.2,0.1), oma=c(6, 7, 1, 1),mgp=c(3,0.4,0))
for (i in 1:length(unique(tmb_scale_pweek$spatioLevelfactor))){
  obsData <- globscale_pweek[which(as.integer(tmb_scale_pweek$spatioLevelfactor)==i),]
  estData <- estAll[which(as.integer(tmb_scale_pweek$spatioLevelfactor)==i),]
  title <-as.character(unique(tmb_scale_pweek$spatioLevelfactor)[i])
  plot(as.numeric(obsData$weekno),log10(obsData$cases_new_smoothed),pch=20,col = "darkblue", cex.axis=0.75)
  if (i %in% c(c(0:10)*8+1)){mtext(side=2, line=1.3, "Log(cases)", font=1,cex=0.6)}
  if (i %in% c(81:87)){mtext(side=1, line=1.1, "Week", font=1,cex=0.6)}
  title(title, line = 0.2, cex.main=0.8)
  lines(obsData$weekno,estData[,1],col="red")
  polygon(c(rev(obsData$weekno), obsData$weekno), c(rev(estData[,1]+1.96*estData[,2]), estData[,1]-1.96*estData[,2]), col = rgb(0.7,0.7,0.7,0.7) , border = NA)
 # lines(obsData$weekno,estData[,1]+1.96*estData[,2],col="red",lty=2)
 # lines(obsData$weekno,estData[,1]-1.96*estData[,2],col="red",lty=2)
}


####################################################
################## Figure S1 ########################
####################################################
####################################################
#plot the estimated coefficient of days (the days after the first confirmed case) for each country

# Convert the estimated random effects into a dataframe that can be fed to ggplot
mu_spatioLevel_days <- as.data.frame(rep_random[rownames(rep_random)=="mu_spatioLevel_days",])

rownames(mu_spatioLevel_days) <-  unique(tmb_scale_pweek$spatioLevelfactor)
colnames(mu_spatioLevel_days) <- c("Estimate", "SE")
mu_spatioLevel_days$Country <- rownames(mu_spatioLevel_days)

beta_days_fixed <- rep_fixed["beta_days","Estimate"]
random_effects <- mu_spatioLevel_days %>%
  mutate(Estimate = Estimate + beta_days_fixed) %>%  
  arrange(Estimate)  %>%  
  mutate(Country = as.character(Country))  

level_order <- random_effects$Country #this vector is used for reorder the countries according to the effect sizes

#ggplot 
plot_days<- ggplot(random_effects, aes(x=factor(Country,level = level_order), Estimate)) +
  geom_pointrange(aes(ymin = Estimate - 1.96*SE, ymax = Estimate + 1.96*SE),color="darkblue") + 
  geom_hline(yintercept=0,color = "red") + coord_flip() +
  theme_bw(base_size=10)+
  labs(y = "Random effects of days (the number of days since first confirmed case)", x="Country")+
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold")
  )

#plot the estimated coefficient of squared days (the days after the first confirmed case) for each country

# Convert the estimated random effects into a dataframe that can be fed to ggplot
mu_spatioLevel_sqdays <- as.data.frame(rep_random[rownames(rep_random)=="mu_spatioLevel_sqdays",])

rownames(mu_spatioLevel_sqdays) <- unique(tmb_scale_pweek$spatioLevelfactor)
colnames(mu_spatioLevel_sqdays) <- c("Estimate", "SE")
mu_spatioLevel_sqdays$Country <- rownames(mu_spatioLevel_sqdays)

random_effects <- mu_spatioLevel_sqdays %>%
  arrange(Estimate)  %>%  
  mutate(Country = as.character(Country))  

level_order <- random_effects$Country #this vector is used for reorder the countries according to the effect sizes

#ggplot 
plot_sqdays<- ggplot(random_effects, aes(x=factor(Country,level = level_order), Estimate)) +
  geom_pointrange(aes(ymin = Estimate - 1.96*SE, ymax = Estimate + 1.96*SE),color="darkblue") + 
  geom_hline(yintercept=0,color = "red") + coord_flip() +
  theme_bw(base_size=10)+
  labs(y = "Random effects of quadratic days", x="Country")+
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold")
  )

####################################################
####################################################
#calculate the accumulated UV and cases for each country
uv_by_country <-  aggregate(cbind(globnorm$weather_sw_radiation_mean, globnorm$cases_new_smoothed),
                            by=list(globnorm$country),sum,na.rm=TRUE)



#############################################################
###########################################################
#predict the number of cases when there is no intervention
#and mobility keeps at average level (population mobility = 0)
#and number of new tests keep same
  nobs <- length(tmb_data$lognewcasessmooth)
  
  temp_lag1 <- tmb_data$temp_lag1
  windsp_lag2 <-  tmb_data$windsp_lag2
  humidity_lag2 <- tmb_data$humidity_lag2
  uv <- tmb_data$uv
  population <- tmb_data$population
  mage <- tmb_data$mage
  gdp <- tmb_data$gdp
  newtests <- tmb_data$newtests #tmb_data$newtests
  mobility <- rep(0,nobs)#tmb_data$mobility #tmb_data$mobility
  contatracing <-rep(-1,nobs)#tmb_data$contatracing
  days <- tmb_data$cases_days_since_first
  days_sq <- tmb_data$cases_days_since_first_sq
  
  beta0 <- BestModel$estimates$beta0
  beta_temp_lag1 <- BestModel$estimates$beta_temp_lag1
  beta_windsp_lag2  <- BestModel$estimates$beta_windsp_lag2 
  beta_humidity_lag2 <-   BestModel$estimates$beta_humidity_lag2 
  beta_uv <- BestModel$estimates$beta_uv 
  beta_population <- BestModel$estimates$beta_population  
  beta_mage <- BestModel$estimates$beta_mage
  beta_gdp <- BestModel$estimates$beta_gdp
  beta_newtests <- BestModel$estimates$beta_newtests
  beta_mobility  <- 0#BestModel$estimates$beta_mobility 
  beta_contatracing <-0 # BestModel$estimates$beta_contatracing
  beta_days <- BestModel$estimates$beta_days 
  beta_days_sq <- BestModel$estimates$beta_days_sq

  mu_spatioLevel <- unname(rep_random[rownames(rep_random)=="mu_spatioLevel",1])
  mu_spatioLevel_mobility <-rep(0,87) #unname(rep_random[rownames(rep_random)=="mu_spatioLevel_mobility",1])
  mu_spatioLevel_days <- unname(rep_random[rownames(rep_random)== "mu_spatioLevel_days",1])
  mu_spatioLevel_sqdays  <- unname(rep_random[rownames(rep_random)== "mu_spatioLevel_sqdays",1])
  
  times <- as.vector(unname(table(tmb_data$spatioLevelfactor)))

  est_no_control <- rep(beta0,nobs) + rep(mu_spatioLevel,times) + beta_temp_lag1*temp_lag1/
                    + beta_windsp_lag2*windsp_lag2 + beta_humidity_lag2*humidity_lag2 /
                    + beta_uv*uv + beta_population*population + beta_mage*mage/
                    + beta_newtests*i/
                    + beta_days*days /
                    + rep(mu_spatioLevel_days,times)*days + rep(mu_spatioLevel_sqdays,times)*days_sq
  
 
  #### plot estimated log(cases) without intervention and observed number of cases on the same plot    
  est_all <- as.data.frame(cbind(est_no_control, tmb_data$spatioLevelfactor)) %>%
            rename(country = V2)

  #unique(data$spatioLevelfactor)
  par(mfrow=c(11,8))
  par(mar=c(0.1,1.2,2.2,0.1), oma=c(6, 7, 1, 1),mgp=c(3,0.4,0))
  for (i in 1:length(unique(globnorm$country))){

    obsData <- globnorm[which(as.integer(tmb_data$spatioLevelfactor)==i),]
    estData <- est_all[which(as.integer(tmb_data$spatioLevelfactor)==i),]
    size <- nrow(estData)

    obs_y <- (log(obsData$cases_new_smoothed)[2:size]-log(obsData$cases_new_smoothed)[1:(size-1)])/
             log(obsData$cases_new_smoothed[1:(size-1)])
    
    est_y <- (estData[,1][2:size]-estData[,1][1:(size-1)])/
              estData[,1][1:(size-1)]
    obs_y[obs_y == "NaN"] <- 0
    obs_y[obs_y == "Inf"] <- 0
    est_y[est_y == "Inf"] <- 0
    title <-as.character(unique(globaldata$country)[i])
    plot(as.numeric(obsData$weekno)[2:size], obs_y,col = "darkblue", cex.axis=0.75, 
         ylim=range(c(est_y,obs_y)),type="l")
    if (i %in% c(c(0:10)*8+1)){mtext(side=2, line=1.3, "Log(cases)", font=1,cex=0.6)}
    if (i %in% c(81:87)){mtext(side=1, line=1.1, "Week", font=1,cex=0.6)}
    title(title, line = 0.2, cex.main=0.8)
    lines(obsData$weekno[2:size],est_y,col="red")
  }
  
  
#####PCA analysis
  
  #creat column "continent"for each study country
  library("dplyr")
  continents  <- read.csv("contries&continents.csv", header=T,sep = ";") 
  head(continents)
  
  countries <- as.data.frame(rownames(random_slop_mobility)) %>%
    dplyr::rename("Country.or.Area" = "rownames(random_slop_mobility)")
  
  count_contin <- merge(countries,continents, by = c("Country.or.Area"), all.x = T) %>%
                  select(Country.or.Area,Continent)
  
  write.csv(count_contin,"count_contin.csv", row.names = FALSE)
  
  count_contin_full  <- read.csv("country_continent_full.csv", header=T)   %>%
                       dplyr::rename("Country"= "Country.or.Area")

#perform pca analysis
  colnames(random_est) = c("random_intercept", "random_mobility","random_days","random_quadratic_days")
  
 randeffcts.pca <- prcomp(random_est, center = TRUE,scale. = TRUE)

#summary
summary(randeffcts.pca)

#load function "ggbiplot"
source("ggbiplot.R")

ggbiplot(randeffcts.pca, labels=count_contin_full$Country, groups = count_contin_full$Continent,ellipse=TRUE)

ggbiplot(randeffcts.pca, choices=c(1,3),labels=count_contin_full$Country, groups = count_contin_full$Continent,ellipse=TRUE)

ggbiplot(randeffcts.pca, choices=c(1,4),labels=count_contin_full$Country, groups = count_contin_full$Continent,ellipse=TRUE)

ggbiplot(randeffcts.pca, choices=c(2,3),labels=count_contin_full$Country, groups = count_contin_full$Continent,ellipse=TRUE)

ggbiplot(randeffcts.pca, choices=c(2,4),labels=count_contin_full$Country, groups = count_contin_full$Continent,ellipse=TRUE)

ggbiplot(randeffcts.pca, choices=c(3,4),labels=count_contin_full$Country, groups = count_contin_full$Continent,ellipse=TRUE)











