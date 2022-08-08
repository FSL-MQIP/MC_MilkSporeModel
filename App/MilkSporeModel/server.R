####################################### 1. Set up environment ###########################################
# load packages
library(shiny)
library(readr)
library(tidyverse)
library(censReg)
library(fitdistrplus)
library(splitstackshape)
library(rmutil)
library(truncnorm)
library(EnvStats)
library(epiR)
library(stats)

# load functions
source("UtilityFunctions.R")

# set seed
set.seed(1)


######################################## 2. Prepare input files ##########################################
# load input files
## (a) AT/ST frequency data 
spore_ATfrequency_file <- "InputFiles/ColdgrowATFreq_CrossSect_NGATs.csv"
spore_ATfreq_import <- read.csv(spore_ATfrequency_file, stringsAsFactors = FALSE, header = TRUE)

## (b) growth parameter data
spore_growth_file <- "InputFiles/GrowthParameters_NoGrowthATs.csv" 
spore_growth_import <-read.csv(spore_growth_file, stringsAsFactors = FALSE)

## (c) initial microbial count data 
spore_init_file <- read.csv("InputFiles/CrossSectional_RawMPN_3.6.17.csv")


# Modify input files
spore_ATfreq_import$temp <- "AT"
spore_ATfreq_import$ClosestAT <- paste(spore_ATfreq_import$temp,spore_ATfreq_import$ClosestAT,sep="_")

spore_growth_import <- spore_growth_import %>%
  .[c(1:5)]%>%
  rename(STorAT = rpoBAT)
spore_growth_import$model_name <- "buchanan" 
spore_growth_import$temp <- "AT"
spore_growth_import$STorAT <- paste(spore_growth_import$temp,spore_growth_import$STorAT,sep="_")


######################################## 3. Organize model parameters ####################################
# (a) Initial contamination: logMPN normal distribution
spore_init_file$log10left <- log10(spore_init_file$left)
spore_init_file$log10right <- log10(spore_init_file$right)
spore_init_file$log10MPN <- log10(spore_init_file$MPN)
cens_data <- spore_init_file[,c("log10left","log10right")]
names(cens_data) <- c("left","right")
spore_fit <- fitdistcens(censdata = cens_data,distr = "norm")

spore_log10MPN_mean <- spore_fit$estimate[1]
spore_log10MPN_sd <- spore_fit$estimate[2]

## (b) Frequency of allelic types
AT_freq <- spore_ATfreq_import$ClosestAT 

## (c) Simulation setting
n_sim <- 100 
n_units <- 10

lot_id <- rep(seq(1, n_sim), each = n_units)
unit_id <- rep(seq(1,n_units), times = n_sim ) 

STorAT <- vector(mode="logical", n_sim  *n_units)

#stage1, storage at facility
t_F <- vector(mode="logical", n_sim  *n_units) 
T_F <- vector(mode="logical", n_sim  *n_units) 
count_F <- vector(mode = "logical", n_sim  *n_units)

#stage2, transport to retail store
t_T <- vector(mode="logical", n_sim  *n_units) 
T_T <- vector(mode="logical", n_sim  *n_units) 
count_T <- vector(mode = "logical", n_sim  *n_units)

#stage3, storage/display at retail store
t_S <- vector(mode="logical", n_sim  *n_units) 
T_S <- vector(mode="logical", n_sim  *n_units) 
count_S <- vector(mode = "logical", n_sim  *n_units)

#stage4, transport from retail store to homes
t_T2 <- vector(mode="logical", n_sim  *n_units) 
T_T2 <- vector(mode="logical", n_sim  *n_units) 
count_T2 <- vector(mode = "logical", n_sim  *n_units)

#stage5, storage at homes
t_H <- vector(mode="logical", n_sim  *n_units) 
T_H <- vector(mode="logical", n_sim  *n_units) 
count_H <- vector(mode = "logical", n_sim  *n_units)

data <- data.frame(lot_id, unit_id,STorAT,t_F, T_F, count_F, 
                   t_T, T_T, count_T, t_S, T_S, count_S, t_T2, T_T2, count_T2,
                   t_H, T_H, count_H)

lagAtNewTemp <- function (newTemp, oldLag, oldTemp = 6, T0 = 1.15) {
  numerator <- oldTemp -T0
  denom <- newTemp - T0
  newLag <- ( (numerator / denom)^2) * oldLag
  return(newLag)
}

######################################### Model #############################################
shinyServer(function(input, output) {
  
#-----------------------------------------------------setting---------------------------------------------
  # get initial count
  
  model_result <- reactive({
    
  spore_log10MPN_samp <-  rnorm(n_sim, input$count_mean, input$count_sd) 

  #Convert spore_log10MPN_samp back to spore_MPN_samp
  spore_MPN_samp <- 10^(spore_log10MPN_samp) 
  
  #Convert MPN for each sample to equivalent in milk unit of interest (1900 mL in half gallon)
  spore_MPN_samp_halfgal <- spore_MPN_samp * 1900
  
  #Sample the MPN_init distribution
  spore_MPN_init<-vector()
  for (i in 1:n_sim){
    spore_MPN_init_samp <-rep(rpois(n_units, spore_MPN_samp_halfgal[i]))
    spore_MPN_init<-c(spore_MPN_init, spore_MPN_init_samp)}
  
  #First convert spore_MPN_init from half-gallon to mLs
  spore_MPN_init[spore_MPN_init<1] = 0
  spore_MPN_init_mL <- spore_MPN_init / 1900
  data$spore_MPN_init_mL <- spore_MPN_init_mL
  
  #Also need to remove 0's from the data and replace with detection limit
  data$spore_MPN_init_mL[data$spore_MPN_init_mL<=0 ] <- 0
  data$spore_MPN_init_mL[data$spore_MPN_init_mL == 0] <- 0.01
  
  #Add spore_log10MPN_init to dataframe
  data$spore_log10MPN_init_mL <- log10(data$spore_MPN_init_mL) 
  
  if (input$mf) {data$spore_log10MPN_init_mL = data$spore_log10MPN_init_mL - 2.2}
  if (input$bf1) {data$spore_log10MPN_init_mL = data$spore_log10MPN_init_mL - 1.4}
  if (input$bf2) {data$spore_log10MPN_init_mL = data$spore_log10MPN_init_mL - 2}


  
  #Sample AT from AT_freq & add to df; no-growth ATs were removed
  AT <- vector()
  for (i in 1:(n_units*n_sim)){
    AT_samp <- sample(AT_freq, 1,replace = T)
    while(AT_samp == "AT_23" || AT_samp == "AT_159"){
      AT_samp <- sample(AT_freq, 1,replace = T)}
    AT[i] = AT_samp
  }
  
  data$STorAT = AT
  
#-----------------------------------------------Stage 1: Storage at facility----------------------------------------------
  ## (a)  Sample the temperature distribution & add to dataframe
  temps_F <- rep(runif(n_sim,min=3.5,max=4.5),each=n_units) #uniform distribution
  if (input$f_reduceT) {temps_F <- rep(runif(n_sim,min=2.5,max=3.5),each=n_units)}
  if (input$f_supercool) {temps_F <- rep(runif(n_sim,min=0.5,max=1.5),each=n_units)}
  data$T_F <- temps_F
  
  ## (b) Sample the storage time (in days) distribution & add to dataframe
  times_F <- rep(runif(n_sim,min=1,max=2),each=n_units) #uniform distribution
  data$t_F <- times_F
  
  ## (c) Determine newLag_F and newMu_F
  for (i in 1:(n_sim*n_units)){
    #Find row in growth parameter data that corresponds to allele sample
    allele_index <- which(spore_growth_import$STorAT == data$STorAT[i])
    
    #Calculate the new growth parameters at new temp using the square root model 
    newLag_F <- lagAtNewTemp(data$T_F[i], spore_growth_import$lag[allele_index])
    newMu_F <-  muAtNewTemp(data$T_F[i], spore_growth_import$mumax[allele_index])
    
    data$newLag_F[i] <- newLag_F
    data$newMu_F[i] <- newMu_F
  }
  
  ## (d) Determine count_F
  for (i in 1:(n_sim *n_units)){
    #Find row in growth parameter data that corresponds to allele sample
    allele_index <- which(spore_growth_import$STorAT == data$STorAT[i]) 
    
    #Calculate the log10N count using the new growth parameters
    data$count_F[i] <- log10N_func(data$t_F[i], data$newLag_F[i],data$newMu_F[i],data$spore_log10MPN_init_mL[i],spore_growth_import$LOG10Nmax[allele_index])
  }
  
  
#-----------------------------------------------Stage 2: Transport from facility to retail-------------------------------
  ## (a)  Sample the temperature distribution & add to dataframe
  temps_T <- rep(rtri(n_sim,min=1.7,max=10.0,mode=4.4),each=n_units) #triangular distribution
  data$T_T <- temps_T
  
  ## (b) Sample the storage time (in days) distribution & add to df
  times_T <- rep(rtri(n_sim,min=1,max=10,mode=5),each=n_units)
  data$t_T <- times_T
  
  ## (c) Determine Lag_T and Mu_T
  ## determine lag & mumax at new temp
  for (i in 1:(n_sim*n_units)){
    #Find row in growth parameter data that corresponds to allele sample
    allele_index <- which(spore_growth_import$STorAT == data$STorAT[i])
    
    #Calculate the new growth parameters at new temp using the square root model 
    newLag_T <- lagAtNewTemp(data$T_T[i], spore_growth_import$lag[allele_index])
    newMu_T <-  muAtNewTemp(data$T_T[i], spore_growth_import$mumax[allele_index])
    
    data$newLag_T[i] <- newLag_T
    data$newMu_T[i] <- newMu_T
  }
  
  ### Determine Lag_T
  data$checkLagPhase_T <- ifelse((data$t_F/data$newLag_F) <1, 1, 0)
  data$adjLag_T <- ((1 - (data$t_F/data$newLag_F))*data$newLag_T) 
  data$Lag_T <- ifelse(data$checkLagPhase_T==0,0,data$adjLag_T)
  
  ### Determine Mu_T
  data$Mu_T <- ifelse(data$T_T >= data$T_F*0.75 & data$T_T <= data$T_F*1.25, data$newMu_F, data$newMu_T)
  data$checkTemp_FtoT <- ifelse(data$T_T >= data$T_F*0.75 & data$T_T <= data$T_F*1.25, "T_F", "T_T") 
  
  ## (d) Determine count_T
  for (i in 1:(n_sim *n_units)){
    #Find row in growth parameter data that corresponds to allele sample
    allele_index <- which(spore_growth_import$STorAT == data$STorAT[i]) 
    
    #Calculate the log10N count using the new growth parameters
    data$count_T[i] <- log10N_func(data$t_T[i], data$Lag_T[i],data$Mu_T[i],data$count_F[i],spore_growth_import$LOG10Nmax[allele_index])
  }
  
#-----------------------------------------Stage 3: Storage at retail store-------------------------------------------------------
  ## (a)  Sample the temperature distribution & add to dataframe
  unif_mean = 2.3
  unif_b = 5.4
  if (input$r_reduceT) {unif_mean = 1.8}
  if (input$r_alarm) {unif_b = 4}
  temps_S <- rep(rtruncnorm(n_sim,a=-1.4,b=unif_b,mean=unif_mean,sd=1.8),each=n_units) #triangular distribution
  data$T_S <- temps_S
  
  ## (b) Sample the retail storage & display time (in days) distribution & add to dataframe
  times_S <- rep(rtruncnorm(n_sim,a=0.042,b=10.0, mean=1.821,sd=3.3),each=n_units)
  data$t_S <- times_S
  
  ## (c) Determine Lag_S and Mu_S
  ## determine lag & mumax at new temp
  for (i in 1:(n_sim*n_units)){
    #Find row in growth parameter data that corresponds to allele sample
    allele_index <- which(spore_growth_import$STorAT == data$STorAT[i])
    
    #Calculate the new growth parameters at new temp using the square root model 
    newLag_S <- lagAtNewTemp(data$T_S[i], spore_growth_import$lag[allele_index])
    newMu_S <-  muAtNewTemp(data$T_S[i], spore_growth_import$mumax[allele_index])
    
    data$newLag_S[i] <- newLag_S
    data$newMu_S[i] <- newMu_S
  }
  
  ## Determine Lag_S
  data$checkLagPhase_S <- ifelse(((data$t_F/data$newLag_F)+(data$t_T/data$Lag_T)) <1, 1, 0)
  data$adjLag_S <- ((1 - (data$t_F/data$newLag_F)-(data$t_T/data$Lag_T))*data$newLag_S) 
  data$Lag_S <- ifelse(data$checkLagPhase_S==0,0,data$adjLag_S)
  
  ## Determine Mu_S
  data$Mu_S <- ifelse(data$T_S >= data$T_T*0.75 & data$T_S <= data$T_T*1.25, data$newMu_T, data$newMu_S)
  data$checkTemp_TtoS <- ifelse(data$T_S >= data$T_T*0.75 & data$T_S <= data$T_T*1.25, "T_T", "T_S") 
  
  ## (d) Determine count_S
  for (i in 1:(n_sim *n_units)){
    #Find row in growth parameter data that corresponds to allele sample
    allele_index <- which(spore_growth_import$STorAT == data$STorAT[i]) 
    
    #Calculate the log10N count using the new growth parameters
    data$count_S[i] <- log10N_func(data$t_S[i], data$Lag_S[i],data$Mu_S[i],data$count_T[i],spore_growth_import$LOG10Nmax[allele_index])
  }
  
#-------------------------------------------------Stage 4: Transport from retail to home--------------------------------------------------
  ## (a)  Sample the temperature distribution & add to df
  temps_T2 <- rep(rtruncnorm(n_sim,a=0,b=10,mean=8.5,sd=1.0),each=n_units) #we made this up
  data$T_T2 <- temps_T2
  
  ## (b) Sample the transportion time (in days) distribution & add to dataframe
  times_T2 <- rep(rtruncnorm(n_sim,a=0.01,b=0.24, mean=0.04,sd=0.02),each=n_units)
  data$t_T2 <- times_T2
  
  ## (c) Determine Lag_S and Mu_S
  ## determine lag & mumax at new temp
  for (i in 1:(n_sim*n_units)){
    #Find row in growth parameter data that corresponds to allele sample
    allele_index <- which(spore_growth_import$STorAT == data$STorAT[i])
    
    #Calculate the new growth parameters at new temp using the square root model 
    newLag_T2 <- lagAtNewTemp(data$T_T2[i], spore_growth_import$lag[allele_index])
    newMu_T2 <-  muAtNewTemp(data$T_T2[i], spore_growth_import$mumax[allele_index])
    
    data$newLag_T2[i] <- newLag_T2
    data$newMu_T2[i] <- newMu_T2
  }
  
  ## Determine Lag_T2
  data$checkLagPhase_T2 <- ifelse(((data$t_F/data$newLag_F)+(data$t_T/data$Lag_T)+(data$t_S/data$Lag_S)) <1, 1, 0)
  data$adjLag_T2 <- (1 - (data$t_F/data$newLag_F)-(data$t_T/data$Lag_T)-(data$t_S/data$Lag_S))*data$newLag_T2 
  data$Lag_T2 <- ifelse(data$checkLagPhase_T2==0,0,data$adjLag_T2)
  
  ## Determine Mu_T2
  data$Mu_T2 <- ifelse(data$T_T2 >= data$T_S*0.75 & data$T_T2 <= data$T_S*1.25, data$newMu_S, data$newMu_T2)
  data$checkTemp_StoT2 <- ifelse(data$T_T2 >= data$T_S*0.75 & data$T_T2 <= data$T_S*1.25, "T_S", "T_T2") 
  
  ## (d) Determine count_T2
  for (i in 1:(n_sim *n_units)){
    #Find row in growth parameter data that corresponds to allele sample
    allele_index <- which(spore_growth_import$STorAT == data$STorAT[i]) 
    
    #Calculate the log10N count using the new growth parameters
    data$count_T2[i] <- log10N_func(data$t_T2[i], data$Lag_T2[i],data$Mu_T2[i],data$count_S[i],spore_growth_import$LOG10Nmax[allele_index])
  }

#-----------------------------------------------------Stage 5: Storage at home----------------------------------------------------------
  ## (a)  Sample the temperature distribution & add to df
  temps <- rep(NA, n_sim)
  
  for (i in 1:(n_sim*n_units)){
    number <- rlaplace(1,m=4.06,s=2.31)
    while (number > 15 | number < -1) {
      number <- rlaplace(1,m=4.06,s=2.31) #make sure that this cannot be >15 or < -1
    }
    temps[i] <- number
  }
  
  data$T_H <- temps 
  
  ## (b) Define t_H as 1 d for all units; this is just the first day of home storage
  data$t_H <- rep(1, each = n_sim*n_units)
  
  ## (c) Determine Lag_H and Mu_H
  ## determine lag & mumax at new temp
  for (i in 1:(n_sim*n_units)){
    #Find row in growth parameter data that corresponds to allele sample
    allele_index <- which(spore_growth_import$STorAT == data$STorAT[i])
    
    #Calculate the new growth parameters at new temp using the square root model 
    newLag_H <- lagAtNewTemp(data$T_H[i], spore_growth_import$lag[allele_index])
    newMu_H <-  muAtNewTemp(data$T_H[i], spore_growth_import$mumax[allele_index])
    
    data$newLag_H[i] <- newLag_H
    data$newMu_H[i] <- newMu_H
  }
  
  ## Determine Lag_H
  data$checkLagPhase_H <- ifelse(((data$t_F/data$newLag_F)+(data$t_T/data$Lag_T)+(data$t_S/data$Lag_S)+(data$t_T2/data$Lag_T2)) <1, 1, 0)
  data$adjLag_H <- (1 - (data$t_F/data$newLag_F)-(data$t_T/data$Lag_T)-(data$t_S/data$Lag_S)-(data$t_T2/data$Lag_T2))*data$newLag_H 
  data$Lag_H <- ifelse(data$checkLagPhase_H==0,0,data$adjLag_H)
  
  ## Determine Mu_H
  data$Mu_H <- ifelse(data$T_H >= data$T_T2*0.75 & data$T_H <= data$T_T2*1.25, data$newMu_T2, data$newMu_H)
  data$checkTemp_T2toH <- ifelse(data$T_H >= data$T_T2*0.75 & data$T_H <= data$T_T2*1.25, "T_T2", "T_H") 
  
  ## (d) Determine count_H (this is the count for day 1 of home storage)
  for (i in 1:(n_sim *n_units)){
    #Find row in growth parameter data that corresponds to allele sample
    allele_index <- which(spore_growth_import$STorAT == data$STorAT[i]) 
    
    #Calculate the log10N count using the new growth parameters
    data$count_H[i] <- log10N_func(data$t_H[i], data$Lag_H[i],data$Mu_H[i],data$count_T2[i],spore_growth_import$LOG10Nmax[allele_index])
  }
  
  ## (e) simulate the rest of storage 
  df2 <- data[c(1:3,17,18,48,49,50,20)]
  start_day <- 2 
  end_day <- input$day

  actual_t_H <- rep(rep(seq(start_day, end_day)), times = n_units*n_sim) 
  sim_day <- rep(rep(seq(start_day-1, end_day-1)), times = n_units*n_sim)

  df3 <- df2[rep(seq_len(nrow(df2)), each = end_day-1), ] 
  
  # add actual_t_H (i.e., actual day of home storage) and sim_day to dataframe
  df3$actual_t_H <- actual_t_H
  df3$sim_day <- sim_day
  
  # make an empty vector for "count" with the length of n_sim * n_units * n_day
  count <- vector(mode = "logical", n_sim * n_units * (end_day-1))
  
  # add empty count to dataframe
  df3$count <- count
  
  ## (d) Determine count for each day of shelf-life for each unit
  for (i in 1:(n_sim*n_units*(end_day-1))){
    #Find row in growth parameter data that corresponds to allele sample
    allele_index <- which(spore_growth_import$STorAT == df3$STorAT[i]) 
    
    #Calculate the log10N count using the new growth parameters
    df3$count[i] <- log10N_func(df3$sim_day[i], df3$Lag_H[i],df3$Mu_H[i],df3$count_H[i],spore_growth_import$LOG10Nmax[allele_index])
  }
  
  per_spoiled = vector()
  for (i in 1:(end_day-1)){
  per_spoiled = c(per_spoiled, sum(df3$count[sim_day==i] > input$threshold)/length(df3$count[sim_day==i]))
  }
  per_spoiled = per_spoiled * 100
  days = c(2:end_day)
  result = cbind(days, per_spoiled) %>% as.data.frame()
  result
})  

#--------------------------------------------------------Display results----------------------------------------------------------
output$plot <- renderPlot({
    
  ggplot(data = model_result(), aes(x = days, y = per_spoiled))+
      geom_line(aes(y=per_spoiled))+
      labs(title="Simulated % of spoiled half-gallon milk containers due to outgrowth of psychrotolerant sporeformers",
           x="Consumer storage (days)",
           y="% of spoiled half-gallon milk containers")+
      theme_classic()+
      theme(plot.title = element_text(size=10))+
      ylim(0,100)+
      scale_x_continuous(breaks = scales::pretty_breaks(10))
  },
  res = 108)

output$rawdata <- renderPrint({
  model_result()
})

output$downloadData <- downloadHandler(
  
  filename = function() { "spoilage_data.csv" },
  
  content = function(file) {
    write.csv(model_result(), file, row.names = FALSE)
  })





})


