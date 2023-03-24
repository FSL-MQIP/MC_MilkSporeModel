####################################### 1. Set up environment ###########################################
# load packages
library(shiny)
library(readr)
library(censReg)
library(fitdistrplus)
library(splitstackshape)
library(rmutil)
library(truncnorm)
library(EnvStats)
library(epiR)
library(stats)
library(tidyverse)

# load functions
source("UtilityFunctions.R")

# set seed
set.seed(1)


######################################## 2. Prepare input files ##########################################
# load input files
## AT frequency data (from Ariel_2017MC_SporeModel on MC-2020 on GitHub)
ATfreq = read.csv("InputFiles/Baseline_ATFreq.csv", stringsAsFactors = FALSE, header = TRUE)

## growth parameter data; make sure this contains growth parameters & growth model name
GrowthParas = read.csv("InputFiles/GrowthParameters_NoGrowthATs.csv" , stringsAsFactors = FALSE)

## initial microbial count data (from Ariel_2017MC_SporeModel on MC-2020 on GitHub)
InitConcData = read.csv("InputFiles/Baseline_InitialSporeCount.csv")


######################################## 3. Organize model parameters ####################################
## (a) Simulation setting
n_sim = 100
n_units = 10

lot_id = rep(seq(1, n_sim), each = n_units)
unit_id = rep(seq(1,n_units), times = n_sim) 
ModelData = tibble(lot_id, unit_id)

## (b) AT frequency 
ATfreq = ATfreq %>% pull(ClosestAT)

GrowthParas = GrowthParas %>%
  select(rpoBAT, lag, mumax, LOG10Nmax)

growthAT = vector()
for (i in 1:(n_units * n_sim)){
  AT_samp = sample(ATfreq, 1,replace = T)
  while(AT_samp == "AT_23" || AT_samp == "AT_159"){
    AT_samp = sample(ATfreq, 1,replace = T)}
  growthAT[i] = AT_samp
}

ModelData$AT = growthAT
ModelData$AT_index = match(ModelData$AT, GrowthParas$rpoBAT)
ModelData$log10Nmax = GrowthParas$LOG10Nmax[ModelData$AT_index]

## (c) Time and temperature distributions 
#stage1, storage at facility
ModelData$T_F = rep(runif(n_sim, min = 3.5, max = 4.5), each = n_units) # temp
ModelData$t_F = rep(runif(n_sim, min = 1, max = 2), each = n_units)     # time

#stage2, transport to retail store
ModelData$T_T = rep(rtri(n_sim, min = 1.7, max = 10.0, mode = 4.4), each = n_units) # temp
ModelData$t_T = rep(rtri(n_sim, min = 1, max = 10, mode = 5), each = n_units)       # time

#stage3, storage/display at retail store
ModelData$T_S = rep(rtruncnorm(n_sim, a = -1.4, b = 5.4, mean = 2.3, sd = 1.8), each = n_units)       # temp
ModelData$t_S = rep(rtruncnorm(n_sim, a = 0.042, b = 10.0, mean = 1.821, sd = 3.3), each = n_units)   # time

#stage4, transport from retail store to homes
ModelData$T_T2 = rep(rtruncnorm(n_sim, a = 0, b = 10, mean = 8.5, sd = 1.0), each = n_units)          # temp
ModelData$t_T2 = rep(rtruncnorm(n_sim, a = 0.01, b = 0.24, mean = 0.04, sd = 0.02), each = n_units)   # time

#stage5, storage at homes
temps = vector()
for (i in 1:(n_sim*n_units)){
  number <- rlaplace(1,m=4.06,s=2.31)
  while (number > 15 | number < -1) {
    number <- rlaplace(1,m=4.06,s=2.31) #make sure that this cannot be >15 or < -1
  }
  temps[i] <- number
}

ModelData$T_H = temps 





######################################### Model #############################################
shinyServer(function(input, output) {
  
#-----------------------------------------------------setting---------------------------------------------
  # get initial count
  
  model_result = reactive({
    
  spore_log10MPN_samp =  rnorm(n_sim, input$count_mean, input$count_sd) 

  ## Convert log10 MPN/mL to MPN/mL
  spore_MPN_samp <- 10^spore_log10MPN_samp
  ## Convert MPN for each 1900 mL milk container
  spore_MPN_samp_halfgal <- spore_MPN_samp * 1900
  #Sample the MPN_init distribution
  spore_MPN_init<-vector()
  for (i in 1:n_sim){
    spore_MPN_init_samp <-rep(rpois(n_units, spore_MPN_samp_halfgal[i]))
    spore_MPN_init<-c(spore_MPN_init, spore_MPN_init_samp)}
  # First convert spore_MPN_init from half-gallon to mLs
  spore_MPN_init[spore_MPN_init < 1] = 0
  spore_MPN_init_mL = spore_MPN_init / 1900
  # Remove 0's from the data and replace with detection limit
  spore_MPN_init_mL[spore_MPN_init_mL <= 0 ] = 0.01
  
  
  ModelData$log10InitConc = log10(spore_MPN_init_mL)
  
#-----------------------------------------------------intervention strategies---------------------------------------------
  # spore reduction 
  if (input$sporeRed == "mf") {ModelData$log10InitConc = ModelData$log10InitConc - 2.2}
  else if (input$sporeRed == "bf1") {ModelData$log10InitConc = ModelData$log10InitConc - 1.4}
  else if (input$sporeRed == "bf2") {ModelData$log10InitConc = ModelData$log10InitConc - 2}
  else {}
  
  # temp control at facility
  if (input$f_intervention == "f_reduceT") {ModelData$T_F = rep(runif(n_sim, min = 2.5, max = 3.5), each = n_units)}
  else if (input$f_intervention == "f_supercool") {ModelData$T_F = rep(runif(n_sim, min = 0.5, max = 1.5), each = n_units)}
  else if (input$f_intervention == "f_reduceVar") {ModelData$T_F = rep(runif(n_sim, min = 3.75, max = 4.25), each = n_units)}
  else {}
    
  # temp control at transportation
  if (input$ftr_intervention == "ftr_alarm") {ModelData$T_T = rep(rtri(n_sim, min = 1.7, max = 6.0, mode = 4.4), each = n_units)}
  else if (input$ftr_intervention == "ftr_opt") {ModelData$t_T = rep(rtri(n_sim, min = 1, max = 7, mode = 5), each = n_units)}
  else {}
  
  # temp control at retail
  if (input$r_intervention == "r_reduceT") {ModelData$T_S = rep(rtruncnorm(n_sim, a = -1.4, b = 5.4, mean = 1.8, sd = 1.8), each = n_units)}
  else if (input$r_intervention == "r_alarm") {ModelData$T_S = rep(rtruncnorm(n_sim, a = -1.4, b = 4, mean = 2.3, sd = 1.8), each = n_units)}
  else if (input$r_intervention == "r_reduceVar") {ModelData$T_S = rep(rtruncnorm(n_sim, a = -1.4, b = 5.4, mean = 2.3, sd = 0.9), each = n_units)}
  else {}
  
#----------------------------------------------- set up growth parameters------------------------------------------------
  # stage1, storage at facility
  ModelData$newLag_F = lagAtNewTemp(ModelData$T_F, GrowthParas$lag[ModelData$AT_index])
  ModelData$newMu_F =  muAtNewTemp(ModelData$T_F, GrowthParas$mumax[ModelData$AT_index])
  
  # stage2, transport to retail store
  ModelData$newLag_T = lagAtNewTemp(ModelData$T_T, GrowthParas$lag[ModelData$AT_index])
  ModelData$newMu_T =  muAtNewTemp(ModelData$T_T, GrowthParas$mumax[ModelData$AT_index])
  
  # stage3, storage/display at retail store
  ModelData$newLag_S = lagAtNewTemp(ModelData$T_S, GrowthParas$lag[ModelData$AT_index])
  ModelData$newMu_S =  muAtNewTemp(ModelData$T_S, GrowthParas$mumax[ModelData$AT_index])
  
  # stage4, transport from retail store to homes
  ModelData$newLag_T2 = lagAtNewTemp(ModelData$T_T2, GrowthParas$lag[ModelData$AT_index])
  ModelData$newMu_T2 =  muAtNewTemp(ModelData$T_T2, GrowthParas$mumax[ModelData$AT_index])
  
  # stage5, storage at homes
  ModelData$newLag_H = lagAtNewTemp(ModelData$T_H, GrowthParas$lag[ModelData$AT_index])
  ModelData$newMu_H =  muAtNewTemp(ModelData$T_H, GrowthParas$mumax[ModelData$AT_index])
  
#-----------------------------------------------Stage 1: Storage at facility----------------------------------------------
  ModelData = ModelData %>%
    rowwise() %>% 
    mutate(predConc_F = log10N_func(t = t_F, lag = newLag_F, mumax = newMu_F,
                                    LOG10N0 = log10InitConc, LOG10Nmax = log10Nmax))
  
  # stage2, transport to retail store
  ModelData = ModelData %>% 
    rowwise() %>% 
    mutate(Lag_T = max(0, 1 - t_F/newLag_F) * newLag_T) %>%  
    mutate(Mu_T = if_else(T_T >= T_F * 0.75 & T_T <= T_F * 1.25, newMu_F, newMu_T)) %>% 
    mutate(predConc_T = log10N_func(t = t_T, lag = Lag_T, mumax = Mu_T,
                                    LOG10N0 = predConc_F, LOG10Nmax = log10Nmax))
  
  # stage3, storage/display at retail store
  ModelData = ModelData %>% 
    rowwise() %>% 
    mutate(Lag_S = max(0, 1 - t_T/Lag_T) * newLag_S) %>%  
    mutate(Mu_S = if_else(T_S >= T_T * 0.75 & T_S <= T_T * 1.25, newMu_T, newMu_S)) %>% 
    mutate(predConc_S = log10N_func(t = t_S, lag = Lag_S, mumax = Mu_S,
                                    LOG10N0 = predConc_T, LOG10Nmax = log10Nmax))
  
  # stage4, transport from retail store to homes
  ModelData = ModelData %>% 
    rowwise() %>% 
    mutate(Lag_T2 = max(0, 1 - t_S/Lag_S) * newLag_T2) %>%  
    mutate(Mu_T2 = if_else(T_T2 >= T_S * 0.75 & T_T2 <= T_S * 1.25, newMu_S, newMu_T2)) %>% 
    mutate(predConc_T2 = log10N_func(t = t_T2, lag = Lag_T2, mumax = Mu_T2,
                                     LOG10N0 = predConc_S, LOG10Nmax = log10Nmax))
  
  # stage5, storage at homes
  ModelData = ModelData %>% 
    slice(rep(1:n(), each = 35)) 
  
  ModelData$t_H = rep(1:35, times = n_sim * n_units)
  
  ModelData = ModelData %>% 
    rowwise() %>% 
    mutate(Lag_H = max(0, 1 - t_T2/Lag_T2) * newLag_H) %>%  
    mutate(Mu_H = if_else(T_H >= T_T2 * 0.75 & T_H <= T_T2 * 1.25, newMu_T2, newMu_H)) %>% 
    mutate(predConc_H = log10N_func(t = t_H, lag = Lag_H, mumax = Mu_H,
                                    LOG10N0 = predConc_T2, LOG10Nmax = log10Nmax))
  
  # Assign the spoilage to each milk container at each day
  ModelData = ModelData %>% 
    rowwise() %>% 
    mutate(spoil = predConc_H > input$threshold)
  
  # Summarise percent of spoiled milk containers per storage day
  spoilData = ModelData %>% 
    group_by(t_H) %>% 
    summarise(perSpoil = 100 * sum(spoil) / (n_sim * n_units)) %>% 
    rename(days = t_H)
  
  spoilData
})  

#--------------------------------------------------------Display results----------------------------------------------------------
output$plot <- renderPlot({
    
  ggplot(data = model_result(), aes(x = days, y = perSpoil))+
      geom_line(aes(y=perSpoil))+
      labs(title="Simulated % of spoiled half-gallon milk containers due to outgrowth of psychrotolerant sporeformers",
           x="Consumer storage (days)",
           y="% of spoiled half-gallon milk containers")+
      theme_classic()+
      theme(plot.title = element_text(size=10))+
      ylim(0,100)+
      scale_x_continuous(breaks = scales::pretty_breaks(10))
  },
  res = 108)

output$shelfLife <- renderText({
  perSpoil = model_result() %>% pull(perSpoil)
  perSpoil_min = min(perSpoil[perSpoil > input$shelfLife_threshold])
  if (perSpoil_min != Inf){
    shelfLife = model_result() %>% filter(perSpoil == perSpoil_min) %>% pull(days)
    paste0("The predicted shelf life is ", shelfLife, " days.")
  } else {
    "The predicted shelf life is beyond 35 days"
  }
  
})

output$downloadData <- downloadHandler(
  
  filename = function() { "spoilage_data.csv" },
  
  content = function(file) {
    write.csv(model_result(), file, row.names = FALSE)
  })





})


