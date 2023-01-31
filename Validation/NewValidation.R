## Check sporeformer distribution of NY and TX milk
TX_df = read.csv("InputFiles/Texas_PSC.csv")
TX_PSC = TX_df %>% filter(test == "SP PSC") %>% pull(total) %>% na.omit()



######################################################################################################
# validation against TX and NY data
TX = read.csv("InputFiles/ValidationData2.csv", header = TRUE)
TX_AT = read.csv("InputFiles/Texas_ATFreq.csv")
TX_ATfreq = TX_AT$AT
TX_ATfreq = ifelse(TX_ATfreq == 3513, 513, TX_ATfreq) #assigning 3513 to 3 only
################################################################# Validating selected ATs #######################3
TX = TX %>% 
  select(where(~ !(all(is.na(.)) | all(. == ""))))

TX.Val.isolate = TX_AT %>% 
  as_tibble() %>% 
  filter(Trial != "SP PSC") %>% 
  mutate(AT = ifelse(AT == 3513, 513, AT)) %>% 
  mutate_at(vars(Trial, HTST_Temp, Storage_Temp, Day), as.numeric) %>% 
  full_join(., y = TX, by = c("Trial", "HTST_Temp", "Storage_Temp", "Day")) %>% 
  na.omit() %>% 
  filter(CVTA == "N") #%>% 
  #dplyr::slice(rep(1:n(), 100))
  
TX.Val.isolate$DayInt_PSC = sample(TX_PSC, nrow(TX.Val.isolate), replace = T)
TX.Val.isolate$DayInt_PSC = ifelse(TX.Val.isolate$DayInt_PSC == 0 , 0.25, TX.Val.isolate$Log10_DayInt)
TX.Val.isolate$Log10_DayInt = log10(TX.Val.isolate$DayInt_PSC )

index = sapply(TX.Val.isolate$AT, function(x) which(x==spore_growth_import$AT))

TX.Val.isolate$Mu = muAtNewTemp(TX.Val.isolate$Storage_Temp, spore_growth_import$mumax[index])
TX.Val.isolate$Lag = lagAtNewTemp(TX.Val.isolate$Storage_Temp, spore_growth_import$lag[index])
TX.Val.isolate$Nmax = spore_growth_import$LOG10Nmax[index]

TX.Microflora = 62.73

TX.Val.isolate$Pred = log10N_func(TX.Val.isolate$Day, TX.Val.isolate$Lag, TX.Val.isolate$Mu, 
                                  TX.Val.isolate$Log10_DayInt, TX.Val.isolate$Nmax)
TX.Val.isolate$PredAdj = log10(10^TX.Val.isolate$Pred + TX.Microflora)
TX.Val.isolate$Actual = log10(TX.Val.isolate$APC + 1)

TX.Val.iso_tidy = pivot_longer(TX.Val.isolate, c(PredAdj, APCavg), names_to = "Type", values_to = "Count") %>%
  mutate(Day = as.factor(Day))
  
TX.Val.iso_T3 = TX.Val.iso_tidy %>% filter(Storage_Temp == 3) %>% as.data.frame()
TX.Val.iso_T6.5 = TX.Val.iso_tidy %>% filter(Storage_Temp == 6.5) %>% as.data.frame()
TX.Val.iso_T10 = TX.Val.iso_tidy %>% filter(Storage_Temp == 10) %>% as.data.frame()

#Visualization
ggplot(data = TX.Val.iso_T3, aes(x=as.factor(Day), y=Count, fill=Type))+
  geom_boxplot( position = position_dodge(0.9),
                show.legend = TRUE)+
  geom_jitter(width = 0.25, aes(color = Type)) +
  scale_fill_manual(values = c("white", "grey"))+
  #facet_wrap(~Type)+
  theme_classic()+
  #scale_x_discrete(breaks = sort(unique(Val2_T3$Day)))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Predicted concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (3C) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))

ggplot(data = TX.Val.iso_T6.5, aes(x=as.factor(Day), y=Count, fill=Type))+
  geom_boxplot( position = position_dodge(0.9),
                show.legend = TRUE)+
  geom_jitter(width = 0.25, aes(color = Type)) +
  scale_fill_manual(values = c("white", "grey"))+
  #facet_wrap(~Type)+
  theme_classic()+
  #scale_x_discrete(breaks = sort(unique(Val2_T3$Day)))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Predicted concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (6.5C) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))

ggplot(data = TX.Val.iso_T10, aes(x=as.factor(Day), y=Count, fill=Type))+
  geom_boxplot( position = position_dodge(0.9),
                show.legend = TRUE)+
  geom_jitter(width = 0.25, aes(color = Type)) +
  scale_fill_manual(values = c("white", "grey"))+
  #facet_wrap(~Type)+
  theme_classic()+
  #scale_x_discrete(breaks = sort(unique(Val2_T3$Day)))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Predicted concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (10C) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))

########################################################################### TX Dataset ###############################################################
# TX Validation data cleaning
nsim = 50
TX = TX[TX$CVTA=="N" &  TX$Day != 17 & TX$Day != 49 & TX$Day != 63,]
TX.Val = tibble(rep(TX$Day, nsim), rep(TX$Storage_Temp,nsim), rep(TX$APC,nsim))
names(TX.Val) = c("Day","Storage_Temp","APC")

TX.Val$AT = sample(TX_ATfreq, nrow(TX.Val), replace = T)
spore_growth_import$AT = c(3,100,15,513,61,45,179,340,17,23,159)
index = sapply(TX.Val$AT, function(x) which(x==spore_growth_import$AT))

TX.Val$Log10_DayInt = rnorm(nrow(TX.Val), -0.72, 0.99)

TX.Val$Mu = muAtNewTemp(TX.Val$Storage_Temp, spore_growth_import$mumax[index])
TX.Val$Lag = lagAtNewTemp(TX.Val$Storage_Temp, spore_growth_import$lag[index])
TX.Val$Nmax = spore_growth_import$LOG10Nmax[index]

TX.Microflora = TX.Val %>% filter(Day == 0) %>% pull(APC) %>% mean()

TX.Val$Pred = log10N_func(TX.Val$Day, TX.Val$Lag, TX.Val$Mu, TX.Val$Log10_DayInt, TX.Val$Nmax)
TX.Val$PredAdj = log10(10^TX.Val$Pred + TX.Microflora)
TX.Val$Actual = log10(TX.Val$APC + 1)

TX.Val$Day = as.factor(TX.Val$Day)
TX.Val_tidy = pivot_longer(TX.Val, c(PredAdj, Actual), names_to = "Type", values_to = "Count")
TX.Val_T3 = TX.Val_tidy %>% filter(Storage_Temp == 3) %>% as.data.frame()
TX.Val_T6.5 = TX.Val_tidy %>% filter(Storage_Temp == 6.5) %>% as.data.frame()
TX.Val_T10 = TX.Val_tidy %>% filter(Storage_Temp == 10) %>% as.data.frame()

#Visualization
ggplot(data = TX.Val_T3, aes(x=Day, y=Count, fill=Type))+
  geom_boxplot( position = position_dodge(0.9),
                show.legend = TRUE)+
  geom_jitter(width = 0.2, aes(color = Type), alpha = 0.3) +
  scale_fill_manual(values = c("white", "grey"))+
  #facet_wrap(~Type)+
  theme_classic()+
  #scale_x_discrete(breaks = sort(unique(Val2_T3$Day)))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Predicted concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (3C) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))

ggplot(data = TX.Val_T6.5, aes(x=Day, y=Count, fill=Type))+
  geom_boxplot( position = position_dodge(0.9),
                show.legend = TRUE)+
  geom_jitter(width = 0.2, aes(color = Type), alpha = 0.3) +
  scale_fill_manual(values = c("white", "grey"))+
  #facet_wrap(~Type)+
  theme_classic()+
  #scale_x_discrete(breaks = sort(unique(Val2_T3$Day)))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Predicted concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (6.5C) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))

ggplot(data = TX.Val_T10, aes(x=Day, y=Count, fill=Type))+
  geom_boxplot( position = position_dodge(0.9),
                show.legend = TRUE)+
  geom_jitter(width = 0.2, aes(color = Type), alpha = 0.3) +
  scale_fill_manual(values = c("white", "grey"))+
  #facet_wrap(~Type)+
  theme_classic()+
  #scale_x_discrete(breaks = sort(unique(Val2_T3$Day)))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Predicted concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (10C) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))

############################################################## NY Dataset ######################################################
NY = read.csv("InputFiles/ValidationData3.csv", header = TRUE)
NY.Val = NY[NY$PPC == "N" & NY$Storage_Temp != "SC",]
NY.Val$Storage_Temp = as.numeric(NY.Val$Storage_Temp)
NY.Val$APC = as.numeric(NY.Val$APC)
NY.Val = data.frame(rep(NY.Val$Source.ID, 100), rep(NY.Val$Day, 100), rep(NY.Val$Storage_Temp,100), rep(NY.Val$APC,100))
names(NY.Val) = c("Source.ID","Day","Storage_Temp","APC")

NY.AT <- vector()
for (i in 1:(nrow(NY.Val))){
  AT_samp <- sample(AT_freq, 1,replace = T)
  while(AT_samp == "AT_23" || AT_samp == "AT_159"){
    AT_samp <- sample(AT_freq, 1,replace = T)}
  NY.AT[i] = AT_samp
}

NY.Val$AT = NY.AT
index = sapply(NY.Val$AT, function(x) which(x==spore_growth_import$STorAT))

NY.Val$Log10_DayInt = rnorm(nrow(NY.Val), -0.72, 0.99)
NY.Val$Log10_DayInt = truncnorm::rtruncnorm(nrow(NY.Val), a = log10(1/200), b = Inf, mean = -0.72, sd = 0.99)

NY.Val$Mu = muAtNewTemp(NY.Val$Storage_Temp, spore_growth_import$mumax[index])
NY.Val$Lag = lagAtNewTemp(NY.Val$Storage_Temp, spore_growth_import$lag[index])
NY.Val$Nmax = spore_growth_import$LOG10Nmax[index]

NY.Microflora = NY.Val %>% filter(Day == 0) %>% pull(APC) %>% mean()
#NY.Val$Microflora = sample(NY.Microflora, nrow(NY.Val), replace = T)

NY.Val$Pred = log10N_func(NY.Val$Day, NY.Val$Lag, NY.Val$Mu, NY.Val$Log10_DayInt, NY.Val$Nmax)

NY.Val$PredAdj = log10(10^NY.Val$Pred + NY.Microflora)
NY.Val$Actual = log10(NY.Val$APC + 1)


NY.Val_tidy = pivot_longer(NY.Val, c(PredAdj, Actual), names_to = "Type", values_to = "Count")
NY.Val_tidy$Day = as.factor(NY.Val_tidy$Day)
NY.Val_T4 = NY.Val_tidy %>% filter(Storage_Temp == 4) %>% as.data.frame()
NY.Val_T6 = NY.Val_tidy %>% filter(Storage_Temp == 6) %>% as.data.frame()

#Visualization
ggplot(data = NY.Val_T4, aes(x=Day, y=Count, fill=Type))+
  geom_boxplot( position = position_dodge(0.9),
                show.legend = TRUE)+
  geom_jitter(width = 0.2, aes(color = Type), alpha = 0.3) +
  scale_fill_manual(values = c("white", "grey"))+
  facet_wrap(~Source.ID)+
  theme_classic()+
  #scale_x_discrete(breaks = sort(unique(Val2_T3$Day)))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Predicted concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (4C) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))

ggplot(data = NY.Val_T6, aes(x=Day, y=Count, fill=Type))+
  geom_boxplot( position = position_dodge(0.9),
                show.legend = TRUE)+
  geom_jitter(width = 0.2, aes(color = Type), alpha = 0.3) +
  scale_fill_manual(values = c("white", "grey"))+
  facet_wrap(~Source.ID)+
  theme_classic()+
  #scale_x_discrete(breaks = sort(unique(Val2_T3$Day)))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Predicted concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (6C) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))


### simulated SC
nsim = 100
NY.Val.SC = NY[NY$PPC == "N" & NY$Storage_Temp == "SC",]
NY.Val.SC$APC = as.numeric(NY.Val.SC$APC)
NY.Val.SC = data.frame(rep(NY.Val.SC$Source.ID, each = nsim),
                       rep(NY.Val.SC$Day, each = nsim), 
                       rep(NY.Val.SC$Storage_Temp, each = nsim), 
                       rep(NY.Val.SC$APC, each = nsim))
names(NY.Val.SC) = c("Source.ID","Day","Temp","APC")

NY.AT.SC <- vector()
for (i in 1:(nrow(NY.Val.SC))){
  AT_samp <- sample(AT_freq, 1,replace = T)
  while(AT_samp == "AT_23" || AT_samp == "AT_159"){
    AT_samp <- sample(AT_freq, 1,replace = T)}
  NY.AT.SC[i] = AT_samp
}

NY.Val.SC$AT = NY.AT.SC
index = sapply(NY.Val.SC$AT, function(x) which(x==spore_growth_import$STorAT))

NY.Val.SC$Log10_DayInt = rnorm(nrow(NY.Val.SC), -0.72, 0.99)

NY.Val.SC$Mu_2 = muAtNewTemp(2, spore_growth_import$mumax[index])
NY.Val.SC$Lag_2 = lagAtNewTemp(2, spore_growth_import$lag[index])
NY.Val.SC$Mu_4 = muAtNewTemp(4, spore_growth_import$mumax[index])
NY.Val.SC$Lag_4 = lagAtNewTemp(4, spore_growth_import$lag[index])
NY.Val.SC$Mu_10 = muAtNewTemp(10, spore_growth_import$mumax[index])
NY.Val.SC$Lag_10 = lagAtNewTemp(10, spore_growth_import$lag[index])
NY.Val.SC$Nmax = spore_growth_import$LOG10Nmax[index]

NY.SC.Microflora = NY.Val.SC %>% filter(Day == 0) %>% pull(APC) %>% mean()

# 41 hours at 4C
NY.Val.SC$Stage1 = log10N_func(41/24, NY.Val.SC$Lag_4, NY.Val.SC$Mu_4, NY.Val.SC$Log10_DayInt, NY.Val.SC$Nmax)
NY.Val.SC$Lag_remain = ifelse(NY.Val.SC$Lag_4 >= 41/24, NY.Val.SC$Lag_4 - 41/24, 0 )
NY.Val.SC$Lag_prop = NY.Val.SC$Lag_remain/NY.Val.SC$Lag_4
  
# 24 hours at 2C
NY.Val.SC$Stage2 = log10N_func(1, NY.Val.SC$Lag_2 * NY.Val.SC$Lag_prop, NY.Val.SC$Mu_2, NY.Val.SC$Stage1, NY.Val.SC$Nmax)
NY.Val.SC$Lag_remain = ifelse(NY.Val.SC$Lag_2 >= 1, NY.Val.SC$Lag_4 - 1, 0 )
NY.Val.SC$Lag_prop = NY.Val.SC$Lag_remain/NY.Val.SC$Lag_2

# 26 mins at 10C
NY.Val.SC$Stage3 = log10N_func(26/60/24, NY.Val.SC$Lag_10 * NY.Val.SC$Lag_prop, NY.Val.SC$Mu_10, NY.Val.SC$Stage2, NY.Val.SC$Nmax)
NY.Val.SC$Lag_remain = ifelse(NY.Val.SC$Lag_10 >= 26/60/24, NY.Val.SC$Lag_10 - 26/60/24, 0 )
NY.Val.SC$Lag_prop = NY.Val.SC$Lag_remain/NY.Val.SC$Lag_10

# remainder time at 4C
NY.Val.SC$Pred = log10N_func(NY.Val.SC$Day-41/24-1-26/60/24, NY.Val.SC$Lag_4 * NY.Val.SC$Lag_prop, 
                               NY.Val.SC$Mu_4, NY.Val.SC$Stage3, NY.Val.SC$Nmax)

NY.Val.SC$PredAdj = log10(10^NY.Val.SC$Pred + NY.SC.Microflora)
NY.Val.SC$Actual = log10(NY.Val.SC$APC + 1)
NY.Val.SC$Day = as.factor(NY.Val.SC$Day)
NY.Val.SC_tidy = pivot_longer(NY.Val.SC, c(PredAdj, Actual), names_to = "Type", values_to = "Count")

#Visualization
ggplot(data = NY.Val.SC_tidy, aes(x=Day, y=Count, fill=Type))+
  geom_boxplot( position = position_dodge(0.9),
                show.legend = TRUE)+
  geom_jitter(width = 0.2, aes(color = Type), alpha = 0.3) +
  scale_fill_manual(values = c("white", "grey"))+
  facet_wrap(~Source.ID)+
  theme_classic()+
  #scale_x_discrete(breaks = sort(unique(Val2_T3$Day)))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Predicted concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (SC) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))
