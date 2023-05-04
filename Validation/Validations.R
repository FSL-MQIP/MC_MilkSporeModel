## Check sporeformer distribution of NY and TX milk
TX_df = read.csv("InputFiles/Texas_PSC.csv")
TX_PSC = TX_df %>% filter(test == "SP PSC") %>% pull(total) %>% na.omit()



######################################################################################################
# validation against TX and NY data
TX = read.csv("InputFiles/ValidationData2.csv", header = TRUE)
TX_AT = read.csv("InputFiles/Texas_ATFreq.csv")
TX_ATfreq = TX_AT$AT
TX_ATfreq = ifelse(TX_ATfreq == 3513, 513, TX_ATfreq) #assigning 3513 to 3 only
TX_ATfreq = TX_ATfreq[!TX_ATfreq %in% c(23, 159)]


########################################################################### TX Dataset ###############################################################
# TX Validation data cleaning
nsim = 10000
TX = TX[TX$CVTA=="N" &  TX$Day != 17 & TX$Day != 49 & TX$Day != 63,]
TX.Val = tibble(rep(TX$Day, nsim), rep(TX$Storage_Temp,nsim), rep(TX$APC,nsim))
names(TX.Val) = c("Day","Storage_Temp","APC")

TX.Val$AT = sample(TX_ATfreq, nrow(TX.Val), replace = T)
GrowthParas$AT = c(3,100,15,513,61,45,179,340,17,23,159)
index = sapply(TX.Val$AT, function(x) which(x==GrowthParas$rpoBAT))

TX.Val$Log10_DayInt = rnorm(nrow(TX.Val), -0.72, 0.99)

TX.Val$Mu = muAtNewTemp(TX.Val$Storage_Temp, GrowthParas$mumax[index])
TX.Val$Lag = lagAtNewTemp(TX.Val$Storage_Temp, GrowthParas$lag[index])
TX.Val$Nmax = GrowthParas$LOG10Nmax[index]

TX.Microflora = TX.Val %>% filter(Day == 0) %>% pull(APC) %>% mean()

TX.Val$Pred = log10N_func(TX.Val$Day, TX.Val$Lag, TX.Val$Mu, TX.Val$Log10_DayInt, TX.Val$Nmax)
TX.Val$PredAdj = log10(10^TX.Val$Pred + TX.Microflora)
TX.Val$Actual = ifelse(TX.Val$APC == 0, log10(2.5), log10(TX.Val$APC))

TX.Val$Day = as.factor(TX.Val$Day)
TX.Val_tidy = pivot_longer(TX.Val, c(PredAdj, Actual), names_to = "Type", values_to = "Count")
TX.Val_T3 = TX.Val_tidy %>% filter(Storage_Temp == 3) %>% as.data.frame()
TX.Val_T6.5 = TX.Val_tidy %>% filter(Storage_Temp == 6.5) %>% as.data.frame()
TX.Val_T10 = TX.Val_tidy %>% filter(Storage_Temp == 10) %>% as.data.frame()


## RMSE
Val2_Actual_3C =  TX.Val_T3 %>% group_by(Day) %>% filter(Type == "Actual") %>% summarise(mean_count = mean(Count))
Val2_Pred_3C =  TX.Val_T3 %>% group_by(Day) %>% filter(Type == "PredAdj") %>% summarise(mean_count = mean(Count))
Val2_RMSE_3C = sqrt(sum((Val2_Actual_3C$mean_count - Val2_Pred_3C$mean_count)^2)/length(Val2_Pred_3C))
print(paste("RMSE for Val2 at 3C is ", Val2_RMSE_3C))

Val2_Actual_6.5C =  TX.Val_T6.5 %>% group_by(Day) %>% filter(Type == "Actual") %>% summarise(mean_count = mean(Count))
Val2_Pred_6.5C =  TX.Val_T6.5 %>% group_by(Day) %>% filter(Type == "PredAdj") %>% summarise(mean_count = mean(Count))
Val2_RMSE_6.5C = sqrt(sum((Val2_Actual_6.5C$mean_count - Val2_Pred_6.5C$mean_count)^2)/length(Val2_Pred_6.5C))
print(paste("RMSE for Val2 at 6.5C is ", Val2_RMSE_6.5C))

Val2_Actual_10C =  TX.Val_T10 %>% group_by(Day) %>% filter(Type == "Actual") %>% summarise(mean_count = mean(Count))
Val2_Pred_10C =  TX.Val_T10 %>% group_by(Day) %>% filter(Type == "PredAdj") %>% summarise(mean_count = mean(Count))
Val2_RMSE_10C = sqrt(sum((Val2_Actual_10C$mean_count - Val2_Pred_10C$mean_count)^2)/length(Val2_Pred_10C))
print(paste("RMSE for Val2 at 10C is ", Val2_RMSE_10C))

Val2_RMSE = sqrt(sum((Val2_Actual_3C$mean_count - Val2_Pred_3C$mean_count)^2, 
                (Val2_Actual_6.5C$mean_count - Val2_Pred_6.5C$mean_count)^2,
                (Val2_Actual_10C$mean_count - Val2_Pred_10C$mean_count)^2)/sum(
                  length(Val2_Pred_3C), length(Val2_Pred_6.5C), length(Val2_Pred_10C)))
print(paste("RMSE for Val2 is ", Val2_RMSE))


Val2_Actual_3C =  TX.Val_T3 %>% filter(as.numeric(as.character(Day)) <= 21) %>% group_by(Day) %>% filter(Type == "Actual") %>% summarise(mean_count = mean(Count))
Val2_Pred_3C =  TX.Val_T3 %>% filter(as.numeric(as.character(Day)) <= 21) %>% group_by(Day) %>% filter(Type == "PredAdj") %>% summarise(mean_count = mean(Count))
Val2_RMSE_3C = sqrt(sum((Val2_Actual_3C$mean_count - Val2_Pred_3C$mean_count)^2)/length(Val2_Pred_3C))
print(paste("RMSE for Val2 at 3C prior to 21 days is ", Val2_RMSE_3C))

Val2_Actual_6.5C =  TX.Val_T6.5 %>% filter(as.numeric(as.character(Day)) <= 21) %>% group_by(Day) %>% filter(Type == "Actual") %>% summarise(mean_count = mean(Count))
Val2_Pred_6.5C =  TX.Val_T6.5 %>% filter(as.numeric(as.character(Day)) <= 21) %>% group_by(Day) %>% filter(Type == "PredAdj") %>% summarise(mean_count = mean(Count))
Val2_RMSE_6.5C = sqrt(sum((Val2_Actual_6.5C$mean_count - Val2_Pred_6.5C$mean_count)^2)/length(Val2_Pred_6.5C))
print(paste("RMSE for Val2 at 6.5C prior to 21 days is ", Val2_RMSE_6.5C))

Val2_Actual_10C =  TX.Val_T10 %>% filter(as.numeric(as.character(Day)) <= 21) %>% group_by(Day) %>% filter(Type == "Actual") %>% summarise(mean_count = mean(Count))
Val2_Pred_10C =  TX.Val_T10 %>% filter(as.numeric(as.character(Day)) <= 21) %>% group_by(Day) %>% filter(Type == "PredAdj") %>% summarise(mean_count = mean(Count))
Val2_RMSE_10C = sqrt(sum((Val2_Actual_10C$mean_count - Val2_Pred_10C$mean_count)^2)/length(Val2_Pred_10C))
print(paste("RMSE for Val2 at 10C prior to 21 days is ", Val2_RMSE_10C))

Val2_RMSE = sqrt(sum((Val2_Actual_3C$mean_count - Val2_Pred_3C$mean_count)^2, 
                     (Val2_Actual_6.5C$mean_count - Val2_Pred_6.5C$mean_count)^2,
                     (Val2_Actual_10C$mean_count - Val2_Pred_10C$mean_count)^2)/sum(
                       length(Val2_Pred_3C), length(Val2_Pred_6.5C), length(Val2_Pred_10C)))
print(paste("RMSE for Val2 prior to 21 days is ", Val2_RMSE))

#Visualization
ggplot(data = TX.Val_T3, aes(x=Day, y=Count, fill=Type))+
  geom_boxplot( position = position_dodge(0.9),
                show.legend = FALSE)+
  #geom_jitter(width = 0.1, aes(color = Type), alpha = 0.2) +
  scale_fill_manual(values = c("white", "grey"))+
  #facet_wrap(~Type)+
  theme_classic()+
  #scale_x_discrete(breaks = sort(unique(Val2_T3$Day)))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (3°C) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("Figures/Fig3.Val2_TX_3C.tiff", height = 3, width = 8, units = "in", dpi = "retina")

ggplot(data = TX.Val_T6.5, aes(x=Day, y=Count, fill=Type))+
  geom_boxplot( position = position_dodge(0.9),
                show.legend = FALSE)+
  #geom_jitter(width = 0.1, aes(color = Type), alpha = 0.2) +
  scale_fill_manual(values = c("white", "grey"))+
  #facet_wrap(~Type)+
  theme_classic()+
  #scale_x_discrete(breaks = sort(unique(Val2_T3$Day)))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (6.5°C) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("Figures/Fig3.Val2_TX_6.5C.tiff", height = 3, width = 8, units = "in", dpi = "retina")

ggplot(data = TX.Val_T10, aes(x=Day, y=Count, fill=Type))+
  geom_boxplot( position = position_dodge(0.9),
                show.legend = FALSE)+
  #geom_jitter(width = 0.1, aes(color = Type), alpha = 0.2) +
  scale_fill_manual(values = c("white", "grey"))+
  #facet_wrap(~Type)+
  theme_classic()+
  #scale_x_discrete(breaks = sort(unique(Val2_T3$Day)))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (10°C) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("Figures/Fig3.Val2_TX_10C.tiff", height = 3, width = 8, units = "in", dpi = "retina")

############################################################## NY Dataset ######################################################
NY = read.csv("InputFiles/ValidationData3.csv", header = TRUE)
NY.Val = NY[NY$PPC == "N" & NY$Storage_Temp != "SC" & NY$Day != 42,]
NY.Val$Storage_Temp = as.numeric(NY.Val$Storage_Temp)
NY.Val$APC = as.numeric(NY.Val$APC)
NY.Val = data.frame(rep(NY.Val$Source.ID, 100), rep(NY.Val$Day, 100), rep(NY.Val$Storage_Temp,100), rep(NY.Val$APC,100))
names(NY.Val) = c("Source.ID","Day","Storage_Temp","APC")

NY.AT <- vector()
for (i in 1:(nrow(NY.Val))){
  AT_samp <- sample(ATfreq, 1,replace = T)
  while(AT_samp == "AT_23" || AT_samp == "AT_159"){
    AT_samp <- sample(ATfreq, 1,replace = T)}
  NY.AT[i] = AT_samp
}

NY.Val$AT = NY.AT
index = sapply(NY.Val$AT, function(x) which(x==GrowthParas$rpoBAT))

NY.Val$Log10_DayInt = rnorm(nrow(NY.Val), -0.72, 0.99)
NY.Val$Log10_DayInt = truncnorm::rtruncnorm(nrow(NY.Val), a = log10(1/200), b = Inf, mean = -0.72, sd = 0.99)

NY.Val$Mu = muAtNewTemp(NY.Val$Storage_Temp, GrowthParas$mumax[index])
NY.Val$Lag = lagAtNewTemp(NY.Val$Storage_Temp, GrowthParas$lag[index])
NY.Val$Nmax = GrowthParas$LOG10Nmax[index]

NY.Microflora.SU1 = NY.Val %>% filter(Day == 0 & Source.ID == "SU1") %>% pull(APC) %>% mean()
NY.Microflora.SU2 = NY.Val %>% filter(Day == 0 & Source.ID == "SU2") %>% pull(APC) %>% mean()
#NY.Val$Microflora = sample(NY.Microflora, nrow(NY.Val), replace = T)

NY.Val$Pred = log10N_func(NY.Val$Day, NY.Val$Lag, NY.Val$Mu, NY.Val$Log10_DayInt, NY.Val$Nmax)


NY.Val$PredAdj = ifelse(NY.Val$Source.ID == "SU1", 
                        log10(10^NY.Val$Pred + NY.Microflora.SU1),
                        log10(10^NY.Val$Pred + NY.Microflora.SU2))
NY.Val$Actual = ifelse(NY.Val$APC == 0, log10(2.5), log10(NY.Val$APC))


NY.Val_tidy = pivot_longer(NY.Val, c(PredAdj, Actual), names_to = "Type", values_to = "Count")
NY.Val_tidy$Day = as.factor(NY.Val_tidy$Day)
NY.Val_T4 = NY.Val_tidy %>% filter(Storage_Temp == 4) %>% as.data.frame()
NY.Val_T6 = NY.Val_tidy %>% filter(Storage_Temp == 6) %>% as.data.frame()

#Visualization
ggplot(data = NY.Val_T4, aes(x=Day, y=Count, fill=Type))+
  geom_boxplot( position = position_dodge(0.9),
                show.legend = FALSE)+
  #geom_jitter(width = 0.1, aes(color = Type), alpha = 0.2) +
  scale_fill_manual(values = c("white", "grey"))+
  #facet_wrap(~Source.ID)+
  theme_classic()+
  #scale_x_discrete(breaks = sort(unique(Val2_T3$Day)))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (4°C) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("Figures/Fig4.Val3_NY_4C.tiff", height = 3, width = 8, units = "in", dpi = "retina")

ggplot(data = NY.Val_T6, aes(x=Day, y=Count, fill=Type))+
  geom_boxplot( position = position_dodge(0.9),
                show.legend = FALSE)+
  #geom_jitter(width = 0.1, aes(color = Type), alpha = 0.2) +
  scale_fill_manual(values = c("white", "grey"))+
  #facet_wrap(~Source.ID)+
  theme_classic()+
  #scale_x_discrete(breaks = sort(unique(Val2_T3$Day)))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (6°C) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("Figures/Fig4.Val3_NY_6C.tiff", height = 3, width = 8, units = "in", dpi = "retina")

### simulated SC
nsim = 10000
NY.Val.SC = NY[NY$PPC == "N" & NY$Storage_Temp == "SC" & NY$Day != 42,]
NY.Val.SC$APC = as.numeric(NY.Val.SC$APC)
NY.Val.SC = data.frame(rep(NY.Val.SC$Source.ID, each = nsim),
                       rep(NY.Val.SC$Day, each = nsim), 
                       rep(NY.Val.SC$Storage_Temp, each = nsim), 
                       rep(NY.Val.SC$APC, each = nsim))
names(NY.Val.SC) = c("Source.ID","Day","Temp","APC")

NY.AT.SC <- vector()
for (i in 1:(nrow(NY.Val.SC))){
  AT_samp <- sample(ATfreq, 1,replace = T)
  while(AT_samp == "AT_23" || AT_samp == "AT_159"){
    AT_samp <- sample(ATfreq, 1,replace = T)}
  NY.AT.SC[i] = AT_samp
}

NY.Val.SC$AT = NY.AT.SC
index = sapply(NY.Val.SC$AT, function(x) which(x==GrowthParas$rpoBAT))

NY.Val.SC$Log10_DayInt = rnorm(nrow(NY.Val.SC), -0.72, 0.99)

NY.Val.SC$Mu_2 = muAtNewTemp(2, GrowthParas$mumax[index])
NY.Val.SC$Lag_2 = lagAtNewTemp(2, GrowthParas$lag[index])
NY.Val.SC$Mu_4 = muAtNewTemp(4, GrowthParas$mumax[index])
NY.Val.SC$Lag_4 = lagAtNewTemp(4, GrowthParas$lag[index])
NY.Val.SC$Mu_10 = muAtNewTemp(10, GrowthParas$mumax[index])
NY.Val.SC$Lag_10 = lagAtNewTemp(10, GrowthParas$lag[index])
NY.Val.SC$Nmax = GrowthParas$LOG10Nmax[index]

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

NY.Val.SC$PredAdj = ifelse(NY.Val.SC$Source.ID == "SU1", 
                           log10(10^NY.Val.SC$Pred + NY.Microflora.SU1),
                           log10(10^NY.Val.SC$Pred + NY.Microflora.SU2))
NY.Val.SC$Actual = ifelse(NY.Val.SC$APC == 0, log10(2.5), log10(NY.Val.SC$APC))
NY.Val.SC$Day = as.factor(NY.Val.SC$Day)
NY.Val.SC_tidy = pivot_longer(NY.Val.SC, c(PredAdj, Actual), names_to = "Type", values_to = "Count")

#Visualization
ggplot(data = NY.Val.SC_tidy, aes(x=Day, y=Count, fill=Type))+
  geom_boxplot( position = position_dodge(0.9),
                show.legend = FALSE)+
  #geom_jitter(width = 0.1, aes(color = Type), alpha = 0.2) +
  scale_fill_manual(values = c("white", "grey"))+
  #facet_wrap(~Source.ID)+
  theme_classic()+
  #scale_x_discrete(breaks = sort(unique(Val2_T3$Day)))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (SC) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("Figures/Fig4.Val3_NY_SC.tiff", height = 3, width = 8, units = "in", dpi = "retina")


## RMSE
Val3_Actual_4C =  NY.Val_T4 %>% group_by(Day) %>% filter(Type == "Actual") %>% summarise(mean_count = mean(Count))
Val3_Pred_4C =  NY.Val_T4 %>% group_by(Day) %>% filter(Type == "PredAdj") %>% summarise(mean_count = mean(Count))
Val3_RMSE_4C = sqrt(sum((Val3_Actual_4C$mean_count - Val3_Pred_4C$mean_count)^2)/length(Val3_Pred_4C))
print(paste("RMSE for Val3 at 4C is ", Val3_RMSE_4C))

Val3_Actual_6C =  NY.Val_T6 %>% group_by(Day) %>% filter(Type == "Actual") %>% summarise(mean_count = mean(Count))
Val3_Pred_6C =  NY.Val_T6 %>% group_by(Day) %>% filter(Type == "PredAdj") %>% summarise(mean_count = mean(Count))
Val3_RMSE_6C = sqrt(sum((Val3_Actual_6C$mean_count - Val3_Pred_6C$mean_count)^2)/length(Val3_Pred_6C))
print(paste("RMSE for Val3 at 6C is ", Val3_RMSE_6C))

Val3_Actual_SC =  NY.Val.SC_tidy %>% group_by(Day) %>% filter(Type == "Actual") %>% summarise(mean_count = mean(Count))
Val3_Pred_SC =  NY.Val.SC_tidy %>% group_by(Day) %>% filter(Type == "PredAdj") %>% summarise(mean_count = mean(Count))
Val3_RMSE_SC = sqrt(sum((Val3_Actual_SC$mean_count - Val3_Pred_SC$mean_count)^2)/length(Val3_Pred_SC))
print(paste("RMSE for Val3 at SC is ", Val3_RMSE_SC))

Val3_RMSE = sqrt(sum((Val3_Actual_4C$mean_count - Val3_Pred_4C$mean_count)^2, 
                     (Val3_Actual_6C$mean_count - Val3_Pred_6C$mean_count)^2,
                     (Val3_Actual_SC$mean_count - Val3_Pred_SC$mean_count)^2)/sum(
                       length(Val3_Pred_4C), length(Val3_Pred_6C), length(Val3_Pred_SC)))
print(paste("RMSE for Val3 is ", Val3_RMSE))

## RMSE less than 21 days
Val3_Actual_4C =  NY.Val_T4 %>% filter(as.numeric(as.character(Day)) <= 21) %>% group_by(Day) %>% filter(Type == "Actual") %>% summarise(mean_count = mean(Count))
Val3_Pred_4C =  NY.Val_T4 %>% filter(as.numeric(as.character(Day)) <= 21) %>% group_by(Day) %>% filter(Type == "PredAdj") %>% summarise(mean_count = mean(Count))
Val3_RMSE_4C = sqrt(sum((Val3_Actual_4C$mean_count - Val3_Pred_4C$mean_count)^2)/length(Val3_Pred_4C))
print(paste("RMSE for Val3 at 4C prior to 21 days is ", Val3_RMSE_4C))

Val3_Actual_6C =  NY.Val_T6 %>% filter(as.numeric(as.character(Day)) <= 21) %>% group_by(Day) %>% filter(Type == "Actual") %>% summarise(mean_count = mean(Count))
Val3_Pred_6C =  NY.Val_T6 %>% filter(as.numeric(as.character(Day)) <= 21) %>% group_by(Day) %>% filter(Type == "PredAdj") %>% summarise(mean_count = mean(Count))
Val3_RMSE_6C = sqrt(sum((Val3_Actual_6C$mean_count - Val3_Pred_6C$mean_count)^2)/length(Val3_Pred_6C))
print(paste("RMSE for Val3 at 6C prior to 21 days is ", Val3_RMSE_6C))

Val3_Actual_SC =  NY.Val.SC_tidy %>% filter(as.numeric(as.character(Day)) <= 21) %>% group_by(Day) %>% filter(Type == "Actual") %>% summarise(mean_count = mean(Count))
Val3_Pred_SC =  NY.Val.SC_tidy %>% filter(as.numeric(as.character(Day)) <= 21) %>% group_by(Day) %>% filter(Type == "PredAdj") %>% summarise(mean_count = mean(Count))
Val3_RMSE_SC = sqrt(sum((Val3_Actual_SC$mean_count - Val3_Pred_SC$mean_count)^2)/length(Val3_Pred_SC))
print(paste("RMSE for Val3 at SC prior to 21 days is ", Val3_RMSE_SC))

Val3_RMSE = sqrt(sum((Val3_Actual_4C$mean_count - Val3_Pred_4C$mean_count)^2, 
                     (Val3_Actual_6C$mean_count - Val3_Pred_6C$mean_count)^2,
                     (Val3_Actual_SC$mean_count - Val3_Pred_SC$mean_count)^2)/sum(
                       length(Val3_Pred_4C), length(Val3_Pred_6C), length(Val3_Pred_SC)))
print(paste("RMSE for Val3 prior to 21 days is ", Val3_RMSE))