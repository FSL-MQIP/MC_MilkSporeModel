#Second validation
#LS AT frequency and validation data

#AT frequency table
LS_AT = read.csv("InputFiles/Texas_ATFreq.csv")
LS_ATfreq = LS_AT$AT
LS_ATfreq = ifelse(LS_ATfreq == 3513, 513, LS_ATfreq) #assigning 3513 to 3 only

#Summary stats
LS_AT_Table = data.frame(LS_AT$Genus, LS_AT$species, LS_AT$AT)
names(LS_AT_Table) = c("Genus","Species","AT")
LS_AT_Table$AT = ifelse(LS_AT_Table$AT == 3513, 513, LS_AT_Table$AT)
LS_AT_Table = LS_AT_Table %>% mutate_all(as.factor)

#Validation dataset
Val2_data = read.csv("InputFiles/ValidationData2.csv")

#Summary stats
Val2_sum = Val2_data %>%  mutate_at(c("CVTA", "Storage_Temp"),as.factor)
Val2_sum %>% count(CVTA == "Y")
Val2_sum %>% filter(Storage_Temp == 10) %>% group_by(Day) %>% summarise(n = n())
Val2_sum %>% filter(Trial == 1) %>% pull(DayInt_SPC) %>% mean()
Val2_sum %>% filter(Trial == 2) %>% pull(DayInt_SPC) %>% mean()
Val2_sum %>% filter(Trial == 3) %>% pull(DayInt_SPC) %>% mean()
Val2_sum %>% filter(Trial == 4) %>% pull(DayInt_SPC) %>% mean()
Val2_sum %>% filter(CVTA == "Y" & Storage_Temp == 10 & Day <= 21) %>% count() 
Val2_sum %>% filter(Storage_Temp == 10 & Day <= 21) %>% count()
Val2_sum %>% filter(CVTA == "Y" & Storage_Temp == 3 & Day <= 21) %>% count() 
Val2_sum %>% filter(Storage_Temp == 3 & Day <= 21) %>% count()
Val2_sum %>% filter(CVTA == "Y" & Storage_Temp == 6.5 & Day <= 21) %>% count() 
Val2_sum %>% filter(Storage_Temp == 6.5 & Day <= 21) %>% count()
Val2_sum %>% filter(Storage_Temp == 3 & Day == 0) %>% pull(APC) %>% mean()
Val2_sum %>% filter(Storage_Temp == 6.5 & Day == 0) %>% pull(APC) %>% mean()
Val2_sum %>% filter(Storage_Temp == 10 & Day == 0) %>% pull(APC) %>% mean()

#Data processing
Val2_data = Val2_data[Val2_data$CVTA=="N" & Val2_data$Day != 0 & Val2_data$Day != 17 & Val2_data$Day != 49 & Val2_data$Day != 63,]
Val2 = data.frame(rep(Val2_data$Day, 1000), rep(Val2_data$DayInt_SPC, 1000), rep(Val2_data$Storage_Temp,1000), rep(Val2_data$APC,1000))
names(Val2) = c("Day","DayInt_SPC","Storage_Temp","APC")

#Setting up growth rate and lag at different temp
Val2$AT = sample(LS_ATfreq, nrow(Val2), replace = T)
spore_growth_import$AT = c(3,100,15,513,61,45,179,340,17,23,159)

index = sapply(Val2$AT, function(x) which(x==spore_growth_import$AT))

#Val2$DayInt_SPC = runif(nrow(Val2), min = 0, max = 160) %>% round()

Val2$Mu = muAtNewTemp(Val2$Storage_Temp, spore_growth_import$mumax[index])
Val2$Lag = lagAtNewTemp(Val2$Storage_Temp, spore_growth_import$lag[index])
Val2$Nmax = spore_growth_import$LOG10Nmax[index]
Val2$Log10_DayInt = log10(Val2$DayInt_SPC + 1)

#Calculating growth 
Val2$Pred = log10N_func(Val2$Day, Val2$Lag, Val2$Mu, Val2$Log10_DayInt, Val2$Nmax)
Val2[which(is.na(Val2$Pred)),]$Pred = 0
Val2$Actual = log10(Val2$APC + 1)

#summary(Val2$Pred - Val2$Actual)
#sqrt(sum((Val2$Pred - Val2$Actual)^2)/nrow(Val2))


#Performance
Val2$Day = as.factor(Val2$Day)
Val2_tidy = pivot_longer(Val2, cols = 10:11, names_to = "Type", values_to = "Count")
Val2_T3 = Val2_tidy %>% filter(Storage_Temp == 3) %>% as.data.frame()
Val2_T6.5 = Val2_tidy %>% filter(Storage_Temp == 6.5) %>% as.data.frame()
Val2_T10 = Val2_tidy %>% filter(Storage_Temp == 10) %>% as.data.frame()

## RMSE for T3
T3D7.P = Val2_T3 %>% filter(Type == "Pred" & Day == 7) %>% pull(Count) %>% mean()
T3D7.A = Val2_T3 %>% filter(Type == "Actual" & Day == 7) %>% pull(Count) %>% mean()
T3D14.P = Val2_T3 %>% filter(Type == "Pred" & Day == 14) %>% pull(Count) %>% mean()
T3D14.A = Val2_T3 %>% filter(Type == "Actual" & Day == 14) %>% pull(Count) %>% mean()
T3D21.P = Val2_T3 %>% filter(Type == "Pred" & Day == 21) %>% pull(Count) %>% mean()
T3D21.A = Val2_T3 %>% filter(Type == "Actual" & Day == 21) %>% pull(Count) %>% mean()
T3D28.P = Val2_T3 %>% filter(Type == "Pred" & Day == 28) %>% pull(Count) %>% mean()
T3D28.A = Val2_T3 %>% filter(Type == "Actual" & Day == 28) %>% pull(Count) %>% mean()
T3D35.P = Val2_T3 %>% filter(Type == "Pred" & Day == 35) %>% pull(Count) %>% mean()
T3D35.A = Val2_T3 %>% filter(Type == "Actual" & Day == 35) %>% pull(Count) %>% mean()
T3D42.P = Val2_T3 %>% filter(Type == "Pred" & Day == 42) %>% pull(Count) %>% mean()
T3D42.A = Val2_T3 %>% filter(Type == "Actual" & Day == 42) %>% pull(Count) %>% mean()
T3.RMSE = sqrt(((T3D7.P-T3D7.A)^2 + 
                 (T3D14.P-T3D14.A)^2 +
                 (T3D21.P-T3D21.A)^2 +
                 (T3D28.P-T3D28.A)^2 +
                 (T3D35.P-T3D35.A)^2 +
                 (T3D42.P-T3D42.A)^2) /6)
T3.21.RMSE = sqrt(((T3D7.P-T3D7.A)^2 + 
                  (T3D14.P-T3D14.A)^2 +
                  (T3D21.P-T3D21.A)^2)/3)

## RMSE for T6.5
T6.5D7.P = Val2_T6.5 %>% filter(Type == "Pred" & Day == 7) %>% pull(Count) %>% mean()
T6.5D7.A = Val2_T6.5 %>% filter(Type == "Actual" & Day == 7) %>% pull(Count) %>% mean()
T6.5D14.P = Val2_T6.5 %>% filter(Type == "Pred" & Day == 14) %>% pull(Count) %>% mean()
T6.5D14.A = Val2_T6.5 %>% filter(Type == "Actual" & Day == 14) %>% pull(Count) %>% mean()
T6.5D21.P = Val2_T6.5 %>% filter(Type == "Pred" & Day == 21) %>% pull(Count) %>% mean()
T6.5D21.A = Val2_T6.5 %>% filter(Type == "Actual" & Day == 21) %>% pull(Count) %>% mean()
T6.5D28.P = Val2_T6.5 %>% filter(Type == "Pred" & Day == 28) %>% pull(Count) %>% mean()
T6.5D28.A = Val2_T6.5 %>% filter(Type == "Actual" & Day == 28) %>% pull(Count) %>% mean()
T6.5D35.P = Val2_T6.5 %>% filter(Type == "Pred" & Day == 35) %>% pull(Count) %>% mean()
T6.5D35.A = Val2_T6.5 %>% filter(Type == "Actual" & Day == 35) %>% pull(Count) %>% mean()
T6.5D42.P = Val2_T6.5 %>% filter(Type == "Pred" & Day == 42) %>% pull(Count) %>% mean()
T6.5D42.A = Val2_T6.5 %>% filter(Type == "Actual" & Day == 42) %>% pull(Count) %>% mean()

T6.5.RMSE = sqrt(((T6.5D7.P-T6.5D7.A)^2 + 
                  (T6.5D14.P-T6.5D14.A)^2 +
                  (T6.5D21.P-T6.5D21.A)^2 +
                  (T6.5D28.P-T6.5D28.A)^2 +
                  (T6.5D35.P-T6.5D35.A)^2 +
                  (T6.5D42.P-T6.5D42.A)^2)/6)
T6.5.21.RMSE = sqrt(((T6.5D7.P-T6.5D7.A)^2 + 
                    (T6.5D14.P-T6.5D14.A)^2 +
                    (T6.5D21.P-T6.5D21.A)^2)/3)

## RMSE for T10
T10D3.P = Val2_T10 %>% filter(Type == "Pred" & Day == 3) %>% pull(Count) %>% mean()
T10D3.A = Val2_T10 %>% filter(Type == "Actual" & Day == 3) %>% pull(Count) %>% mean()
T10D5.P = Val2_T10 %>% filter(Type == "Pred" & Day == 5) %>% pull(Count) %>% mean()
T10D5.A = Val2_T10 %>% filter(Type == "Actual" & Day == 5) %>% pull(Count) %>% mean()
T10D7.P = Val2_T10 %>% filter(Type == "Pred" & Day == 7) %>% pull(Count) %>% mean()
T10D7.A = Val2_T10 %>% filter(Type == "Actual" & Day == 7) %>% pull(Count) %>% mean()
T10D9.P = Val2_T10 %>% filter(Type == "Pred" & Day == 9) %>% pull(Count) %>% mean()
T10D9.A = Val2_T10 %>% filter(Type == "Actual" & Day == 9) %>% pull(Count) %>% mean()
T10D11.P = Val2_T10 %>% filter(Type == "Pred" & Day == 11) %>% pull(Count) %>% mean()
T10D11.A = Val2_T10 %>% filter(Type == "Actual" & Day == 11) %>% pull(Count) %>% mean()
T10D13.P = Val2_T10 %>% filter(Type == "Pred" & Day == 13) %>% pull(Count) %>% mean()
T10D13.A = Val2_T10 %>% filter(Type == "Actual" & Day == 13) %>% pull(Count) %>% mean()
T10D15.P = Val2_T10 %>% filter(Type == "Pred" & Day == 15) %>% pull(Count) %>% mean()
T10D15.A = Val2_T10 %>% filter(Type == "Actual" & Day == 15) %>% pull(Count) %>% mean()


T10.RMSE = sqrt(((T10D3.P-T10D3.A)^2 + 
                  (T10D5.P-T10D5.A)^2 +
                  (T10D7.P-T10D7.A)^2 +
                  (T10D9.P-T10D9.A)^2 +
                  (T10D11.P-T10D11.A)^2 +
                  (T10D13.P-T10D13.A)^2 +
                  (T10D15.P-T10D15.A)^2)/7)
## Overall
(All.RMSE = sqrt(((T3D7.P-T3D7.A)^2 + 
                  (T3D14.P-T3D14.A)^2 +
                  (T3D21.P-T3D21.A)^2 +
                  (T3D28.P-T3D28.A)^2 +
                  (T3D35.P-T3D35.A)^2 +
                  (T3D42.P-T3D42.A)^2 +
                  (T6.5D7.P-T6.5D7.A)^2 + 
                  (T6.5D14.P-T6.5D14.A)^2 +
                  (T6.5D21.P-T6.5D21.A)^2 +
                  (T6.5D28.P-T6.5D28.A)^2 +
                  (T6.5D35.P-T6.5D35.A)^2 +
                  (T6.5D42.P-T6.5D42.A)^2+
                  (T10D3.P-T10D3.A)^2 + 
                  (T10D5.P-T10D5.A)^2 +
                  (T10D7.P-T10D7.A)^2 +
                  (T10D9.P-T10D9.A)^2 +
                  (T10D11.P-T10D11.A)^2 +
                  (T10D13.P-T10D13.A)^2 +
                  (T10D15.P-T10D15.A)^2/22)))


#Visualization
ggplot(data = Val2_T3, aes(x=Day, y=Count, fill=Type))+
  geom_boxplot( position = position_dodge(0.9),
                show.legend = FALSE)+
  scale_fill_manual(values = c("white", "grey"))+
  #facet_wrap(~Type)+
  theme_classic()+
  #scale_x_discrete(breaks = sort(unique(Val2_T3$Day)))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Predicted concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (3°C) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))

ggplot(data = Val2_T6.5, aes(x=Day, y=Count, fill=Type))+
  geom_boxplot( position = position_dodge(0.9),
                show.legend = F)+
  scale_fill_manual(values = c("white", "grey"))+
  #facet_wrap(~Type)+
  theme_classic()+
  #scale_x_discrete(breaks = sort(unique(Val2_T3$Day)))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Predicted concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (6.5°C) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))

ggplot(data = Val2_T10, aes(x=Day, y=Count, fill=Type))+
  geom_boxplot( position = position_dodge(0.9),
                show.legend = FALSE)+
  scale_fill_manual(values = c("white", "grey"))+
  #facet_wrap(~Type)+
  theme_classic()+
  #scale_x_discrete(breaks = sort(unique(Val2_T3$Day)))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Predicted concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (10°C) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))
