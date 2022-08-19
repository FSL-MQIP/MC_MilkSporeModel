third_val = read.csv("InputFiles/ValidationData3.csv", header = TRUE)

# Summary stats
third_val %>% count(PPC == "Y")
third_val %>% filter(Storage_Temp == 4 & Day == 0) %>% pull(APC) %>% as.numeric() %>% mean()
third_val %>% filter(Storage_Temp == 6 & Day == 0) %>% pull(APC) %>% as.numeric() %>% mean()
third_val %>% filter(Storage_Temp == "SC" & Day == 0) %>% pull(APC) %>% as.numeric() %>% mean()
third_val %>% filter(Storage_Temp == 4 & Day <= 21) %>% count()
third_val %>% filter(Storage_Temp == 4 & Day <= 21 & PPC != "Y" & APC != "TNTC") %>% count()
third_val %>% filter(Storage_Temp == 6 & Day <= 21) %>% count()
third_val %>% filter(Storage_Temp == 6 & Day <= 21 & PPC != "Y" & APC != "TNTC") %>% count()
third_val %>% filter(Storage_Temp == "SC" & Day <= 21) %>% count()
third_val %>% filter(Storage_Temp == "SC" & Day <= 21 & PPC != "Y" & APC != "TNTC") %>% count()

# Subset validation data
third_val = third_val[which(third_val$PPC != "Y" & third_val$APC != "TNTC"),]
third_val_T4D0 = third_val[which(third_val$Day==0 & third_val$Storage_Temp==4), 5] %>% as.numeric() %>% log10()
third_val_T6D0 = third_val[which(third_val$Day==0 && third_val$Storage_Temp==6), 5] %>% as.numeric() %>% log10()
third_val_SCD0 = third_val[which(third_val$Day==0 & third_val$Storage_Temp=="SC"), 5] %>% as.numeric() %>% log10()

third_val_T4D7 = third_val[which(third_val$Day==7 & third_val$Storage_Temp==4), 5] %>% as.numeric() %>% log10()
third_val_T6D7 = third_val[which(third_val$Day==7 & third_val$Storage_Temp==6), 5] %>% as.numeric() %>% log10()
third_val_SCD7 = third_val[which(third_val$Day==7 & third_val$Storage_Temp=="SC"), 5] %>% as.numeric() %>% log10()

third_val_T4D14 = third_val[which(third_val$Day==14 & third_val$Storage_Temp==4), 5] %>% as.numeric() %>% log10()
third_val_T6D14 = third_val[which(third_val$Day==14 & third_val$Storage_Temp==6), 5] %>% as.numeric() %>% log10()
third_val_SCD14 = third_val[which(third_val$Day==14 & third_val$Storage_Temp=="SC"), 5] %>% as.numeric() %>% log10()

third_val_T4D21 = third_val[which(third_val$Day==21 & third_val$Storage_Temp==4), 5] %>% as.numeric() %>% log10()
third_val_T6D21 = third_val[which(third_val$Day==21 & third_val$Storage_Temp==6), 5] %>% as.numeric() %>% log10()
third_val_SCD21 = third_val[which(third_val$Day==21 & third_val$Storage_Temp=="SC"), 5] %>% as.numeric() %>% log10()

third_val_T4D28 = third_val[which(third_val$Day==28 & third_val$Storage_Temp==4), 5] %>% as.numeric() %>% log10()
third_val_T6D28 = third_val[which(third_val$Day==28 & third_val$Storage_Temp==6), 5] %>% as.numeric() %>% log10()
third_val_SCD28 = third_val[which(third_val$Day==28 & third_val$Storage_Temp=="SC"), 5] %>% as.numeric() %>% log10()

third_val_T4D35 = third_val[which(third_val$Day==35 & third_val$Storage_Temp==4), 5] %>% as.numeric() %>% log10()
third_val_T6D35 = third_val[which(third_val$Day==35 & third_val$Storage_Temp==6), 5] %>% as.numeric() %>% log10()
third_val_SCD35 = third_val[which(third_val$Day==35 & third_val$Storage_Temp=="SC"), 5] %>% as.numeric() %>% log10()

third_val_T4D42 = third_val[which(third_val$Day==42 & third_val$Storage_Temp==4), 5] %>% as.numeric() %>% log10()
third_val_T6D42 = third_val[which(third_val$Day==42 & third_val$Storage_Temp==6), 5] %>% as.numeric() %>% log10()
third_val_SCD42 = third_val[which(third_val$Day==42 & third_val$Storage_Temp=="SC"), 5] %>% as.numeric() %>% log10()

################################################################## day 7 #########################################################
# Day 7 Temp = 4C
#third_val_T4D0 = third_val[c(1,2,5,6,29,30), 5] %>% as.numeric() %>% log10()
third_val_T4D0 = third_val %>% 
                    filter(Storage_Temp == 4 & Day == 0) %>% 
                    pull(APC) %>% as.numeric() %>% log10()
T4D0 = rep(third_val_T4D0, 1000)
sim_T4D7 = vector()
newLag_T4 = lagAtNewTemp(4, spore_growth_import$lag)
newMu_T4 =  muAtNewTemp(4, spore_growth_import$mumax)

AT = SampleAT(T4D0)

for (i in 1:length(T4D0)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_T4D7[i] = log10N_func(7, newLag_T4[allele_index],newMu_T4[allele_index],T4D0[i],spore_growth_import$LOG10Nmax[allele_index])
}
#summary(sim_T4D7)
#summary(third_val_T4D7)
#boxplot(sim_T4D7, third_val_T4D7, main = "4C at Day 7",names = c("Pred","Actual"))

# Day 7 Temp = 6C
#third_val_T6D0 = third_val[c(11,12,15,16,39,40), 5] %>% as.numeric() %>% log10()
third_val_T6D0 = third_val %>% 
                    filter(Storage_Temp == 6 & Day == 0) %>% 
                    pull(APC) %>% as.numeric() %>% log10()
T6D0 = rep(third_val_T6D0, 1000)
sim_T6D7 = vector()
newLag_T6 = lagAtNewTemp(6, spore_growth_import$lag)
newMu_T6 =  muAtNewTemp(6, spore_growth_import$mumax)

AT = SampleAT(T6D0)

for (i in 1:length(T6D0)){
  AT = sample(AT_freq, length(T6D0),T)
  allele_index = which(spore_growth_import$STorAT == AT[i])
  
  sim_T6D7[i] = log10N_func(7, newLag_T6[allele_index],newMu_T6[allele_index],T6D0[i],spore_growth_import$LOG10Nmax[allele_index])
}
#summary(sim_T6D7)
#summary(third_val_T6D7)
#boxplot(sim_T6D7, third_val_T6D7, main = "6C at Day 7",names = c("Pred","Actual"))

# Day 7 Temp = SC
#third_val_SCD0 = third_val[c(19,20,21,22,47,48),5] %>% as.numeric() %>% log10()
third_val_SCD0 = third_val %>% 
                    filter(Storage_Temp == "SC" & Day == 0) %>% 
                    pull(APC) %>% as.numeric() %>% log10()
SCD0 = rep(third_val_SCD0, 1000)
AT = SampleAT(SCD0)

newLag_T2 = lagAtNewTemp(2, spore_growth_import$lag)
newMu_T2 =  muAtNewTemp(2, spore_growth_import$mumax)
newLag_T10 = lagAtNewTemp(10, spore_growth_import$lag)
newMu_T10 =  muAtNewTemp(10, spore_growth_import$mumax)

sim_SCD7 = vector()
sim_SCS1 = vector()
remain_lag = vector()
remain_prop = vector()
for (i in 1:length(SCD0)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCS1[i] = log10N_func(41/24, newLag_T4[allele_index],newMu_T4[allele_index],SCD0[i],spore_growth_import$LOG10Nmax[allele_index])
  remain_lag[i] = ifelse(newLag_T4[allele_index]>=41/24, newLag_T4[allele_index]-41/24,0)
  remain_prop[i] = remain_lag[i]/newLag_T4[allele_index]
}

sim_SCS2 = vector()
for (i in 1:length(sim_SCS1)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCS2[i] = log10N_func(1, newLag_T2[allele_index]*remain_prop[i],newMu_T2[allele_index],sim_SCS1[i],spore_growth_import$LOG10Nmax[allele_index])
  remain_lag[i] = ifelse(newLag_T2[allele_index]*remain_prop[i]>=1, newLag_T2[allele_index]*remain_prop[i]-1,0)
  remain_prop[i] = remain_lag[i]/newLag_T4[allele_index]
}

sim_SCS3 = vector()
for (i in 1:length(sim_SCS2)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCS3[i] = log10N_func(26/60/24, newLag_T10[allele_index]*remain_prop[i], newMu_T10[allele_index],sim_SCS2[i],spore_growth_import$LOG10Nmax[allele_index])
  remain_lag[i] = ifelse(newLag_T10[allele_index]*remain_prop[i]>=26/60/24, newLag_T10[allele_index]*remain_prop[i]-26/60/24,0)
  remain_prop[i] = remain_lag[i]/newLag_T4[allele_index]
}

for (i in 1:length(sim_SCS3)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCD7[i] = log10N_func(7-41/24-1-26/60/24, newLag_T4[allele_index]*remain_prop[i], newMu_T4[allele_index],sim_SCS3[i],spore_growth_import$LOG10Nmax[allele_index])
}

#summary(sim_SCD7)
#summary(third_val_SCD7)
#boxplot(sim_SCD7, third_val_SCD7, main = "SC at Day 7",names = c("Pred","Actual"))

# Combined day 7
sim_D7 = c(sim_T4D7,sim_T6D7, sim_SCD7)
third_val_D7 = c(third_val_T4D7,third_val_T6D7, third_val_SCD7)
#ks.test(sim_D7,third_val_D7)
boxplot(sim_D7, third_val_D7, main = "Validation at Day 7",names = c("Pred","Actual"))

# Visualization
pred_D7 = cbind(rep("pred", length(sim_D7)),sim_D7)
pred_D7 = as.data.frame(pred_D7)
names(pred_D7) = c("type","counts")
actual_D7 = cbind(rep("actual", length(third_val_D7)), third_val_D7)
actual_D7 = as.data.frame(actual_D7)
names(actual_D7) = c("type","counts")
val_D7 = rbind(pred_D7,actual_D7)
val_D7$counts = as.numeric(val_D7$counts)

ggplot(data = val_D7, aes(x=type,y=counts))+
  geom_boxplot()+
  coord_flip()+
  theme_bw()

ggplot(data = val_D7, aes(x=counts))+
  stat_ecdf(aes(linetype=type))+
  theme_bw()+
  labs(x=expression(paste(LOG[10],CFU/mL)), y= "Density")



################################################################## day 14 #########################################################
# Day 14 Temp = 4C
#third_val_T4D0 = third_val[c(3,4,31,32), 5] %>% as.numeric() %>% log10()
#T4D0 = rep(third_val_T4D0, 1000)
#AT = SampleAT(T4D0)

sim_T4D14 = vector()
for (i in 1:length(T4D0)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_T4D14[i] = log10N_func(14, newLag_T4[allele_index],newMu_T4[allele_index],T4D0[i],spore_growth_import$LOG10Nmax[allele_index])
}
#summary(sim_T4D14)
#summary(third_val_T4D14)
#boxplot(sim_T4D14, third_val_T4D14, main = "4C at Day 14",names = c("Pred","Actual"))


# Day 14 Temp = 6C
#third_val_T6D0 = third_val[c(11,12,39,40), 5] %>% as.numeric() %>% log10()
#T6D0 = rep(third_val_T6D0, 1000)
#AT = SampleAT(T6D0)

sim_T6D14 = vector()
for (i in 1:length(T6D0)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  
  sim_T6D14[i] = log10N_func(14, newLag_T6[allele_index],newMu_T6[allele_index],T6D0[i],spore_growth_import$LOG10Nmax[allele_index])
}
#summary(sim_T6D14)
#summary(third_val_T6D14)
#boxplot(sim_T6D14, third_val_T6D14, main = "6C at Day 14",names = c("Pred","Actual"))

# Day 14 Temp = SC
#third_val_SCD0 = third_val[c(19,20,21,22,23,24,47,48),5] %>% as.numeric() %>% log10()
#SCD0 = rep(third_val_SCD0, 1000)
#AT = SampleAT(SCD0)

sim_SCD14 = vector()
sim_SCS1 = vector()
remain_lag = vector()
remain_prop = vector()
for (i in 1:length(SCD0)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCS1[i] = log10N_func(41/24, newLag_T4[allele_index],newMu_T4[allele_index],SCD0[i],spore_growth_import$LOG10Nmax[allele_index])
  remain_lag[i] = ifelse(newLag_T4[allele_index]>=41/24, newLag_T4[allele_index]-41/24,0)
  remain_prop[i] = remain_lag[i]/newLag_T4[allele_index]
}

sim_SCS2 = vector()
for (i in 1:length(sim_SCS1)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCS2[i] = log10N_func(1, newLag_T2[allele_index]*remain_prop[i],newMu_T2[allele_index],sim_SCS1[i],spore_growth_import$LOG10Nmax[allele_index])
  remain_lag[i] = ifelse(newLag_T2[allele_index]*remain_prop[i]>=1, newLag_T2[allele_index]*remain_prop[i]-1,0)
  remain_prop[i] = remain_lag[i]/newLag_T4[allele_index]
}

sim_SCS3 = vector()
for (i in 1:length(sim_SCS2)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCS3[i] = log10N_func(26/60/24, newLag_T10[allele_index]*remain_prop[i], newMu_T10[allele_index],sim_SCS2[i],spore_growth_import$LOG10Nmax[allele_index])
  remain_lag[i] = ifelse(newLag_T10[allele_index]*remain_prop[i]>=26/60/24, newLag_T10[allele_index]*remain_prop[i]-26/60/24,0)
  remain_prop[i] = remain_lag[i]/newLag_T4[allele_index]
}

for (i in 1:length(sim_SCS3)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCD14[i] = log10N_func(14-41/24-1-26/60/24, newLag_T4[allele_index]*remain_prop[i], newMu_T4[allele_index],sim_SCS3[i],spore_growth_import$LOG10Nmax[allele_index])
}

#summary(sim_SCD14)
#summary(third_val_SCD14)
#boxplot(sim_SCD14, third_val_SCD14, main = "SC at Day 14",names = c("Pred","Actual"))

# Combined day 14
sim_D14 = c(sim_T4D14,sim_T6D14, sim_SCD14)
third_val_D14 = c(third_val_T4D14,third_val_T6D14, third_val_SCD14)
#ks.test(sim_D14,third_val_D14)
boxplot(sim_D14, third_val_D14, main = "Validation at Day 14",names = c("Pred","Actual"))

# Visualization
pred_D14 = cbind(rep("pred", length(sim_D14)),sim_D14)
pred_D14 = as.data.frame(pred_D14)
names(pred_D14) = c("type","counts")
actual_D14 = cbind(rep("actual", length(third_val_D14)), third_val_D14)
actual_D14 = as.data.frame(actual_D14)
names(actual_D14) = c("type","counts")
val_D14 = rbind(pred_D14,actual_D14)
val_D14$counts = as.numeric(val_D14$counts)

ggplot(data = val_D14, aes(x=type,y=counts))+
  geom_boxplot()+
  coord_flip()+
  theme_bw()

ggplot(data = val_D14, aes(x=counts))+
  stat_ecdf(aes(linetype=type))+
  theme_bw()+
  labs(x=expression(paste(LOG[10],CFU/mL)), y= "Density")

################################################################## day 21 #########################################################
# Day 21 Temp = 4C
#third_val_T4D0 = third_val[c(1,2,3,4,31,32), 5] %>% as.numeric() %>% log10()
#T4D0 = rep(third_val_T4D0, 1000)
#AT = SampleAT(T4D0)

sim_T4D21 = vector()
for (i in 1:length(T4D0)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_T4D21[i] = log10N_func(21, newLag_T4[allele_index],newMu_T4[allele_index],T4D0[i],spore_growth_import$LOG10Nmax[allele_index])
}
#summary(sim_T4D21)
#boxplot(sim_T4D21, third_val_T4D21, main = "4C at Day 21",names = c("Pred","Actual"))

# Day 21 Temp = 6C
#third_val_T6D0 = third_val[c(9,10,11,12,39,40), 5] %>% as.numeric() %>% log10()
#T6D0 = rep(third_val_T6D0, 100)
#AT = SampleAT(T6D0)

sim_T6D21 = vector()
for (i in 1:length(T6D0)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_T6D21[i] = log10N_func(21, newLag_T6[allele_index],newMu_T6[allele_index],T6D0[i],spore_growth_import$LOG10Nmax[allele_index])
}
#summary(sim_T6D21)
#summary(third_val_T6D21)
#boxplot(sim_T6D21, third_val_T6D21, main = "6C at Day 21",names = c("Pred","Actual"))


# Day 21 Temp = SC
#third_val_SCD0 = third_val[c(19,20,21,22,47,48),5] %>% as.numeric() %>% log10()
#SCD0 = rep(third_val_SCD0, 1000)
#AT = SampleAT(SCD0)

sim_SCD21 = vector()
sim_SCS1 = vector()
remain_lag = vector()
remain_prop = vector()
for (i in 1:length(SCD0)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCS1[i] = log10N_func(41/24, newLag_T4[allele_index],newMu_T4[allele_index],SCD0[i],spore_growth_import$LOG10Nmax[allele_index])
  remain_lag[i] = ifelse(newLag_T4[allele_index]>=41/24, newLag_T4[allele_index]-41/24,0)
  remain_prop[i] = remain_lag[i]/newLag_T4[allele_index]
}

sim_SCS2 = vector()
for (i in 1:length(sim_SCS1)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCS2[i] = log10N_func(1, newLag_T2[allele_index]*remain_prop[i],newMu_T2[allele_index],sim_SCS1[i],spore_growth_import$LOG10Nmax[allele_index])
  remain_lag[i] = ifelse(newLag_T2[allele_index]*remain_prop[i]>=1, newLag_T2[allele_index]*remain_prop[i]-1,0)
  remain_prop[i] = remain_lag[i]/newLag_T4[allele_index]
}

sim_SCS3 = vector()
for (i in 1:length(sim_SCS2)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCS3[i] = log10N_func(26/60/24, newLag_T10[allele_index]*remain_prop[i], newMu_T10[allele_index],sim_SCS2[i],spore_growth_import$LOG10Nmax[allele_index])
  remain_lag[i] = ifelse(newLag_T10[allele_index]*remain_prop[i]>=26/60/24, newLag_T10[allele_index]*remain_prop[i]-26/60/24,0)
  remain_prop[i] = remain_lag[i]/newLag_T4[allele_index]
}

for (i in 1:length(sim_SCS3)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCD21[i] = log10N_func(21-41/24-1-26/60/24, newLag_T4[allele_index]*remain_prop[i], newMu_T4[allele_index],sim_SCS3[i],spore_growth_import$LOG10Nmax[allele_index])
}

#summary(sim_SCD21)
#summary(third_val_SCD21)
#boxplot(sim_SCD21, third_val_SCD21, main = "SC at Day 21",names = c("Pred","Actual"))

# Combined day 21
sim_D21 = c(sim_T4D21,sim_T6D21, sim_SCD21)
third_val_D21 = c(third_val_T4D21,third_val_T6D21, third_val_SCD21)
#ks.test(sim_D21,third_val_D21)
boxplot(sim_D21, third_val_D21, main = "Validation at Day 21",names = c("Pred","Actual"))

# Visualization
pred_D7 = cbind(rep("pred", length(sim_D7)),sim_D7)
pred_D7 = as.data.frame(pred_D7)
names(pred_D7) = c("type","counts")
actual_D7 = cbind(rep("actual", length(third_val_D7)), third_val_D7)
actual_D7 = as.data.frame(actual_D7)
names(actual_D7) = c("type","counts")
val_D7 = rbind(pred_D7,actual_D7)
val_D7$counts = as.numeric(val_D7$counts)

ggplot(data = val_D7, aes(x=type,y=counts))+
  geom_boxplot()+
  coord_flip()+
  theme_bw()

ggplot(data = val_D7, aes(x=counts))+
  stat_ecdf(aes(linetype=type))+
  theme_bw()+
  labs(x=expression(paste(LOG[10],CFU/mL)), y= "Density")

################################################################## day 28 #########################################################
# Day 28 Temp = 4C


# Day 28 Temp = 6C
#third_val_T6D0 = third_val[c(9,10,11,12,15,16), 5] %>% as.numeric() %>% log10()
#T6D0 = rep(third_val_T6D0, 100)
#AT = SampleAT(T6D0)

sim_T6D28 = vector()
for (i in 1:length(T6D0)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_T6D28[i] = log10N_func(28, newLag_T6[allele_index],newMu_T6[allele_index],T6D0[i],spore_growth_import$LOG10Nmax[allele_index])
}
#summary(sim_T6D28)
#summary(third_val_T6D28)
#boxplot(sim_T6D28, third_val_T6D28, main = "6C at Day 28",names = c("Pred","Actual"))


# Day 28 Temp = SC
#third_val_SCD0 = third_val[c(19,20,21,22,23,24),5] %>% as.numeric() %>% log10()
#SCD0 = rep(third_val_SCD0, 1000)
#AT = SampleAT(SCD0)

sim_SCD28 = vector()
sim_SCS1 = vector()
remain_lag = vector()
remain_prop = vector()
for (i in 1:length(SCD0)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCS1[i] = log10N_func(41/24, newLag_T4[allele_index],newMu_T4[allele_index],SCD0[i],spore_growth_import$LOG10Nmax[allele_index])
  remain_lag[i] = ifelse(newLag_T4[allele_index]>=41/24, newLag_T4[allele_index]-41/24,0)
  remain_prop[i] = remain_lag[i]/newLag_T4[allele_index]
}

sim_SCS2 = vector()
for (i in 1:length(sim_SCS1)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCS2[i] = log10N_func(1, newLag_T2[allele_index]*remain_prop[i],newMu_T2[allele_index],sim_SCS1[i],spore_growth_import$LOG10Nmax[allele_index])
  remain_lag[i] = ifelse(newLag_T2[allele_index]*remain_prop[i]>=1, newLag_T2[allele_index]*remain_prop[i]-1,0)
  remain_prop[i] = remain_lag[i]/newLag_T4[allele_index]
}

sim_SCS3 = vector()
for (i in 1:length(sim_SCS2)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCS3[i] = log10N_func(26/60/24, newLag_T10[allele_index]*remain_prop[i], newMu_T10[allele_index],sim_SCS2[i],spore_growth_import$LOG10Nmax[allele_index])
  remain_lag[i] = ifelse(newLag_T10[allele_index]*remain_prop[i]>=26/60/24, newLag_T10[allele_index]*remain_prop[i]-26/60/24,0)
  remain_prop[i] = remain_lag[i]/newLag_T4[allele_index]
}

for (i in 1:length(sim_SCS3)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCD28[i] = log10N_func(28-41/24-1-26/60/24, newLag_T4[allele_index]*remain_prop[i], newMu_T4[allele_index],sim_SCS3[i],spore_growth_import$LOG10Nmax[allele_index])
}

#summary(sim_SCD28)
#summary(third_val_SCD28)
#boxplot(sim_SCD28, third_val_SCD28, main = "SC at Day 7",names = c("Pred","Actual"))

# Combined day 28
sim_D28 = c(sim_T6D28, sim_SCD28)
third_val_D28 = c(third_val_T6D28, third_val_SCD28)
#ks.test(sim_D28,third_val_D28)
boxplot(sim_D28, third_val_D28, main = "Validation at Day 28",names = c("Pred","Actual"))

# Visualization
pred_D7 = cbind(rep("pred", length(sim_D7)),sim_D7)
pred_D7 = as.data.frame(pred_D7)
names(pred_D7) = c("type","counts")
actual_D7 = cbind(rep("actual", length(third_val_D7)), third_val_D7)
actual_D7 = as.data.frame(actual_D7)
names(actual_D7) = c("type","counts")
val_D7 = rbind(pred_D7,actual_D7)
val_D7$counts = as.numeric(val_D7$counts)

ggplot(data = val_D7, aes(x=type,y=counts))+
  geom_boxplot()+
  coord_flip()+
  theme_bw()

ggplot(data = val_D7, aes(x=counts))+
  stat_ecdf(aes(linetype=type))+
  theme_bw()+
  labs(x=expression(paste(LOG[10],CFU/mL)), y= "Density")


################################################################## day 35 #########################################################
# Day 35 Temp = 4C
#third_val_T4D0 = third_val[c(1,2), 5] %>% as.numeric() %>% log10()
#T4D0 = rep(third_val_T4D0, 100)
#AT = SampleAT(T4D0)

sim_T4D35 = vector()
for (i in 1:length(T4D0)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_T4D35[i] = log10N_func(35, newLag_T4[allele_index],newMu_T4[allele_index],T4D0[i],spore_growth_import$LOG10Nmax[allele_index])
}
#summary(sim_T4D35)
#summary(third_val_T4D35)
#boxplot(sim_T4D35, third_val_T4D35, main = "4C at Day 35",names = c("Pred","Actual"))

# Day 35 Temp = 6C
#third_val_T6D0 = third_val[c(39,40), 5] %>% as.numeric() %>% log10()
#T6D0 = rep(third_val_T6D0, 100)
#AT = SampleAT(T6D0)

sim_T6D35 = vector()
for (i in 1:length(T6D0)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_T6D35[i] = log10N_func(35, newLag_T6[allele_index],newMu_T6[allele_index],T6D0[i],spore_growth_import$LOG10Nmax[allele_index])
}
#summary(sim_T6D35)
#summary(third_val_T6D35)
#boxplot(sim_T6D35, third_val_T6D35, main = "6C at Day 35",names = c("Pred","Actual"))


# Day 35 Temp = SC
#third_val_SCD0 = third_val[c(19,20,21,22,47,48),5] %>% as.numeric() %>% log10()
#SCD0 = rep(third_val_SCD0, 1000)
#AT = SampleAT(SCD0)

sim_SCD35 = vector()
sim_SCS1 = vector()
remain_lag = vector()
remain_prop = vector()
for (i in 1:length(SCD0)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCS1[i] = log10N_func(41/24, newLag_T4[allele_index],newMu_T4[allele_index],SCD0[i],spore_growth_import$LOG10Nmax[allele_index])
  remain_lag[i] = ifelse(newLag_T4[allele_index]>=41/24, newLag_T4[allele_index]-41/24,0)
  remain_prop[i] = remain_lag[i]/newLag_T4[allele_index]
}

sim_SCS2 = vector()
for (i in 1:length(sim_SCS1)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCS2[i] = log10N_func(1, newLag_T2[allele_index]*remain_prop[i],newMu_T2[allele_index],sim_SCS1[i],spore_growth_import$LOG10Nmax[allele_index])
  remain_lag[i] = ifelse(newLag_T2[allele_index]*remain_prop[i]>=1, newLag_T2[allele_index]*remain_prop[i]-1,0)
  remain_prop[i] = remain_lag[i]/newLag_T4[allele_index]
}

sim_SCS3 = vector()
for (i in 1:length(sim_SCS2)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCS3[i] = log10N_func(26/60/24, newLag_T10[allele_index]*remain_prop[i], newMu_T10[allele_index],sim_SCS2[i],spore_growth_import$LOG10Nmax[allele_index])
  remain_lag[i] = ifelse(newLag_T10[allele_index]*remain_prop[i]>=26/60/24, newLag_T10[allele_index]*remain_prop[i]-26/60/24,0)
  remain_prop[i] = remain_lag[i]/newLag_T4[allele_index]
}

for (i in 1:length(sim_SCS3)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCD35[i] = log10N_func(35-41/24-1-26/60/24, newLag_T4[allele_index]*remain_prop[i], newMu_T4[allele_index],sim_SCS3[i],spore_growth_import$LOG10Nmax[allele_index])
}

#summary(sim_SCD35)
#summary(third_val_SCD35)
#boxplot(sim_SCD35, third_val_SCD35, main = "SC at Day 7",names = c("Pred","Actual"))

# Combined day 35
sim_D35 = c(sim_T4D35,sim_T6D35, sim_SCD35)
third_val_D35 = c(third_val_T4D35,third_val_T6D35, third_val_SCD35)
#ks.test(sim_D35,third_val_D35)
boxplot(sim_D35, third_val_D35, main = "Validation at Day 7",names = c("Pred","Actual"))

# Visualization
pred_D7 = cbind(rep("pred", length(sim_D7)),sim_D7)
pred_D7 = as.data.frame(pred_D7)
names(pred_D7) = c("type","counts")
actual_D7 = cbind(rep("actual", length(third_val_D7)), third_val_D7)
actual_D7 = as.data.frame(actual_D7)
names(actual_D7) = c("type","counts")
val_D7 = rbind(pred_D7,actual_D7)
val_D7$counts = as.numeric(val_D7$counts)

ggplot(data = val_D7, aes(x=type,y=counts))+
  geom_boxplot()+
  coord_flip()+
  theme_bw()

ggplot(data = val_D7, aes(x=counts))+
  stat_ecdf(aes(linetype=type))+
  theme_bw()+
  labs(x=expression(paste(LOG[10],CFU/mL)), y= "Density")

################################################################## day 42 #########################################################
# Day 42 Temp = 4C
#third_val_T4D0 = third_val[c(1,2,3,4,7,8,25,26), 5] %>% as.numeric() %>% log10()
#T4D0 = rep(third_val_T4D0, 1000)
#AT = SampleAT(T4D0)

sim_T4D42 = vector()
for (i in 1:length(T4D0)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_T4D42[i] = log10N_func(42, newLag_T4[allele_index],newMu_T4[allele_index],T4D0[i],spore_growth_import$LOG10Nmax[allele_index])
}
#summary(sim_T4D42)
#summary(third_val_T4D42)
#boxplot(sim_T4D42, third_val_T4D42, main = "4C at Day 42",names = c("Pred","Actual"))

# Day 42 Temp = 6C


# Day 42 Temp = SC
#third_val_SCD0 = third_val[c(19,20,41,42),5] %>% as.numeric() %>% log10()
#SCD0 = rep(third_val_SCD0, 1000)
#AT = SampleAT(SCD0)

sim_SCD42 = vector()
sim_SCS1 = vector()
remain_lag = vector()
remain_prop = vector()
for (i in 1:length(SCD0)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCS1[i] = log10N_func(41/24, newLag_T4[allele_index],newMu_T4[allele_index],SCD0[i],spore_growth_import$LOG10Nmax[allele_index])
  remain_lag[i] = ifelse(newLag_T4[allele_index]>=41/24, newLag_T4[allele_index]-41/24,0)
  remain_prop[i] = remain_lag[i]/newLag_T4[allele_index]
}

sim_SCS2 = vector()
for (i in 1:length(sim_SCS1)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCS2[i] = log10N_func(1, newLag_T2[allele_index]*remain_prop[i],newMu_T2[allele_index],sim_SCS1[i],spore_growth_import$LOG10Nmax[allele_index])
  remain_lag[i] = ifelse(newLag_T2[allele_index]*remain_prop[i]>=1, newLag_T2[allele_index]*remain_prop[i]-1,0)
  remain_prop[i] = remain_lag[i]/newLag_T4[allele_index]
}

sim_SCS3 = vector()
for (i in 1:length(sim_SCS2)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCS3[i] = log10N_func(26/60/24, newLag_T10[allele_index]*remain_prop[i], newMu_T10[allele_index],sim_SCS2[i],spore_growth_import$LOG10Nmax[allele_index])
  remain_lag[i] = ifelse(newLag_T10[allele_index]*remain_prop[i]>=26/60/24, newLag_T10[allele_index]*remain_prop[i]-26/60/24,0)
  remain_prop[i] = remain_lag[i]/newLag_T4[allele_index]
}

for (i in 1:length(sim_SCS3)){
  allele_index = which(spore_growth_import$STorAT == AT[i])
  sim_SCD42[i] = log10N_func(42-41/24-1-26/60/24, newLag_T4[allele_index]*remain_prop[i], newMu_T4[allele_index],sim_SCS3[i],spore_growth_import$LOG10Nmax[allele_index])
}

#summary(sim_SCD42)
#summary(third_val_SCD42)
#boxplot(sim_SCD42, third_val_SCD42, main = "SC at Day 7",names = c("Pred","Actual"))


# Combined day 42
sim_D42 = c(sim_T4D42, sim_SCD42)
third_val_D42 = c(third_val_T4D42, third_val_SCD42)
#ks.test(sim_D42,third_val_D42)
boxplot(sim_D42, third_val_D42, main = "Validation at Day 42",names = c("Pred","Actual"))

# Visualization
pred_D7 = cbind(rep("pred", length(sim_D7)),sim_D7)
pred_D7 = as.data.frame(pred_D7)
names(pred_D7) = c("type","counts")
actual_D7 = cbind(rep("actual", length(third_val_D7)), third_val_D7)
actual_D7 = as.data.frame(actual_D7)
names(actual_D7) = c("type","counts")
val_D7 = rbind(pred_D7,actual_D7)
val_D7$counts = as.numeric(val_D7$counts)

ggplot(data = val_D7, aes(x=type,y=counts))+
  geom_boxplot()+
  coord_flip()+
  theme_bw()

ggplot(data = val_D7, aes(x=counts))+
  stat_ecdf(aes(linetype=type))+
  theme_bw()+
  labs(x=expression(paste(LOG[10],CFU/mL)), y= "Density")


#################################################################################
sum_stats = matrix(c(ks.test(sim_D7,third_val_D7)$statistic, ks.test(sim_D7,third_val_D7)$p.value, median(sim_D7),median(third_val_D7),mean(sim_D7),mean(third_val_D7),sd(sim_D7),sd(third_val_D7),
                     ks.test(sim_D14,third_val_D14)$statistic, ks.test(sim_D14,third_val_D14)$p.value, median(sim_D14),median(third_val_D14),mean(sim_D14),mean(third_val_D14),sd(sim_D14),sd(third_val_D14),
                     ks.test(sim_D21,third_val_D21)$statistic, ks.test(sim_D21,third_val_D21)$p.value, median(sim_D21),median(third_val_D21),mean(sim_D21),mean(third_val_D21),sd(sim_D21),sd(third_val_D21),
                     ks.test(sim_D28,third_val_D28)$statistic, ks.test(sim_D28,third_val_D28)$p.value, median(sim_D28),median(third_val_D28),mean(sim_D28),mean(third_val_D28),sd(sim_D28),sd(third_val_D28),
                     ks.test(sim_D35,third_val_D35)$statistic, ks.test(sim_D35,third_val_D35)$p.value, median(sim_D35),median(third_val_D35),mean(sim_D35),mean(third_val_D35),sd(sim_D35),sd(third_val_D35),
                     ks.test(sim_D42,third_val_D42)$statistic, ks.test(sim_D42,third_val_D42)$p.value, median(sim_D42),median(third_val_D42),mean(sim_D42),mean(third_val_D42),sd(sim_D42),sd(third_val_D42)),
                   ncol = 8, byrow = TRUE)
colnames(sum_stats) = c("D","p-value","pred_med","actual_med","pred_mean","actual_mean","pred_sd","actual_sd")
rownames(sum_stats) = c("day7","day14","day21","day28","day35","day42")
(metric_mean = sqrt(((mean(sim_D7)-mean(third_val_D7))^2+
              (mean(sim_D14)-mean(third_val_D14))^2+
              (mean(sim_D21)-mean(third_val_D21))^2+
              (mean(sim_D28)-mean(third_val_D28))^2+
              (mean(sim_D35)-mean(third_val_D35))^2+
              (mean(sim_D42)-mean(third_val_D42))^2)/6))
(metric_med = sqrt(((median(sim_D7)-median(third_val_D7))^2+
                      (median(sim_D14)-median(third_val_D14))^2+
                      (median(sim_D21)-median(third_val_D21))^2+
                      (median(sim_D28)-median(third_val_D28))^2+
                      (median(sim_D35)-median(third_val_D35))^2+
                      (median(sim_D42)-median(third_val_D42))^2)/6))
(metric_sd = sqrt(((sd(sim_D7)-sd(third_val_D7))^2+
                      (sd(sim_D14)-sd(third_val_D14))^2+
                      (sd(sim_D21)-sd(third_val_D21))^2+
                      (sd(sim_D28)-sd(third_val_D28))^2+
                      (sd(sim_D35)-sd(third_val_D35))^2+
                      (sd(sim_D42)-sd(third_val_D42))^2)/6))
(metric_D = mean(sum_stats[,"D"]))

############################# temperature specific RMSE##################################
## T4
(Val3_T4.RMSE = sqrt((mean(sim_T4D7) - mean(third_val_T4D7))^2+
                       (mean(sim_T4D14) - mean(third_val_T4D14))^2+
                       (mean(sim_T4D21) - mean(third_val_T4D21))^2+
                       (mean(sim_T4D42) - mean(third_val_T4D42))^2)/5)
(Val3_T4.21.RMSE = sqrt((mean(sim_T4D7) - mean(third_val_T4D7))^2+
                       (mean(sim_T4D14) - mean(third_val_T4D14))^2+
                       (mean(sim_T4D21) - mean(third_val_T4D21))^2)/3)
## T6
(Val3_T6.RMSE = sqrt((mean(sim_T6D7) - mean(third_val_T6D7))^2+
                       (mean(sim_T6D14) - mean(third_val_T6D14))^2+
                       (mean(sim_T6D21) - mean(third_val_T6D21))^2+
                       (mean(sim_T6D28) - mean(third_val_T6D28))^2+
                       (mean(sim_T6D35) - mean(third_val_T6D35))^2)/5)
(Val3_T6.21.RMSE = sqrt((mean(sim_T6D7) - mean(third_val_T6D7))^2+
                          (mean(sim_T6D14) - mean(third_val_T6D14))^2+
                          (mean(sim_T6D21) - mean(third_val_T6D21))^2)/3)

## SC
(Val3_SC.RMSE = sqrt((mean(sim_SCD7) - mean(third_val_SCD7))^2+
                       (mean(sim_SCD14) - mean(third_val_SCD14))^2+
                       (mean(sim_SCD21) - mean(third_val_SCD21))^2+
                       (mean(sim_SCD28) - mean(third_val_SCD28))^2+
                       (mean(sim_SCD35) - mean(third_val_SCD35))^2+
                       (mean(sim_SCD42) - mean(third_val_SCD42))^2)/6)
(Val3_SC.21.RMSE = sqrt((mean(sim_SCD7) - mean(third_val_SCD7))^2+
                          (mean(sim_SCD14) - mean(third_val_SCD14))^2+
                          (mean(sim_SCD21) - mean(third_val_SCD21))^2)/3)

#########################################################################
T4_num_list = sapply(list(sim_T4D7, sim_T4D14, sim_T4D21, sim_T4D35, sim_T4D42, third_val_T4D7, third_val_T4D14, third_val_T4D21, third_val_T4D35,third_val_T4D42), length)
T4_lcm = mLCM(T4_num_list)
T6_num_list = sapply(list(sim_T6D7, sim_T6D14, sim_T6D21, sim_T6D28, sim_T6D35, third_val_T6D7, third_val_T6D14, third_val_T6D21, third_val_T6D28,third_val_T6D35), length)
T6_lcm = mLCM(T6_num_list)
SC_num_list = sapply(list(sim_SCD7, sim_SCD14, sim_SCD21, sim_SCD28, sim_SCD35, sim_SCD42, third_val_SCD7, third_val_SCD14, third_val_SCD21, third_val_SCD28,third_val_SCD35,third_val_SCD42), length)
SC_lcm = mLCM(SC_num_list)

auto_rep = function(x, l){
  rep_num = l / length(x)
  x = rep(x, rep_num)
}


## Validation at 4C
sim_D7 = rbind(auto_rep(sim_T4D7, T4_lcm), auto_rep(sim_T6D7, T6_lcm), auto_rep(sim_SCD7, SC_lcm))
sim_D14 = rbind(auto_rep(sim_T4D14, T4_lcm), auto_rep(sim_T6D14, T6_lcm), auto_rep(sim_SCD14, SC_lcm))
sim_D21 = rbind(auto_rep(sim_T4D21, T4_lcm), auto_rep(sim_T6D21, T6_lcm), auto_rep(sim_SCD21, SC_lcm))
sim_D28 = rbind(auto_rep(sim_T6D28, T6_lcm), auto_rep(sim_SCD28, SC_lcm))
sim_D35 = rbind(auto_rep(sim_T4D35, T4_lcm), auto_rep(sim_T6D35, T6_lcm), auto_rep(sim_SCD35, SC_lcm))
sim_D42 = rbind(auto_rep(sim_T4D42, T4_lcm), auto_rep(sim_SCD42, SC_lcm))




third_val_T4_data = cbind(auto_rep(third_val_T4D7, T4_lcm),auto_rep(sim_T4D7, T4_lcm), 
                        auto_rep(third_val_T4D14, T4_lcm),auto_rep(sim_T4D14, T4_lcm),
                        auto_rep(third_val_T4D21, T4_lcm),auto_rep(sim_T4D21, T4_lcm),
                        auto_rep(third_val_T4D35, T4_lcm),auto_rep(sim_T4D35, T4_lcm),
                        auto_rep(third_val_T4D42, T4_lcm),auto_rep(sim_T4D42, T4_lcm))

third_val_T4_melt = melt(third_val_T4_data)
third_val_T4_melt$Var1 = c(rep(7, 96000),rep(14, 96000), rep(21, 96000), rep(35,96000), rep(42, 96000))

ggplot(data=third_val_T4_melt, aes(x=Var1, y=value,group=Var2))+
  geom_boxplot(fill=c("white","grey","white","grey","white","grey","white","grey",
                    "white","grey"))+
  theme_classic()+
  scale_x_continuous(breaks = c(7,14,21,28,35,42))+
  geom_hline(yintercept = 4.3, linetype=2)+
  labs(y="Spore-former concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (4°C) storage (Days)")+
       #title = "Validation for milks stored at 4°C")+
  theme(axis.title.y = element_text(size=12),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5))



## Validation at 6C
third_val_T6_data = cbind( auto_rep(third_val_T6D7, T6_lcm),auto_rep(sim_T6D7, T6_lcm),
                          auto_rep(third_val_T6D14, T6_lcm),auto_rep(sim_T6D14, T6_lcm), 
                          auto_rep(third_val_T6D21, T6_lcm), auto_rep(sim_T6D21, T6_lcm),
                           auto_rep(third_val_T6D28, T6_lcm),auto_rep(sim_T6D28, T6_lcm),
                          auto_rep(third_val_T6D35, T6_lcm),auto_rep(sim_T6D35, T6_lcm))

third_val_T6_melt = melt(third_val_T6_data)
third_val_T6_melt$Var1 = c(rep(7, 96000),rep(14, 96000), rep(21, 96000), rep(28, 96000),rep(35, 96000))

ggplot(data=third_val_T6_melt, aes(x=Var1,y=value,group=Var2))+
  geom_boxplot(fill=c("white","grey","white","grey","white","grey","white","grey",
                      "white","grey"))+
  theme_classic()+
  geom_hline(yintercept = 4.3, linetype=2)+
  scale_x_continuous(breaks = c(7,14,21,28,35,42))+
  labs(y="Spore-former concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (6°C) storage (Days)")+
  theme(axis.title.y = element_text(size=12),
        plot.title = element_text(hjust = 0.5))



## Validation at SC
third_val_SC_data = cbind( auto_rep(third_val_SCD7, SC_lcm),auto_rep(sim_SCD7, SC_lcm),
                           auto_rep(third_val_SCD14, SC_lcm),auto_rep(sim_SCD14, SC_lcm),
                           auto_rep(third_val_SCD21, SC_lcm),auto_rep(sim_SCD21, SC_lcm),
                           auto_rep(third_val_SCD28, SC_lcm),auto_rep(sim_SCD28, SC_lcm),
                          auto_rep(third_val_SCD35, SC_lcm),auto_rep(sim_SCD35, SC_lcm), 
                           auto_rep(third_val_SCD42, SC_lcm),auto_rep(sim_SCD42, SC_lcm))

third_val_SC_melt = melt(third_val_SC_data)
third_val_SC_melt$Var1 = c(rep(7, 96000),rep(14, 96000), rep(21, 96000), rep(28, 96000),rep(35, 96000),rep(42, 96000))

ggplot(data=third_val_SC_melt, aes(x=Var1,y=value,group=Var2))+
  geom_boxplot(fill=c("white","grey","white","grey","white","grey","white","grey",
                      "white","grey","white","grey"))+
  theme_classic()+
  geom_hline(yintercept = 4.3, linetype=2)+
  scale_x_continuous(breaks = c(7,14,21,28,35,42))+
  labs(y="Spore-former concentrations (log10 cfu/mL)",
       x="Duration of refrigerated (supply chain temperature) storage (Days)")
  theme(axis.title.y = element_text(size=12),
        plot.title = element_text(hjust = 0.5))
