# Initial validation
# Use validation dataset from previous model (Buehler et al., 2018) to validate the modification.

for (i in 1:(n_sim *n_units)){
  #Find row in growth parameter data that corresponds to allele sample
  allele_index <- which(spore_growth_import$STorAT == data$STorAT[i]) 
  
  #Calculate the log10N count using the new growth parameters
  data$val_count[i] <- log10N_func(14, spore_growth_import$lag[allele_index],spore_growth_import$mumax[allele_index],data$spore_log10MPN_init_mL[i],spore_growth_import$LOG10Nmax[allele_index])
}
pred_counts_14 = data$val_count[data$val_count>1]


for (i in 1:(n_sim *n_units)){
  #Find row in growth parameter data that corresponds to allele sample
  allele_index <- which(spore_growth_import$STorAT == data$STorAT[i]) 
  
  #Calculate the log10N count using the new growth parameters
  data$val_count[i] <- log10N_func(21, spore_growth_import$lag[allele_index],spore_growth_import$mumax[allele_index],data$spore_log10MPN_init_mL[i],spore_growth_import$LOG10Nmax[allele_index])
}
pred_counts_21 = data$val_count[data$val_count>1]

# load validation dataset
VSL_import_14 <- read.csv("InputFiles/ValidationData1_d14.csv")
VSL_import_21 <- read.csv("InputFiles/ValidationData1_d21.csv")
VSL_counts_14 <- VSL_import_14$log_VSL_D14
VSL_counts_21 <- VSL_import_21$log_VSL
ks.test(pred_counts_14, VSL_counts_14)
ks.test(pred_counts_21, VSL_counts_21)

#summary(pred_counts)
#sd(pred_counts)
#summary(VSL_counts)
#sd(VSL_counts)

pred_counts_14 = cbind(rep("Simulated", length(pred_counts_14)),pred_counts_14)
pred_counts_14 = as.data.frame(pred_counts_14)
names(pred_counts_14) = c("Type","counts")
VSL_counts_14 = cbind(rep("Actual", length(VSL_counts_14)), VSL_counts_14)
VSL_counts_14 = as.data.frame(VSL_counts_14)
names(VSL_counts_14) = c("Type","counts")
val_data_14 = rbind(pred_counts_14,VSL_counts_14)
val_data_14$counts = as.numeric(val_data_14$counts)

ggplot(data = val_data_14, aes(x=Type,y=counts))+
  geom_boxplot()+
  coord_flip()+
  theme_bw()+
  theme(axis.title.y = element_blank())

ggplot(data = val_data_14, aes(x=counts))+
  stat_ecdf(aes(linetype=Type))+
  theme_bw()+
  theme(legend.position = c(0.8, 0.2))+
  labs(x=expression(paste(LOG[10],CFU/mL)), y= "Density")
