Author: Luke Qian, Sarah I. Murphy, Tim Lott

Title: Development of a predictive model for milk spoilage along the supply chain due to psychrotolerant spore-formers 

In this project, a Monte Carlo simulation model was developed to predict the percentage of spoiled milk containers in a supply chain from processing facility to consumer. The sensitivity analysis was performed to identify important model parameters. The what-if scenarios were used to evaluate the effectiveness of various intervention stratgies including spore removal and temperature control. A user-friendly interface was developed in R Shiny for this model, which is deployed via Shinyapp.io and accessible at https://lukeqian.shinyapps.io/MilkSporeModel/.

The main model file is named "MilkSporeModel.RMD"

The R script named "UtilityFunctions.R" includes all the functions needed for primary and secondary growth moodel and will be called inside the main model file.

The model inputs including allelic type frequency, initial spore concentrations, microbial growth characteristics are located in "InputFiles".

The script codes for model validation are located in "Validation" file. These scripts will be called and run inside the main model file. 

The "App" folder contains the server and ui codes for R Shiny app developed for this model. 
