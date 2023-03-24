library(shiny)
library(shinydashboard)
library(shinycssloaders)


## Application title
header = dashboardHeader(title = "Milk Spoilage",
                         tags$li(a(href = 'https://foodsafety.foodscience.cornell.edu/',
                                   img(src = 'cornell_reduced_white.svg', height = 42, width = 144),
                                   style = 'padding-top:5px; padding-bottom:0px;'), 
                                 class = "dropdown"),
                         tags$li(a(href = 'https://foodsafety.foodscience.cornell.edu/digital-dairy/',
                                   img(src = 'mqip.jfif', height = 42, width = 93),
                                   style = 'padding-top:5px; padding-bottom:0px;'), 
                                 class = "dropdown"))

  
## Sidebar concent
sidebar = dashboardSidebar(
             sidebarMenu(
              menuItem("Introduction", tabName = "introduction"),
              menuItem("Model", tabName = "model"),
              menuItem("Output", tabName = "output")))
  
# Body content
body =  dashboardBody(
    tags$head(tags$style(HTML('
        /* logo */
        .skin-blue .main-header .logo {
                              background-color: #800000;
                              }
        /* logo when hovered */
        .skin-blue .main-header .logo:hover {
                              background-color: #800000;
                              }

        /* navbar (rest of the header) */
        .skin-blue .main-header .navbar {
                              background-color: #B31B1B;
                              }        

        /* active selected tab in the sidebarmenu */
        .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                              background-color: #F7F7F7;
                              }

        /* other links in the sidebarmenu */
        .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                              color: #B31B1B;
                              font-weight: bold;
                              }

        /* toggle button when hovered  */                    
         .skin-blue .main-header .navbar .sidebar-toggle:hover{
                              background-color: #800000;
                              }
        .skin-blue .sidebar-menu > li.active > a,
        .skin-blue .sidebar-menu > li:hover > a {
                    border-left-color: #B31B1B;
                              }            
         '))),
    tabItems(
        #First tab content
        tabItem(tabName = "introduction",
              box(
                width = 8,
                h3("Model Overview"),
                p("This is a web interface for a Monte-Carlo simulation model that aims to predict the spoilage 
                  of half-gallon fluid milk due to psychrotolerant spore-forming bacteria. This model allows the user to 
                  1) input the concentration of psychrotolerant spore-formers in pasteurized half-gallon milk, 2) define a spoilage threshold 
                  level, and 3) select various intervention strategies. As a result, the model will simulate the 
                  percent of spoiled half-gallon milk containers over 35 days of consumer storage. Click the", strong('Model'), "tab to start!"),
                h3("Developers"),
                p("This web app is developed by researchers from", strong("Milk Quality Improvement Lab (MQIP)") ,"in the 
                  department of food science at", strong("Cornell University"),". Our lab aims to leverage digital tools to 
                  facilitate management decision-making in reducing spoilage of dairy products. For more information 
                  about we do, please visit our",
                    a("lab website",
                      href = "https://foodsafety.foodscience.cornell.edu/"), ("and"), 
                    a("digital dairy central hub", 
                      href = "https://foodsafety.foodscience.cornell.edu/digital-dairy/")),
                h3("Funding"),
                p("Funding for this project was provided by the Foundation for Food and Agriculture Research 
                  (FFAR, Washington, DC; award no. CA18-SS-0000000206"),
                h3("Disclaimers"),
                p("The model for predicting spoilage in fluid milk (the 'Model') is provided 'as is'
                    at no cost to the users. The predicted values provided by the Model are only
                    estimations and not actual values obtained in a real milk supply chain. 
                    Actual values will differ from the Model predictions due to variations 
                    such as milk type, microbial load, pressure and other factors not covered by the Model.
                    By using the Model the user is accepting economic and all other risks associated 
                    with it and remain responsible for compliance with all regulatory, food safety 
                    and other requirements related to their cheesemaking activities."),
                br(),
                p(em("CORNELL UNIVERSITY EXPRESSLY DISCLAIMS ANY AND ALL WARRANTIES IN CONNECTION WITH THE MODEL,
                   INCLUDING BUT NOT LIMITED TO WARRANTIES OF FITNESS FOR A PARTICULAR PURPOSE, WARRANTIES OF
                   MERCHANTABILITY, AND WARRANTIES OF NON-INFRINGEMENT.")),
                br(),
                p("By using the Model you acknowledge and agree to the foregoing and waive any and all claims
                    of any kind and description against Cornell University, its officers, employees and agents, 
                    in connection with your use of the Model and any decisions you may make as a result of using 
                    the Model."),
                h3("Contacts"),
                p(strong("Luke Qian"), "- PhD student. Email:cq87@cornell.edu", 
                  a("LinkedIn", 
                  href = "https://www.linkedin.com/in/luke-qian-751a8813a/")),
                p(strong("Aljosa Trmcic"), "- Dairy Extension Associate. Email:at543@cornell.edu", 
                  a("LinkedIn", 
                    href = "https://www.linkedin.com/in/aljosa-trmcic-b8767720/")),
                p(strong("Martin Wiedmann"), "- Gellert Family Professor in Food Safety. Email:mw16@cornell.edu",
                  a("LinkedIn", 
                  href = "https://www.linkedin.com/in/martin-wiedmann-6a731810/"))
                )),
    
        #Second tab content
        tabItem(tabName = "model",
              fluidRow(
                  box(
                    width = 4,
                    title = "User input",
                    solidHeader = TRUE,
                    collapsible = FALSE,
                  #Style
                  tags$head(tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #B31B1B; border-color: #B31B1B;}
                                    .btn {background: #B31B1B; border-color: #B31B1B;}
                                    .btn:hover{background: #000000; border-color: #000000;}
                                    .btn.active{background: #B31B1B; border-color: #000000;}
                                    .btn.visited{background: #B31B1B; border-color: #B31B1B;}
                                    .btn.btn-primary:focus{background-color: #F7F7F7; color:#000000; border-color: #000000;}
                                    .box-header {color:#ffffff;background:#B31B1B}
                                    .box-header{border-bottom-color:#B31B1B;
                                                                border-left-color:#B31B1B;
                                                                border-right-color:#B31B1B;
                                                                border-top-color:#B31B1B;}"
                    ))),
                      
                  #User input
                  
                      numericInput("count_mean", "What is the average (mean) spore concentration (log10 MPN/mL) in milk?", value = -0.72),
                      numericInput("count_sd", "What is the standard deviation of spore concentration (log10 MPN/mL) in milk?", value = 0.99),
                      selectInput("threshold","Spoilage threshold", c("US regulation limit (Pasteurized Milk Ordinance): 20,000 CFU/mL" = log10(20000), 
                                                                      "Limit for consumer to detect the defects: 1,000,000 CFU/mL" = log10(1000000))),
                      h5("The shelf life is defined as the last day at which a certain percentage (%) of spoiled milk containers exceeding the spoilage threshold."),
                      sliderInput("shelfLife_threshold", "Set this percentage (%)", value = 50, min = 0, max = 100),
                      radioButtons("sporeRed", HTML("Select intervention strategies that reduce bacteria load"),
                                   c("No intervention" = "none",
                                     "Microfiltration (2.2log reduction)" = "mf",
                                     "Bactofugation single-pass (1.4log reduction)" = "bf1",
                                     "Bactofugation double-pass (2log reduction)" = "bf2")),
                      h5(strong("Select intervention strategies that improve the temperature control")),
                      radioButtons("f_intervention", HTML("<em>Facility storage</em>"),
                                         c("No intervention" = "none",
                                           "Set storage temperature from 4 to 3 \u00B0C" = "f_reduceT",
                                           "Implement extreme cooling (reduce temperature from 4 to 1 \u00B0C)" = "f_supercool",
                                           "Improve the cooling system to reduce temperature variability" = "f_reduceVar")),
                      radioButtons("ftr_intervention", HTML("<em>Facility-to-retail transportation</em>"),
                                         c("No intervention" = "none",
                                           "Set a temperature alarm system in delivery truncks" = "ftr_alarm",
                                           "Optimize distribution routes to shorten delivery time" = "ftr_opt")),
                      radioButtons("r_intervention", HTML("<em>Retail storage</em>"),
                                         c("No intervention" = "none",
                                           "Reduce average (mean) storage temperature from 2.3 to 1.8 \u00B0C" = "r_reduceT",
                                         "Set a temperature alarm system to limit temperature below 4 \u00B0C" = "r_alarm",
                                         "Improve refrigeration system to reduce temperature variability" = "r_reduceVar")),
                      submitButton("Submit" ,icon("refresh"))), 
                  
    
                  box(
                    width = 6,
                    title = "Simulated outcomes",
                    solidHeader = TRUE,
                    collapsible = FALSE,
                    height = "100%",
                    #Histogram and predicted percent spoilage
                    withSpinner(plotOutput("plot"), color = "#B31B1B"),
                    h3(textOutput("shelfLife"))
                    )
              )
        ),
        
        tabItem(tabName = "output",
                
                box(
                  width = 4,
                  title = "Download data",
                  solidHeader = TRUE,
                  collapsible = FALSE,
                  p("Click the button below to download the data as a csv file"),
                  downloadButton("downloadData", "Download", class = "downloadButt")
                    ),
                
                # Style
                tags$head(tags$style(HTML(".downloadButt:hover{background:#000000; color:#FFFFFF;border-color: #000000}
                                 .downloadButt{background:#B31B1B; color:#FFFFFF;border-color: #B31B1B}
                                 .box-header {color:#ffffff;background:#B31B1B}
                                 .box-header{border-bottom-color:#B31B1B;
                                                  border-left-color:#B31B1B;
                                                  border-right-color:#B31B1B;
                                                  border-top-color:#B31B1B;}"   )))
          
          
          
          
        )
    )
)
  

dashboardPage(header, sidebar, body)