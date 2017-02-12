# This code implements an automated GUI score calculator using the Shiny package in R.
#
# INSTRUCTIONS:
# (1) install the shiny package from CRAN
# (2) ensure the survival curve files from the following website are downloaded to the working directory for the survival curve tab to work: https://stanfordmedicine.box.com/s/y0dd7j64v4irwfzfbtyjslnok5qvzz9o

library(shiny)

ui = navbarPage("SPRINT Challenge",
  
  tabPanel("Risk Calculator",
  
    fluidPage(
      
      titlePanel("Clinical prediction score for benefits/risks of intensive blood pressure treatment"),
      
      fluidRow(
        column(3,
          numericInput("age", label = "Age (years)", value = 60),
          selectInput("gender", label = "Gender", choices = list(Male=0, Female=1)),
          radioButtons("black", label = "Black?", choices = list(No=0, Yes=1), inline = TRUE),
          radioButtons("hisp", label = "Hispanic?", choices = list(No=0, Yes=1), inline = TRUE)
        ),
        column(3,
          numericInput("sysbp", label = "Systolic blood pressure (mm Hg)", value = 140),
          numericInput("totchol", label = "Total cholesterol (mg/dL)", value = 190),
          numericInput("hdlchol", label = "HDL cholesterol (mg/dL)", value = 50)
        ),
        column(3,
          numericInput("bptreat", label = "Number of blood pressure medications (0 or more)", value = 0),
          radioButtons("cursmoke", label = "Currently smoking tobacco?", choices = list(No=0, Yes=1), inline = TRUE),
          radioButtons("aspirin", label = "Taking daily aspirin?", choices = list(No=0, Yes=1), inline = TRUE),
          radioButtons("statin", label = "On statin?", choices = list(No=0, Yes=1), inline = TRUE)
        ),
        column(3,
          radioButtons("diabetes", label = "Type II diabetes mellitus?", choices = list(No=0, Yes=1), inline = TRUE),
          numericInput("sercreat", label = "Serum creatinine (mg/dL)", value = 1.1),
          numericInput("uralbcreat", label = "Urine albumin/creatinine ratio (mg/g)", value = 43)
        )
      ),
      
      hr(),
      
      fluidRow(
        column(3,
               h2("Score:", align = "center"),
               h4("(0 = most risk, 3 = most benefit)", align = "center")
        ),
        column(3,
               h2(textOutput("score"), align = "center")
        ),
        column(6,
               textOutput("score_interp")
        )
      ),
      
      hr(),
      
      fluidRow(
        column(3,
               h3("Number needed to treat (NNT):", align = "center")
        ),
        column(3,
               h3(textOutput("nnt"), align = "center")
        ),
        column(6,
               "Number of patients with this score necessary to treat with intensive treatment to prevent 
                1 cardiovascular event over 5 years (including nonfatal myocardial infarction, stroke, or 
                cardiovascular death)."
        )
      ),
      
      hr(),
      
      fluidRow(
        column(3,
               h3("Number needed to harm (NNH):", align = "center")
        ),
        column(3,
               h3(textOutput("nnh"), align = "center")
        ),
        column(6,
               "Number of patients with this score needed to treat with intensive treatment to cause 
                1 serious adverse event over 5 years (including hypotension, syncope, electrolyte abnormalities, 
                and acute kidney injury or acute renal failure that is fatal/life-threatening, 
                results in clinically significant/persistent disability, requires/prolongs a hospitalization, 
                or causes otherwise clinically significant hazard/harm)."
        )
      ),
      
      hr(),
      "Note: This calculator is intended for informational purposes only, and has not been prospectively 
       evaluated for impact on clinical practice or patient outcomes. Contact: Sanjay Basu, basus@stanford.edu"
    
    )
    
  ),
  
  tabPanel("Survival Curves",
   
   tags$head(tags$style(
     type="text/css",
     "#surv_cvd img {max-width: 100%; width: 100%; height: auto}"
   )),
   
   tags$head(tags$style(
     type="text/css",
     "#surv_sae img {max-width: 100%; width: 100%; height: auto}"
   )),  
                   
   fluidRow(
     column(2),
     column(4,
        h3(textOutput("risk_score"), align = "center")
     ),
     column(4,
        h3(textOutput("diabetes"), align = "center")
     ),
     column(2)
   ),
   
   hr(),
   
   fluidRow(
     column(6,
        imageOutput("surv_cvd")  
     ),
     column(6,
        imageOutput("surv_sae")  
     )
   )
           
  ),
  
  tabPanel("Summary Statistics",
           
     h4("Risk model was derived from SPRINT trial data and validated against both SPRINT and ACCORD-BP 
        trial data. Summary statistics for both trials are presented below:"),
     
     br(),
     
     fluidRow(
       column(12,
              dataTableOutput('summary')
       )
     )
     
  )

)

server = function(input, output) {
  # Access input values with input$*
  # Save output objects to output$*
  # Build objects with render*({ code })
  pred_score = reactive({0.02268000 + 
                        -0.00042530 * as.numeric(input$age) + 
                        -0.05566000 * as.numeric(input$gender) + 
                         0.01166000 * as.numeric(input$black) + 
                         0.01210000 * as.numeric(input$hisp) + 
                         0.00014520 * as.numeric(input$sysbp) + 
                        -0.00395600 * as.numeric(input$bptreat) + 
                        -0.03623000 * as.numeric(input$cursmoke) +
                        -0.00783500 * as.numeric(input$aspirin) + 
                        -0.01896000 * as.numeric(input$statin) + 
                        -0.03999000 * as.numeric(input$sercreat) + 
                         0.00003832 * as.numeric(input$uralbcreat) +
                         0.00023030 * as.numeric(input$totchol) + 
                         0.00002647 * as.numeric(input$hdlchol)})
  
  score = reactive({
    if (pred_score() < -0.041838836) {score = 0}
    else if (pred_score() < -0.020509052) {score = 1}
    else if (pred_score() < 0.001459577) {score = 2}
    else {score = 3}
    score
  })
  
  output$score = renderText({ score() })
  
  output$risk_score = renderText({ paste("Risk score:", score()) })
  
  output$diabetes = renderText({
    if (input$diabetes == 1) {"Diabetes: Yes"}
    else {"Diabetes: No"}
  })
  
  output$score_interp = renderText({
    if (input$diabetes == 0) {
      if (score() == 0) {interp = "Among SPRINT trial participants (patients without diabetes), 
                                   participants with this score had a non-significant 0.4% decrease 
                                   in cardiovascular events and deaths (95% CI: 2.3% decrease to 1.6% increase, 
                                   P = 0.70), and a significant 3.3% increase in serious adverse events 
                                   (0.7% to 6.0% increase, P = 0.01) from intensive treatment over 5 years."}
      
      else if (score() == 1) {interp = "Among SPRINT trial participants (patients without diabetes), 
                                   participants with this score had a non-significant 1.3% decrease 
                                   in cardiovascular events and deaths (95% CI: 3.1% decrease to 0.5% increase,
                                   P = 0.16) and a significant 4.1% increase in serious adverse events 
                                   (1.8% to 6.4% increase, P <0.001) from intensive treatment over 5 years."}
      
      else if (score() == 2) {interp = "Among SPRINT trial participants (patients without diabetes), 
                                   participants with this score had a non-significant 0.3% decrease 
                                   in cardiovascular events and deaths (95% CI: 1.8% decrease to 1.2% increase,
                                   P = 0.67), and a non-significant 1.1% increase in serious adverse events 
                                   (0.8% decrease to 2.9% increase, P = 0.27) from intensive treatment over 5 years."}
      
      else {interp = "Among SPRINT trial participants (patients without diabetes), participants with this score
                                   had a non-significant 0.3% decrease in cardiovascular events and deaths 
                                   (95% CI: 1.8% decrease to 1.2% increase, P = 0.67), and a non-significant 
                                   1.1% increase in serious adverse events (0.8% decrease to 2.9% increase, 
                                   P = 0.27) from intensive treatment over 5 years."}
    } else {
      if (score() == 0) {interp = "Among ACCORD-BP trial participants (patients with type II diabetes), 
                                   participants with this score had a non-significant 0.2% decrease in 
                                   cardiovascular events and deaths (95% CI: 3.6% decrease to 3.1% increase, 
                                   P = 0.89) and a significant 5.6% increase in serious adverse events 
                                   (2.0% to 9.3% increase, P = 0.002) from intensive treatment over 5 years."}
      
      else if (score() == 1) {interp = "Among ACCORD-BP trial participants (patients with type II diabetes), 
                                   participants with this score had a non-significant 0.5% increase in 
                                   cardiovascular events and deaths (95% CI: 2.8% decrease to 3.9% increase, 
                                   P = 0.76), and a significant 5.6% increase in serious adverse events 
                                   (1.8% to 9.4% increase, P = 0.004) from intensive treatment over 5 years."}
      
      else if (score() == 2) {interp = "Among ACCORD-BP trial participants (patients with type II diabetes), 
                                   participants with this score had a non-significant 0.9% decrease in 
                                   cardiovascular events and deaths (95% CI: 4.1% decrease to 2.4% increase,
                                   P = 0.61), and a significant 9.1% increase in serious adverse events 
                                   (5.8% to 12.4% increase, P<0.001) from intensive treatment over 5 years."}
      
      else {interp = "Among ACCORD-BP trial participants (patients with type II diabetes), participants with
                                   this score had a significant 5.7% decrease in cardiovascular events and 
                                   deaths (95% CI: 2.0% to 9.4% decrease, P=0.003), and a significant 7.8% 
                                   increase in serious adverse events (3.9% to 11.6% increase, P<0.001) 
                                   from intensive treatment over 5 years."}
    }
    interp
  })
  
  output$nnt = renderText({
    if (input$diabetes == 0) {
      if (score() == 0) {nnt = 267}
      else if (score() == 1) {nnt = 77}
      else if (score() == 2) {nnt = 303}
      else {nnt = 30}
    } else {
      if (score() == 0) {nnt = 418}
      else if (score() == 1) {nnt = "No benefit"}
      else if (score() == 2) {nnt = 116}
      else {nnt = 18}
    }
    nnt
  })
  
  output$nnh = renderText({
    if (input$diabetes == 0) {
      if (score() == 0) {nnh = 30}
      else if (score() == 1) {nnh = 24}
      else if (score() == 2) {nnh = 95}
      else {nnh = 77}
    } else {
      if (score() == 0) {nnh = 18}
      else if (score() == 1) {nnh = 18}
      else if (score() == 2) {nnh = 11}
      else {nnh = 13}
    }
    nnh
  })
  
  output$summary = renderDataTable({
    rows = c("Age (years)", "Gender", "Black", "Hispanic", "Systolic blood pressure (mm Hg)",
             "Number of current blood pressure medications (0 or more)",
             "Currently smoking tobacco", "Taking daily aspirin", "On statin",
             "Serum creatinine (mg/dL)", "Urine albumin/creatinine ratio (mg/g)",
             "Total cholesterol (mg/dL)", "High-density lipoprotein (HDL) cholesterol (mg/dL)",
             "Type II diabetes mellitus")
    sprint = c("68 (48 to 90)", "35.6% female", "29.9% black", "10.5% Hispanic", 
               "139.7 (72.0 to 231.0)", "1.8 (0 to 6)", "13.3% currently smoking", "50.8% on aspirin",
               "43.3% on statin", "1.1 (0.4 to 4.0)", "42.6 (1.4 to 5000.0)", "190.1 (71.0 to 438.0)",
               "52.9 (17.0 to 161.0)", "0% (none) with type II diabetes")
    accord = c("62 (44 to 79)", "47.7% female", "24.1% black", "7.0% Hispanic",
               "139.2 (90.0 to 202.0)", "1.7 (0 to 6)", "13.2% currently smoking", "51.9% on aspirin",
               "63.9% on statin", "0.9 (0.1 to 3.5)", "10.1 (0.2 to 990.0)", "192.8 (74.0 to 530.0)",
               "46.5 (5.0 to 155.0)", "100% (all) with type II diabetes")
    table = data.frame(rows, sprint, accord)
    colnames(table) = c("", "Mean (range) in SPRINT trial (N = 9,361)", 
                        "Mean (range) in ACCORD-BP trial (N = 4,733)")
    table
  }, options = list(searching = FALSE, paging = FALSE))
  
  output$surv_cvd = renderImage({
    if (input$diabetes == 0) {
      if (score() == 0) {list(src = "score0nodmcvd.png")}
      else if (score() == 1) {list(src = "score1nodmcvd.png")}
      else if (score() == 2) {list(src = "score2nodmcvd.png")}
      else {list(src = "score3nodmcvd.png")}
    } else {
      if (score() == 0) {list(src = "score0wdmcvd.png")}
      else if (score() == 1) {list(src = "score1wdmcvd.png")}
      else if (score() == 2) {list(src = "score2wdmcvd.png")}
      else {list(src = "score3wdmcvd.png")}
    }
  }, deleteFile = FALSE)
  
  output$surv_sae = renderImage({
    if (input$diabetes == 0) {
      if (score() == 0) {list(src = "score0nodmsae.png")}
      else if (score() == 1) {list(src = "score1nodmsae.png")}
      else if (score() == 2) {list(src = "score2nodmsae.png")}
      else {list(src = "score3nodmsae.png")}
    } else {
      if (score() == 0) {list(src = "score0wdmsae.png")}
      else if (score() == 1) {list(src = "score1wdmsae.png")}
      else if (score() == 2) {list(src = "score2wdmsae.png")}
      else {list(src = "score3wdmsae.png")}
    }
  }, deleteFile = FALSE)
}

shinyApp(ui = ui, server = server)

