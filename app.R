library(shiny)

ui = fluidPage(
  
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
           h3("Number needed to harm (NNH)", align = "center")
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
}

shinyApp(ui = ui, server = server)
