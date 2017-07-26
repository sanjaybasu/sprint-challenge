library(shiny)

ui = navbarPage("Intensive BP Rx",
                
                tabPanel("Risk Calculator",
                         
                         fluidPage(
                           
                           titlePanel("Risk Calculator for Benefit and Harm from Intensive Blood Pressure Treatment"),
                           
                           fluidRow(
                             column(3,
                                    numericInput("age", label = "Age (years)", value = 65),
                                    selectInput("gender", label = "Gender", choices = list(Male=0, Female=1)),
                                    radioButtons("black", label = "Black?", choices = list(No=0, Yes=1), inline = TRUE),
                                    radioButtons("hisp", label = "Hispanic?", choices = list(No=0, Yes=1), inline = TRUE)
                             ),
                             column(3,
                                    numericInput("sysbp", label = "Systolic blood pressure (mm Hg)", value = 140),
                                    numericInput("diabp", label = "Diastolic blood pressure (mm Hg)", value = 90),
                                    numericInput("bptreat", label = "Number of blood pressure medications (0 or more)", value = 0),
                                    radioButtons("aspirin", label = "Taking daily aspirin?", choices = list(No=0, Yes=1), inline = TRUE)
                             ),
                             column(3,
                                    numericInput("totchol", label = "Total cholesterol (mg/dL)", value = 190),
                                    numericInput("hdlchol", label = "HDL cholesterol (mg/dL)", value = 50),
                                    numericInput("trig", label = "Triglycerides (mg/dL)", value = 120),
                                    radioButtons("statin", label = "On statin?", choices = list(No=0, Yes=1), inline = TRUE)
                             ),
                             column(3,
                                    numericInput("sercreat", label = "Serum creatinine (mg/dL)", value = 1.1),
                                    numericInput("bmi", label = "Body Mass index (kg/m^2)", value = 30),
                                    radioButtons("diabetes", label = "Type II diabetes mellitus?", choices = list(No=0, Yes=1), inline = TRUE),
                                    radioButtons("cursmoke", label = "Current smoker?", choices = list(No=0, Yes=1), inline = TRUE),
                                    radioButtons("forsmoke", label = "Former smoker?", choices = list(No=0, Yes=1), inline = TRUE)
                             )
                           ),

                           hr(),

                           fluidRow(
                             column(3,
                                    h3("Absolute risk reduction (ARR):", align = "center")
                             ),
                             column(3,
                                    h3(textOutput("arr"), align = "center")
                             ),
                             column(6,
                                    "Absolute reduction in risk of a cardiovascular event over 5 years (including nonfatal myocardial infarction, stroke, acute coronary syndrome not resulting in myocardial infarction, congestive heart failure, or
                                    cardiovascular death)."
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
                                    1 cardiovascular event over 5 years (including nonfatal myocardial infarction, stroke, acute coronary syndrome not resulting in myocardial infarction, congestive heart failure, or
                                    cardiovascular death)."
                             )
                           ),
                           
                           hr(),
                           
                           
                           fluidRow(
                             column(3,
                                    h3("Absolute risk increase (ARI):", align = "center")
                             ),
                             column(3,
                                    h3(textOutput("ari"), align = "center")
                             ),
                             column(6,
                                    "Absolute increase in risk of a serious adverse event over 5 years (including hypotension, syncope, electrolyte abnormalities, bradycardia,
                                    and acute kidney injury or acute renal failure that is fatal/life-threatening,
                                    results in clinically significant/persistent disability, requires/prolongs a hospitalization,
                                    or causes otherwise clinically significant hazard/harm)."
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
                                    1 serious adverse event over 5 years (including hypotension, syncope, electrolyte abnormalities, bradycardia,
                                    and acute kidney injury or acute renal failure that is fatal/life-threatening,
                                    results in clinically significant/persistent disability, requires/prolongs a hospitalization,
                                    or causes otherwise clinically significant hazard/harm)."
                             )
                           ),

                           hr(),
                           "Note: This calculator is intended for informational purposes only, and has not been prospectively
                           evaluated for impact on clinical practice or patient outcomes. Calculations must be re-checked and 
                           should not be used alone to guide patient care, nor should they substitute for clinical judgment. 
                           Contact: Sanjay Basu, basus@stanford.edu"

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

                         ),

                tabPanel("Disclaimers",

                         h5("This website contains clinical tools and data intended for use by healthcare professionals. These tools do not give professional advice; physicians and other healthcare professionals who use these tools or data should exercise their own clinical judgment as to the information they provide. Consumers who use the tools or data do so at their own risk. Individuals with any type of medical condition are specifically cautioned to seek professional medical advice before beginning any sort of health treatment. For medical concerns, including decisions about medications and other treatments, users should always consult their physician or other qualified healthcare professional.

                            Our content developers have carefully tried to create its content to conform to the standards of professional practice that prevailed at the time of development. However, standards and practices in medicine change as new data become available and the individual medical professional should consult a variety of sources.

                            The contents of the Site, such as text, graphics and images are for informational purposes only. We do not recommend or endorse any specific tests, physicians, products, procedures, opinions, or other information that may be mentioned on the Site.

                            While information on this site has been obtained from sources believed to be reliable, neither we nor our content providers warrant the accuracy of the information contained on this site.

                            We do not give medical advice, nor do we provide medical or diagnostic services. Medical information changes rapidly. Neither we nor our content providers guarantee that the content covers all possible uses, directions, precautions, drug interactions, or adverse effects that may be associated with any therapeutic treatments.

                            Your reliance upon information and content obtained by you at or through this site is solely at your own risk. Neither we nor our content providers assume any liability or responsibility for damage or injury (including death) to you, other persons or property arising from any use of any product, information, idea or instruction contained in the content or services provided to you.

                            We cannot and will not be held legally, financially, or medically responsible for decisions made using these calculators, equations, and algorithms, and this Site is for the use of medical professionals only."),

                         br(),

                         h5("This manuscript was prepared using SPRINT_POP and ACCORD research materials obtained from the NHLBI Biologic Specimen and Data Repository Information Coordinating Center and does not necessary reflect the opinions or views of the SPRINT_POP, the ACCORD, or the NHLBI."),


                         br(),
                         h5("Financial support for this calculator and website and its related materials was provided in part by grants from the National Institute On Minority Health And Health Disparities of the National Institutes of Health under Award Numbers DP2MD010478 and U54MD010724; the National Heart, Lung, And Blood Institute of the National Institutes of Health under Award Number K08HL121056; the National Institute of Diabetes, Digestive and Kidney Diseases of The National Institutes of Health under Award Numbers P60DK20572 and K23DK109200; and the Department of Veterans Affairs HSR&D Service under Award Numbers IIR11-088 and CDA13-021. The funding agreement ensured the authors’ independence in designing the calculations, interpreting the data, writing, and publishing the results. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health or the Department of Veterans Affairs, or of any of the authors' affiliated institutions.")


                 )         
                
                )

server = function(input, output) {
  # Access input values with input$*
  # Save output objects to output$*
  # Build objects with render*({ code })
  netben = reactive({
    round(max(0,((1-((0.881*as.numeric(input$diabetes)+0.943*(1-as.numeric(input$diabetes)))^exp(0.06*as.numeric(input$age)+
                                                                                      -0.117*as.numeric(input$gender)+
                                                                                      -0.058*as.numeric(input$black)+
                                                                                      -0.309*as.numeric(input$hisp)+
                                                                                      0.008*as.numeric(input$sysbp)+
                                                                                      0.002*as.numeric(input$diabp)+
                                                                                      0.169*as.numeric(input$bptreat)+
                                                                                      0.761*as.numeric(input$cursmoke)+
                                                                                      0.139*as.numeric(input$forsmoke)+
                                                                                      0.129*as.numeric(input$aspirin)+
                                                                                      0.157*as.numeric(input$statin)+
                                                                                      0.581*as.numeric(input$sercreat)+
                                                                                      0.004*as.numeric(input$totchol)+
                                                                                      -0.013*as.numeric(input$hdlchol)+
                                                                                      0.0004*as.numeric(input$trig)+
                                                                                      0.01*as.numeric(input$bmi)-
                                                                                      (2.109712*as.numeric(input$diabetes)+6.766*(1-as.numeric(input$diabetes))))))-
                           (1-((0.881*as.numeric(input$diabetes)+0.943*(1-as.numeric(input$diabetes)))^exp(0.06*as.numeric(input$age)+
                                                                                      -0.117*as.numeric(input$gender)+
                                                                                      -0.058*as.numeric(input$black)+
                                                                                      -0.309*as.numeric(input$hisp)+
                                                                                      0.008*as.numeric(input$sysbp)+
                                                                                      0.002*as.numeric(input$diabp)+
                                                                                      0.169*as.numeric(input$bptreat)+
                                                                                      0.761*as.numeric(input$cursmoke)+
                                                                                      0.139*as.numeric(input$forsmoke)+
                                                                                      0.129*as.numeric(input$aspirin)+
                                                                                      0.157*as.numeric(input$statin)+
                                                                                      0.581*as.numeric(input$sercreat)+
                                                                                      0.004*as.numeric(input$totchol)+
                                                                                      -0.013*as.numeric(input$hdlchol)+
                                                                                      0.0004*as.numeric(input$trig)+
                                                                                      0.01*as.numeric(input$bmi)+
                                                                                        1.4+
                                                                                        -0.012*as.numeric(input$age)+
                                                                                        -0.098*as.numeric(input$black)+
                                                                                        -0.009*as.numeric(input$diabp)+
                                                                                        0.207*as.numeric(input$cursmoke)+
                                                                                        -0.0009*as.numeric(input$hdlchol)+
                                                                                        -0.001*as.numeric(input$trig)-
                                                                                      (2.109712*as.numeric(input$diabetes)+6.766*(1-as.numeric(input$diabetes)))))))),digits=4)
    
    
    })
    
    nntcalc = reactive({
      if (netben() <= 0) {nntcalc = "No predicted benefit"}
      else {nntcalc = round(1/netben(),digits=0)}
      nntcalc
    })
  


  output$arr = renderText({ netben() })
  

  output$nnt = renderText({ nntcalc() })
  
  
  netharm = reactive({
    round(max(0,((1-((0.887*as.numeric(input$diabetes)+0.897*(1-as.numeric(input$diabetes)))^exp(0.0330300*as.numeric(input$age)+
                                                                                                   0.1439*as.numeric(input$gender)+
                                                                                                   -0.5452*as.numeric(input$hisp)+
                                                                                                   0.0100800*as.numeric(input$sysbp)+
                                                                                                   -0.0078490*as.numeric(input$diabp)+
                                                                                                   0.1817000*as.numeric(input$bptreat)+
                                                                                                   0.4838000*as.numeric(input$cursmoke)+
                                                                                                   0.0912500*as.numeric(input$forsmoke)+
                                                                                                   0.0472900*as.numeric(input$aspirin)+
                                                                                                   -0.1355000*as.numeric(input$statin)+
                                                                                                   0.7796000*as.numeric(input$sercreat)+
                                                                                                   -0.0055050*as.numeric(input$totchol)+
                                                                                                   0.0080720*as.numeric(input$hdlchol)+
                                                                                                   0.0000717*as.numeric(input$trig)+
                                                                                                   -0.8032000+
                                                                                                   -0.0166300*as.numeric(input$gender)+
                                                                                                   0.0943600*as.numeric(input$cursmoke)+
                                                                                                   0.2864000*as.numeric(input$statin)+
                                                                                                   0.0371600*as.numeric(input$sercreat)+
                                                                                                   0.0044300*as.numeric(input$totchol)+
                                                                                                   0.0009118*as.numeric(input$trig)-
                                                                                                   (3.979551*as.numeric(input$diabetes)+4.343*(1-as.numeric(input$diabetes))))))-
                   (1-((0.887*as.numeric(input$diabetes)+0.897*(1-as.numeric(input$diabetes)))^exp(0.0330300*as.numeric(input$age)+
                                                                                                     0.1439*as.numeric(input$gender)+
                                                                                                     -0.5452*as.numeric(input$hisp)+
                                                                                                     0.0100800*as.numeric(input$sysbp)+
                                                                                                     -0.0078490*as.numeric(input$diabp)+
                                                                                                     0.1817000*as.numeric(input$bptreat)+
                                                                                                     0.4838000*as.numeric(input$cursmoke)+
                                                                                                     0.0912500*as.numeric(input$forsmoke)+
                                                                                                     0.0472900*as.numeric(input$aspirin)+
                                                                                                     -0.1355000*as.numeric(input$statin)+
                                                                                                     0.7796000*as.numeric(input$sercreat)+
                                                                                                     -0.0055050*as.numeric(input$totchol)+
                                                                                                     0.0080720*as.numeric(input$hdlchol)+
                                                                                                     0.0000717*as.numeric(input$trig)-
                                                                                                     (3.979551*as.numeric(input$diabetes)+4.343*(1-as.numeric(input$diabetes)))))))),digits=4)
    
    
  })
  
  nnhcalc = reactive({
    if (netharm() <= 0) {nnhcalc = "No predicted benefit"}
    else {nnhcalc = round(1/netharm(),digits=0)}
    nnhcalc
  })
  
  
  
  output$ari = renderText({ netharm() })
  
  
  output$nnh = renderText({ nnhcalc() })
  
  
  
  
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
  

        }

shinyApp(ui = ui, server = server)
