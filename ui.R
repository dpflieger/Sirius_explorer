fluidPage(
    # theme
    theme = shinytheme("cerulean"),
    withMathJax("$$\\require{mhchem}$$"),
    # Application title
    navbarPage("Sirius JSON explorer",
               # 1st tabpanel ------------------------------------------------------------
               tabPanel("Search",
                        fluidRow(
                            column(3, shinyDirButton(id = "dir_select_input", 
                                                     label = "Select a sirius folder", 
                                                     title = "Select your sirius folder !")),
                            column(4, uiOutput("view_sirius_dir"))
                        ), 
                        fluidRow(
                            column(2, numericInput(inputId = "mz_integer_input", "mz (a number)", value = NULL)), 
                            column(3, numericInput(inputId = "mz_approximation_input", "+/- mz approximiation", value = 0.003, step = 0.001))
                        ), # fluidRow
                        h3("Precursors"),
                        wellPanel(style = "background-color: #fff; border-color: #2c3e50;",
                                  DT::dataTableOutput("precursor.dt.display")
                        ), # wellPanel
                        h3("Fragments"),
                        wellPanel(style = "background-color: #fff; border-color: #2c3e50;",
                                  DT::dataTableOutput("fragments.dt.display")
                        ), # wellPanel
                        h3("Losses"),
                        wellPanel(style = "background-color: #fff; border-color: #2c3e50;",
                                  DT::dataTableOutput("losses.dt.display")
                        ) # wellPanel
               ) # tabPanel Search
               
    ) # navbarPage
) 