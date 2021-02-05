dashboardPage(
    title = "Sirius JSON explorer",

    dashboardHeader(title = "Sirius JSON explorer",
                    tags$li(class = "dropdown", tags$a((htmlOutput("view_sirius_dir"))))
                    ), # dashboardHeader
    dashboardSidebar(
        sidebarMenu(
            shinyDirButton(id = "dir_select", label = "Select a sirius folder", title = "Select your sirius folder !", ),
            menuItem("Search", tabName = "search"),
            numericInput(inputId = "mz", "Mass-to-charge ratio (m/z)", value = NULL, step = 0.001),
            numericInput(inputId = "ppm_approximation", "Mass accuracy (ppm)", value = 3, step = 1)
        ) # sideBarMenu
    ), # dashboardSidebar
    
    # Application title
    dashboardBody(
        #tags$style("html, body {overflow: visible !important;"),
        tabItems(
            tabItem(tabName = "search", 
                    fluidPage(
                        fluidRow(
                            box(DT::dataTableOutput("precursor.dt.display"), title = "Formula identification", collapsible = TRUE, width = 12, status = "warning")
                            ), # fluidRow
                        fluidRow(
                            box(DT::dataTableOutput("fragments.dt.display"), title = "Fragments", collapsible = TRUE, width = 12, status = "primary")
                            ), # fluidRow
                        fluidRow(
                            box(DT::dataTableOutput("losses.dt.display"), title = "Losses", collapsible = TRUE, width = 12, status = "info")
                        ), # fluidRow
                        fluidRow(
                            box(visNetworkOutput("network"),  title = "Trees", collapsible = TRUE, width = 12)
                            )
                    ) # fluidPage
            ) # tabItem
        ), # tabItems
        withMathJax("$$\\require{mhchem}$$"),
    ) # dashboardBody
) # dashboardPage