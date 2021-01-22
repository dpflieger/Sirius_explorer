dashboardPage(
    title = "Sirius JSON explorer",
    dashboardHeader(title = "Sirius JSON explorer",
                    dropdownMenuOutput("notifications")
                    ), # dashboardHeader
    dashboardSidebar(
        sidebarMenu(
            shinyDirButton(id = "dir_select_input", label = "Select a sirius folder", title = "Select your sirius folder !"),
            menuItem("Search", tabName = "search"),
            numericInput(inputId = "mz", "mz (a numeric)", value = NULL, step = 0.001),
            numericInput(inputId = "mz_approximation", "+/- mz approximation", value = 0.003, step = 0.001)
        ) # sideBarMenu
    ), # dashboardSidebar
    
    # Application title
    dashboardBody(
        tags$style("html, body {overflow: visible !important;"),
        withMathJax("$$\\require{mhchem}$$"),
        tabItems(
            tabItem(tabName = "search", 
                    fluidPage(
                        fluidRow(
                            box(DT::dataTableOutput("precursor.dt.display"), title = "Formula identification", collapsible = TRUE, width = 12, status = "primary")
                            ), # fluidRow
                        fluidRow(
                            box(DT::dataTableOutput("fragments.dt.display"), title = "Fragments", collapsible = TRUE, width = 12, status = "primary")
                            ), # fluidRow
                        fluidRow(
                            box(DT::dataTableOutput("losses.dt.display"), title = "Losses", collapsible = TRUE, width = 12, status = "primary")
                        ) # fluidRow
                    ) # fluidPage
            ) # tabItem
        ) # tabItems
    ) # dashboardBody
) # dashboardPage