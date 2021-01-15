library(shiny)
library(shinyFiles)
library(shinythemes)
library(jsonlite)
library(DT)
library(data.table)
library(stringr)

shinyServer(function(input, output, session) {
    
    shinyDirChoose(input, 'dir_select_input', roots=c(home = "~"), session = session, defaultRoot = "home")
    #dirname(file.choose())
    
    sirius_dir <- reactive({
        return(parseDirPath(c(home = "~"), input$dir_select_input))
    })
    
    output$view_sirius_dir <- renderUI({
        req(input$dir_select_input)
        tagList(renderText(sirius_dir()))
    })
    
    #dir = "/home/dpflieger/4_Tissue+Exudates-sirius301120_38"
    #sirius_dir = "/home/dpflieger/sirius"
    # 
    formula.dt <- reactive({
        req(sirius_dir())
        dt <- try(fread(file.path(sirius_dir(), "formula_identifications.tsv")))
        dt[, json_path := file.path(sirius_dir(), id, "trees", paste0(precursorFormula, "_", gsub(" ", "", adduct), ".json"))]
    })
    
    rds_objects <- reactiveValues(
        fragments.dt = NULL, 
        losses.dt = NULL
    )
    
    observeEvent(input$dir_select_input, {
        
        fragments.rds <- file.path(sirius_dir(), "fragments.rds")
        losses.rds <- file.path(sirius_dir(), "losses.rds")
            
        if(isTRUE(file.exists(losses.rds))) {
            #message("Found losses.rds: ", losses.rds())
            rds_objects$losses.dt <- readRDS(losses.rds)
        }
        
        if(isTRUE(file.exists(fragments.rds))) {
            #message("Found fragments.rds: ", fragments.rds())
            rds_objects$fragments.dt <- readRDS(fragments.rds)
        }
        
        
        if( !isTRUE(file.exists(losses.rds)) || !isTRUE(file.exists(fragments.rds)) ) {
            #message("Loading data...")
            x <- setNames(formula.dt()$json_path, formula.dt()$json_path)
            withProgress(message = 'Loading ', value = 0, {
                json_data <- lapply(x, function(json_path) {
                    #Sys.sleep(0.05)
                    #message("Reading ", json_path)
                    incProgress(1/length(x), message = paste0("Loading JSON: ", which(x == json_path), "/", length(x)))
                    if(file.exists(json_path)) {
                        res <- jsonlite::fromJSON(json_path, flatten = TRUE)
                        list(fragments = setDT(res$fragments), losses = setDT(res$losses))
                    }
                })
            })
            
            rds_objects$losses.dt <- rbindlist(sapply(json_data, "[[", "losses"), fill = TRUE, idcol = "JSON")
            
            fragments.dt <- rbindlist(sapply(json_data, "[[", "fragments"), fill = TRUE, idcol = "JSON")
            # Select wanted columns
            fragments.dt <- fragments.dt[, c("JSON", "id", "molecularFormula", "massDeviation", "score", "mz"), drop = FALSE]
            # Using MathJax we can display correctly molecular formula
            fragments.dt[, molecularFormula := paste0("$$\\ce{", molecularFormula, " }$$")]
            fragments.dt[, score := signif(score, 6)]
            fragments.dt[, massDeviation := signif(as.numeric(stringr::str_extract(massDeviation, '\\d.\\d+')), 6)]
            
            rds_objects$fragments.dt <- fragments.dt
            
            saveRDS(object = rds_objects$losses.dt, file = losses.rds)
            saveRDS(object = rds_objects$fragments.dt, file = fragments.rds)
            
        }
    })
 
    # print(fragments.dt)
    output$precursor.dt.display <- renderDT({
        req(formula.dt())
        message("Rendering fragments...")
        datatable(formula.dt(),
                  rownames = FALSE,
                  # colnames= c("Rank" = "rank", 
                  #             "molecularFormula" = "molecularFormula",
                  #             "# of tRNAs" = "N"),
                  filter = list(position = 'top', clear = TRUE),
                  selection = 'single',
                  options = list(
                      dom = "Qrtip",
                      search = list(regex = TRUE, caseInsensitive = FALSE),
                      columnDefs = list(list(className = 'dt-center', targets = c(1))),
                      pageLength = 6, rowCallback=JS("function( settings ) { MathJax.Hub.Queue(['Typeset',MathJax.Hub]);}"))
        )
    })
    
    # MathJax does not work correctly on DT when the data does not stay static, 
    # so atm we need a callback to the lib everytime it refreshes.
    # Ugly but works tho ¯\_(ツ)_/¯
    output$fragments.dt.display <- renderDT({
        req(rds_objects$fragments.dt)
        message("Rendering fragments...")
        datatable(rds_objects$fragments.dt, 
                  rownames = FALSE,
                  # colnames= c("Rank" = "rank", 
                  #             "molecularFormula" = "molecularFormula",
                  #             "# of tRNAs" = "N"),
                  filter = list(position = 'top', clear = TRUE),
                  selection = 'single',
                  options = list(
                      dom = "Qrtip",
                      search = list(regex = TRUE, caseInsensitive = FALSE),
                      columnDefs = list(list(className = 'dt-center', targets = c(1))),
                      pageLength = 6, rowCallback=JS("function( settings ) { MathJax.Hub.Queue(['Typeset',MathJax.Hub]);}"))
        )
    })
    
    output$losses.dt.display <- renderDT({
        req(rds_objects$losses.dt)
        message("Rendering losses...")
        s = input$fragments.dt.display_rows_selected
        if (length(s)) {
            selected_row = rds_objects$fragments.dt[s, , drop = FALSE]
            datatable(rds_objects$losses.dt[source %in% selected_row$id & JSON %in% selected_row$JSON],
                      rownames = FALSE,
                      options = list(
                          dom = "rti"
                      ))
        }
        
    })
})