shinyServer(function(input, output, session) {
    
    # Get all hardrives names (Windows/linux/MacOS friendly)
    volumes <- getVolumes()() # yeah weird design
    
    # Select a directory (via the ShinyFiles package)
    shinyDirChoose(input, 'dir_select', roots = volumes, session = session)

# Check directory ---------------------------------------------------------
    
    # Reload last session directory
    if(isTRUE(file.exists("latest_dir.rds"))) {
        print("Loading the latest directory")
        tmp <- readRDS("latest_dir.rds")
        if(dir.exists(tmp)) sirius_dir <- reactive({ tmp })
        
    } else {
        # Change input$dir from vector to character
        sirius_dir <- reactive({
            parseDirPath(volumes, input$dir_select)
        })
        
    }
    
    observeEvent(sirius_dir(), {
        saveRDS(sirius_dir(), "latest_dir.rds")
    })
    
    
    # Render selected directory
    output$view_sirius_dir <- renderUI({
        req(sirius_dir())
        HTML(sirius_dir())
    })
    
# Load formula_identifications file ---------------------------------------
    formula.dt <- reactive({
        req(sirius_dir())
        dt <- try(fread(file.path(sirius_dir(), "formula_identifications.tsv")))
        dt[, json_path := file.path(sirius_dir(), id, "trees", paste0(precursorFormula, "_", gsub(" ", "", adduct), ".json"))]
        # Extract ID from the last column
        dt[, c("info", "useless", "FBMM_id") := tstrsplit(id, "_", fixed=TRUE)]
    })
    
    rds_objects <- reactiveValues(
        fragments.dt = NULL, 
        losses.dt = NULL
    )
    

# Load data ---------------------------------------------------------------
    observeEvent(sirius_dir(), {
        
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
            # Setup a progress to see the loading state of the JSON files
            withProgress(message = 'Loading ', value = 0, {
                json_data <- lapply(x, function(json_path) {
                    incProgress(1/length(x), message = paste0("Loading JSON: ", which(x == json_path), "/", length(x)))
                    if(file.exists(json_path)) {
                        res <- jsonlite::fromJSON(json_path, flatten = TRUE)
                        
                        # Create named vector for fragments mz
                        fragments_mz <- res$fragments$mz
                        names(fragments_mz) <- res$fragments$id
                        setDT(res$losses)
                       
                        # Get the mz of source and target from the fragments_mz named vector
                        if(length(fragments_mz) > 1){
                          res$losses[, source_mz := fragments_mz[as.character(source)]]
                          res$losses[, target_mz := fragments_mz[as.character(target)]]
                          # Compute the delta mz
                          res$losses[, delta_mz := source_mz - target_mz]
                        }
                        list(fragments = setDT(res$fragments), losses = setDT(res$losses))
                    }
                })
            })
            
            rds_objects$losses.dt <- rbindlist(sapply(json_data, "[[", "losses"), fill = TRUE, idcol = "JSON")
            
            fragments.dt <- rbindlist(sapply(json_data, "[[", "fragments"), fill = TRUE, idcol = "JSON")
            # Select wanted columns
            fragments.dt <- fragments.dt[, c("JSON", "id", "molecularFormula", "massDeviation", "score", "mz"), drop = FALSE]
            # Using MathJax we can display correctly molecular formula
            #fragments.dt[, molecularFormula := paste0("$$\\ce{", molecularFormula, " }$$")]
            
            fragments.dt[, score := signif(score, 6)]
            fragments.dt[, massDeviation := signif(as.numeric(stringr::str_extract(massDeviation, '\\d.\\d+')), 6)]
            fragments.dt[, c("info", "useless", "FBMM_id") := tstrsplit(basename(dirname(dirname(JSON))), "_", fixed=TRUE)]
            fragments.dt[, useless := NULL]
            
            rds_objects$fragments.dt <- fragments.dt
            
            saveRDS(object = rds_objects$losses.dt, file = losses.rds)
            saveRDS(object = rds_objects$fragments.dt, file = fragments.rds)
        }
    })
 
    # Reactive fragments table
    fragments.dt <- reactive({
        dt <- rds_objects$fragments.dt
        # if user uses filter
        if(!is.na(input$mz))
            dt <- dt[mz <= input$mz + (input$mz * input$ppm_approximation * 1e-06) & mz >= input$mz - (input$mz * input$ppm_approximation * 1e-06)]
        return(dt)
    }) 

# Render formula  ---------------------------------------------------------
    output$precursor.dt.display <- renderDT({
        req(formula.dt())
        message("Rendering formula identifications file...")
        
        dt <- formula.dt()
        
        s <- input$fragments.dt.display_rows_selected
        if (length(s)) {
          selected_row <- fragments.dt()[s, , drop = FALSE]
          dt <- dt[FBMM_id %in% selected_row$FBMM_id, ]
        }

        dt_f <- dt[, c("rank", "molecularFormula", "IsotopeScore", "ionMass", "id"), with = FALSE]
        setcolorder(dt_f, c("id", "rank"))
        
        datatable(dt_f,
                  rownames = FALSE,
                  filter = list(position = 'top', clear = TRUE),
                  selection = 'single',
                  options = list(
                      dom = "Qrtip",
                      search = list(regex = TRUE, caseInsensitive = FALSE),
                      columnDefs = list(list(className = 'dt-center', targets = c(1))),
                      pageLength = 6
                      #rowCallback=JS("function( settings ) { MathJax.Hub.Queue(['Typeset',MathJax.Hub]);}")
                      )
        )
    })

# Render fragments --------------------------------------------------------
    
    # MathJax does not work correctly on DT when the data does not stay static, 
    # so atm we need a callback to the lib everytime it refreshes.
    # Ugly but works tho ¯\_(ツ)_/¯
    output$fragments.dt.display <- renderDT({
        req(fragments.dt())
        message("Rendering fragments...")
        
        dt <- fragments.dt()

        datatable(dt[, c("JSON", "id", "molecularFormula", "massDeviation", "mz", "FBMM_id"), with = FALSE], 
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
                      pageLength = 6
                      #rowCallback=JS("function( settings ) { MathJax.Hub.Queue(['Typeset',MathJax.Hub]);}")
                      )
        )
    })

# Render losses data ------------------------------------------------------
    output$losses.dt.display <- renderDT({
        req(rds_objects$losses.dt)
        message("Rendering losses...")
        # Extract the selected row (number) in the fragments.dt
        s = input$fragments.dt.display_rows_selected
        if (length(s) != 0) {
            # Get the complete row from the table
            selected_row = fragments.dt()[s, , drop = FALSE]
            
            # # Match elements with the other table
            # datatable(rds_objects$losses.dt[JSON %in% selected_row$JSON & source %in% selected_row$id , c("JSON", "source", "target", "molecularFormula", "score", "delta_mz", "source_mz", "target_mz"), with=FALSE],
            #           rownames = FALSE,
            #           options = list(
            #               dom = "rti"
            #           ))            
            # Match elements with the other table
            datatable(rds_objects$losses.dt[JSON %in% selected_row$JSON, c("JSON", "source", "target", "molecularFormula", "score", "delta_mz", "source_mz", "target_mz"), with=FALSE],
                      rownames = FALSE,
                      options = list(
                          dom = "rti"
                      ))
        }
    })
    

# visNetwork --------------------------------------------------------------
# https://datastorm-open.github.io/visNetwork/shiny.html

    output$network <- renderVisNetwork({
      req(rds_objects$losses.dt, input$fragments.dt.display_rows_selected)
      s = input$fragments.dt.display_rows_selected
      # Get the complete row from the table
      selected_row = fragments.dt()[s, , drop = FALSE]

      nodes <- rds_objects$fragments.dt[JSON %in% selected_row$JSON, c("JSON", "id", "molecularFormula")]
      setnames(nodes, c("molecularFormula"), c("label"))
      edges <- rds_objects$losses.dt[JSON %in% selected_row$JSON, c("JSON", "source", "target", "molecularFormula")]
      setnames(edges, c("source", "target", "molecularFormula"), c("from", "to", "label"))
      
      t <- visNetwork(nodes, edges, width = "120") %>% 
        visIgraphLayout(layout = "layout_as_tree", randomSeed = 1337)
      print("Making tree")
      t$x$nodes$y <- -t$x$nodes$y # reverse the tree, makes it TOPDOWN
      t$x$nodes$shape <- "box" # shape to box
      
      return(t)
    })
    
})