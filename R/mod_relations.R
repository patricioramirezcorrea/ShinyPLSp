# ==============================================================================
# MODULE: STRUCTURAL RELATIONS
# ==============================================================================

mod_relations_ui <- function(id) {
  ns <- NS(id)
  
  nav_panel("Structural relations",
            layout_columns(
              col_widths = c(4, 8),
              card(
                full_screen = TRUE, 
                card_header("Structural model definition"),
                card_body(
                  uiOutput(ns("ui_dep_var")), 
                  uiOutput(ns("ui_indep_vars")),
                  div(class = "action-bar", 
                      actionButton(ns("btn_save_relations"), "Save relations", class = "btn-success btn-compact"), 
                      actionButton(ns("btn_delete_relation"), "Delete selected relation", class = "btn-outline-danger btn-compact")
                  ),
                  tags$hr(), 
                  div(class = "small-muted", "Reciprocal paths and any structural cycles are blocked."), 
                  tags$hr(),
                  h6("Defined relations"), 
                  DTOutput(ns("relation_table"))
                )
              ),
              card(
                full_screen = TRUE, 
                card_header("Model view"),
                card_body(
                  checkboxInput(ns("show_items_model"), "Show indicators in diagram", value = FALSE),
                  DiagrammeR::grVizOutput(ns("model_plot"), height = "520px"), 
                  tags$hr(),
                  h6("Model syntax"), 
                  verbatimTextOutput(ns("model_syntax"))
                )
              )
            )
  )
}

mod_relations_server <- function(id, constructs_rv, relations_rv, result_rv) {
  moduleServer(id, function(input, output, session) {
    
    # Estado interno del módulo
    relation_selected_idx <- reactiveVal(NULL)
    proxy_relations <- dataTableProxy("relation_table", session = session)
    
    # ----------------------------------------------------------------------------
    # REACTIVE: MODEL SYNTAX
    # ----------------------------------------------------------------------------
    model_lavaan <- reactive({
      get_model_syntax(constructs_rv(), relations_rv())
    })
    
    # ----------------------------------------------------------------------------
    # UI DINÁMICA
    # ----------------------------------------------------------------------------
    output$ui_dep_var <- renderUI({
      clist <- constructs_rv()
      if (length(clist) == 0) {
        helpText("Define at least one construct first.")
      } else {
        # session$ns es vital aquí para que el input sea detectado por este módulo
        selectInput(session$ns("dep_var"), "Dependent variable", choices = names(clist))
      }
    })
    
    output$ui_indep_vars <- renderUI({
      clist <- constructs_rv()
      if (length(clist) == 0) {
        helpText("Define at least one construct first.")
      } else {
        allc <- names(clist)
        dep <- if (!is.null(input$dep_var)) input$dep_var else allc[1]
        selectInput(session$ns("indep_vars"), "Independent variables", choices = setdiff(allc, dep), selected = character(0), multiple = TRUE)
      }
    })
    
    observeEvent(input$dep_var, {
      clist <- constructs_rv()
      if (length(clist) == 0) return(NULL)
      indep_choices <- setdiff(names(clist), input$dep_var)
      current_sel <- input$indep_vars[input$indep_vars %in% indep_choices]
      updateSelectInput(session, "indep_vars", choices = indep_choices, selected = current_sel)
    })
    
    # ----------------------------------------------------------------------------
    # LÓGICA DE GUARDADO Y ELIMINACIÓN DE RELACIONES
    # ----------------------------------------------------------------------------
    observeEvent(input$btn_save_relations, {
      req(input$dep_var, input$indep_vars)
      dep <- input$dep_var
      indep <- unique(input$indep_vars)
      
      if (length(indep) == 0) { showNotification("Select at least one independent variable.", type = "error"); return(NULL) }
      if (dep %in% indep) { showNotification("A variable cannot explain itself.", type = "error"); return(NULL) }
      
      df <- relations_rv()
      added <- 0
      for (orig in indep) {
        if (edge_exists(df, orig, dep)) next
        if (edge_exists(df, dep, orig)) { showNotification(paste0("Reciprocal relation blocked: ", orig, " -> ", dep), type = "error", duration = 6); next }
        if (would_create_cycle(df, orig, dep)) { showNotification(paste0("Relation blocked (creates cycle): ", orig, " -> ", dep), type = "error", duration = 6); next }
        df <- rbind(df, data.frame(origin = orig, target = dep, stringsAsFactors = FALSE))
        added <- added + 1
      }
      
      if (added == 0) { showNotification("No new admissible relations were added.", type = "message"); return(NULL) }
      relations_rv(df)
      result_rv(NULL)
      relation_selected_idx(NULL)
      selectRows(proxy_relations, NULL)
      updateSelectInput(session, "indep_vars", selected = character(0))
    })
    
    observeEvent(input$relation_table_rows_selected, { 
      relation_selected_idx(input$relation_table_rows_selected) 
    })
    
    observeEvent(input$btn_delete_relation, {
      idx <- input$relation_table_rows_selected
      idx0 <- relation_selected_idx()
      if (is.null(idx) && is.null(idx0)) { showNotification("Select a relation to delete.", type = "error"); return(NULL) }
      df <- relations_rv()
      delidx <- if (!is.null(idx0)) idx0 else idx
      if (nrow(df) >= delidx) {
        df <- df[-delidx, , drop = FALSE]
        relations_rv(df)
        result_rv(NULL)
        relation_selected_idx(NULL)
        selectRows(proxy_relations, NULL)
      }
    })
    
    # ----------------------------------------------------------------------------
    # OUTPUTS GRÁFICOS Y TABLAS
    # ----------------------------------------------------------------------------
    output$relation_table <- renderDT({ 
      df <- relations_rv()
      if (nrow(df) == 0) empty_dt("No structural relations defined yet.") else standard_dt(df) 
    })
    
    output$model_syntax <- renderText({ 
      model_lavaan() 
    })
    
    output$model_plot <- DiagrammeR::renderGrViz({ 
      plot_graph(clist = constructs_rv(), df_rel = relations_rv(), res = NULL, show_items = isTRUE(input$show_items_model), estimated = FALSE) 
    })
    
    # Retornamos el reactivo model_lavaan para que el global server pueda usarlo al correr cSEM
    return(
      list(
        model_lavaan = model_lavaan
      )
    )
  })
}