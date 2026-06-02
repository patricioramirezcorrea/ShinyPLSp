# ==============================================================================
# MODULE: CONSTRUCTS
# ==============================================================================

mod_constructs_ui <- function(id) {
  ns <- NS(id)
  
  nav_panel("Constructs",
            layout_columns(
              col_widths = c(5, 7),
              card(
                full_screen = TRUE, 
                card_header("Latent construct definition"),
                card_body(
                  textInput(ns("construct_name"), "Construct name", value = ""),
                  selectInput(ns("construct_type"), "Measurement model type", choices = c("Reflective / common factor" = "Common factor", "Composite / emergent" = "Composite"), selected = "Common factor"),
                  div(class = "small-muted", "Common factors are translated as '=~' and composites as '<~'. You can change the measurement model type of an existing construct."),
                  layout_columns(
                    col_widths = c(5, 2, 5),
                    div(tags$b("Available indicators"), selectInput(ns("items_available"), NULL, choices = NULL, multiple = TRUE, selectize= FALSE, size = 12, width = "100%")),
                    div(style = "padding-top: 42px;", actionButton(ns("btn_move_right"), ">>", class = "btn-outline-primary btn-compact"), actionButton(ns("btn_move_left"), "<<", class = "btn-outline-secondary btn-compact")),
                    div(tags$b("Assigned indicators"), selectInput(ns("items_selected"), NULL, choices = NULL, multiple = TRUE, selectize= FALSE, size = 12, width = "100%"))
                  ),
                  div(class = "action-bar mt-3", actionButton(ns("btn_save_construct"), "Save / update construct", class = "btn-success btn-compact"), actionButton(ns("btn_delete_construct"), "Delete selected construct", class = "btn-outline-danger btn-compact"))
                )
              ),
              card(
                full_screen = TRUE, 
                card_header("Defined constructs"), 
                card_body(DTOutput(ns("construct_table")))
              )
            )
  )
}

mod_constructs_server <- function(id, analysis_data, constructs_rv, relations_rv, result_rv, items_available_rv, items_selected_rv, construct_selected_idx, relation_selected_idx) {
  moduleServer(id, function(input, output, session) {
    
    # Proxy para limpiar selecciones de tabla local
    proxy_constructs <- dataTableProxy("construct_table", session = session)
    
    # ----------------------------------------------------------------------------
    # SINCRONIZACIÓN DE UI CON EL ESTADO GLOBAL (SOLUCIÓN AL BUG)
    # ----------------------------------------------------------------------------
    observeEvent(items_available_rv(), {
      updateSelectInput(session, "items_available", choices = items_available_rv(), selected = character(0))
    }, ignoreInit = FALSE)
    
    observeEvent(items_selected_rv(), {
      updateSelectInput(session, "items_selected", choices = items_selected_rv(), selected = character(0))
    }, ignoreInit = FALSE)
    
    # ----------------------------------------------------------------------------
    # HELPER LOCAL: REFRESH CONSTRUCT EDITOR
    # ----------------------------------------------------------------------------
    refresh_construct_editor <- function(current_items = character(0), name = "", type = "Common factor", original = NULL) {
      if (is.null(analysis_data())) return(NULL)
      all_cols <- colnames(analysis_data())
      clist <- constructs_rv()
      construct_selected_idx(NULL)
      if (!is.null(original) && !is.null(clist[[original]])) clist[[original]] <- NULL
      
      used_elsewhere <- unique(unlist(lapply(clist, `[[`, "items"), use.names = FALSE))
      available <- sort(setdiff(all_cols, union(used_elsewhere, current_items)))
      
      # Al actualizar estos RVs, los observeEvent de arriba actualizarán la UI
      items_available_rv(available)
      items_selected_rv(current_items)
      
      updateTextInput(session, "construct_name", value = name)
      updateSelectInput(session, "construct_type", selected = type)
    }
    
    # ----------------------------------------------------------------------------
    # EVENTOS DE INTERFAZ REACTIVA
    # ----------------------------------------------------------------------------
    observeEvent(items_selected_rv(), {
      sel <- items_selected_rv()
      current_name <- trimws(input$construct_name)
      if (length(sel) < 2) return(NULL)
      if (nzchar(current_name)) return(NULL)
      prefixes <- sub("[0-9].*$", "", sel)
      if (length(unique(prefixes)) == 1) {
        pref <- unique(prefixes)
        if (nzchar(pref)) updateTextInput(session, "construct_name", value = pref)
      }
    })
    
    observeEvent(input$btn_move_right, {
      sel <- input$items_available
      if (length(sel) == 0) return(NULL)
      new_selected <- unique(c(items_selected_rv(), sel))
      new_available <- setdiff(items_available_rv(), sel)
      
      items_selected_rv(new_selected)
      items_available_rv(sort(new_available))
    })
    
    observeEvent(input$btn_move_left, {
      sel <- input$items_selected
      if (length(sel) == 0) return(NULL)
      new_selected <- setdiff(items_selected_rv(), sel)
      new_available <- sort(unique(c(items_available_rv(), sel)))
      
      items_selected_rv(new_selected)
      items_available_rv(new_available)
    })
    
    observeEvent(input$btn_save_construct, {
      req(analysis_data())
      name <- sanitize_name(input$construct_name)
      items <- items_selected_rv()
      ctype <- input$construct_type
      clist <- constructs_rv()
      
      if (!nzchar(name)) { showNotification("Please specify a construct name.", type = "error"); return(NULL) }
      if (length(items) == 0) { showNotification("Please assign at least one indicator.", type = "error"); return(NULL) }
      
      if (!is.null(clist[[name]])) {
        old_items <- clist[[name]]$items
        old_type <- clist[[name]]$type
        other_constructs <- clist
        other_constructs[[name]] <- NULL
        used_elsewhere <- unique(unlist(lapply(other_constructs, `[[`, "items"), use.names = FALSE))
        if (any(items %in% used_elsewhere)) { showNotification("Indicators assigned to another construct.", type = "error"); return(NULL) }
        clist[[name]] <- list(items = items, type = ctype)
        if (!identical(old_items, items) || old_type != ctype) result_rv(NULL)
      } else {
        used <- unique(unlist(lapply(clist, `[[`, "items"), use.names = FALSE))
        if (any(items %in% used)) { showNotification("Indicators assigned to another construct.", type = "error"); return(NULL) }
        clist[[name]] <- list(items = items, type = ctype)
        result_rv(NULL)
      }
      
      constructs_rv(clist)
      refresh_construct_editor(current_items = character(0), name = "", type = "Common factor")
      construct_selected_idx(NULL)
      selectRows(proxy_constructs, NULL)
    })
    
    observeEvent(input$construct_table_rows_selected, {
      idx <- input$construct_table_rows_selected
      if (is.null(idx) || length(idx) == 0) return(NULL)
      clist <- constructs_rv()
      name <- names(clist)[idx]
      refresh_construct_editor(current_items = clist[[name]]$items, name = name, type = clist[[name]]$type, original = name)
    })
    
    observeEvent(input$btn_delete_construct, {
      idx <- input$construct_table_rows_selected
      clist <- constructs_rv()
      if (is.null(idx) || length(idx) == 0) { showNotification("Select a construct to delete.", type = "error"); return(NULL) }
      name <- names(clist)[idx]
      if (is.null(name) || is.null(clist[[name]])) return(NULL)
      
      clist[[name]] <- NULL
      constructs_rv(clist)
      
      # Cascading delete of structural relations
      df_rel <- relations_rv()
      if (nrow(df_rel) > 0) {
        df_rel <- df_rel[!(df_rel$origin == name | df_rel$target == name), , drop = FALSE]
        relations_rv(df_rel)
      }
      
      result_rv(NULL)
      construct_selected_idx(NULL)
      relation_selected_idx(NULL)
      
      refresh_construct_editor(current_items = character(0), name = "", type = "Common factor")
      selectRows(proxy_constructs, NULL)
    })
    
    # ----------------------------------------------------------------------------
    # RENDER DE TABLA
    # ----------------------------------------------------------------------------
    output$construct_table <- renderDT({
      standard_dt(constructs_to_df(constructs_rv()))
    })
    
  })
}