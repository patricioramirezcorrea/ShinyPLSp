# ==============================================================================
# MODULE: DATA
# ==============================================================================

mod_data_ui <- function(id) {
  ns <- NS(id)
  
  nav_panel("Data",
            layout_columns(
              col_widths = c(3, 3, 3, 3),
              metric_box("Rows", ns("metric_n_rows")),
              metric_box("Columns", ns("metric_n_cols")),
              metric_box("Applied omission code", ns("metric_omit_code")),
              metric_box("Removed / imputed cases", ns("metric_missing_action"))
            ),
            card(
              card_header("Data preview used for estimation"), 
              card_body(tableOutput(ns("data_preview")))
            )
  )
}

mod_data_server <- function(id, analysis_data_aug_rv, analysis_bundle, omission_code, missing_treatment) {
  moduleServer(id, function(input, output, session) {
    
    # ----------------------------------------------------------------------------
    # METRICS LOGIC (Extraída y adaptada de get_metrics original)
    # ----------------------------------------------------------------------------
    output$metric_n_rows <- renderText({
      bun <- analysis_bundle()
      if (is.null(bun) || is.null(bun$data)) 0 else nrow(bun$data)
    })
    
    output$metric_n_cols <- renderText({
      bun <- analysis_bundle()
      if (is.null(bun) || is.null(bun$data)) 0 else ncol(bun$data)
    })
    
    output$metric_omit_code <- renderText({
      val <- omission_code()
      ifelse(nzchar(trimws(val)), trimws(val), "None")
    })
    
    output$metric_missing_action <- renderText({
      bun <- analysis_bundle()
      mt <- missing_treatment()
      dta <- if (is.null(bun)) NULL else bun$data
      
      if (!is.null(bun)) {
        if (mt == "listwise") {
          paste0("removed ", bun$removed)
        } else if (mt %in% c("mean", "median")) {
          paste0("imputed ", bun$imputed, "; removed ", bun$removed)
        } else {
          paste0("pending NAs: ", sum(is.na(dta)))
        }
      } else {
        "0"
      }
    })
    
    # ----------------------------------------------------------------------------
    # DATA PREVIEW LOGIC
    # ----------------------------------------------------------------------------
    output$data_preview <- renderTable({
      req(analysis_data_aug_rv())
      head(analysis_data_aug_rv(), 10)
    })
    
  })
}