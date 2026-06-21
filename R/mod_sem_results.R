# ==============================================================================
# MODULE: SEM RESULTS
# ==============================================================================

mod_sem_results_ui <- function(id) {
  ns <- NS(id)

  nav_panel("SEM results",
            layout_columns(
              col_widths = c(3, 3, 3, 3),
              metric_box("Constructs", ns("metric_constructs")),
              metric_box("Relations", ns("metric_relations")),
              metric_box("Estimated paths", ns("metric_paths")),
              metric_box("Reported R-squared", ns("metric_r2"))
            ),
            navset_card_pill(
              nav_panel("Model fit", DTOutput(ns("model_fit_table"))),
              nav_panel("Construct types", DTOutput(ns("construct_type_table"))),
              nav_panel("R-squared", DTOutput(ns("r2_table"))),
              nav_panel("Paths", DTOutput(ns("path_table"))),
              nav_panel("Loadings", DTOutput(ns("loading_table"))),
              nav_panel("Measurement quality", DTOutput(ns("mm_quality_table"))),
              nav_panel("HTMT", DTOutput(ns("htmt_table"))),
              nav_panel("Fornell-Larcker", DTOutput(ns("fornell_table"))),
              nav_panel("Hypothesis table",         DTOutput(ns("hypothesis_table"))),
              nav_panel("Indirect & Total effects", DTOutput(ns("indirect_table")))
            ),

            card(
              full_screen = TRUE,
              card_header("Estimated model diagram"),
              card_body(
                checkboxInput(ns("show_items_estimated"), "Show indicators in diagram", value = FALSE),
                DiagrammeR::grVizOutput(ns("estimated_model_plot"), height = "520px")
              )
            ),
            card(
              card_header("Detailed results export"),
              card_body(downloadButton(ns("download_excel"), "Export detailed Excel", class = "btn-success btn-compact"))
            ),

            card(
              card_header("Model status"),
              card_body(verbatimTextOutput(ns("summary_text")))
            ),

            card(
              card_header("Full model statistics (raw output)"),
              card_body(
                div(style = "max-height: 400px; overflow-y: auto;",
                    verbatimTextOutput(ns("raw_sem_summary")))
              )
            )
  )
}

mod_sem_results_server <- function(id, result_rv, constructs_rv, relations_rv, analysis_data, analysis_bundle, raw_data_rv, resample_method, n_boot, project_name) {
  moduleServer(id, function(input, output, session) {

    render_dt_with_pval_color <- function(df, p_col_name = "P_value", digits = 3) {
      if (is.null(df) || nrow(df) == 0) return(empty_dt("No data available."))
      if (!is.null(digits)) df <- format_df(df, digits)
      dt <- datatable(df, selection = "single", rownames = FALSE,
                      options = list(dom = "tip", pageLength = 100, scrollX = TRUE))
      p_cols <- names(df)[grepl("p.*value|^p$|P_value|_P_value", names(df), ignore.case = TRUE)]
      for (col in p_cols) {
        col_idx <- which(names(df) == col)
        if (length(col_idx) == 1) {
          dt <- DT::formatStyle(dt, col,
                                backgroundColor = DT::styleInterval(0.05, c("#d4edda", "white")),
                                color           = DT::styleInterval(0.05, c("#155724", "black")))
        }
      }
      dt
    }

    # ----------------------------------------------------------------------------
    # METRICS
    # ----------------------------------------------------------------------------
    output$metric_constructs <- renderText({
      length(constructs_rv())
    })

    output$metric_relations <- renderText({
      nrow(relations_rv())
    })

    output$metric_paths <- renderText({
      res <- result_rv()
      if (is.null(res) || is.null(extract_paths(res))) 0 else nrow(extract_paths(res))
    })

    output$metric_r2 <- renderText({
      res <- result_rv()
      if (is.null(res) || is.null(extract_r2(res))) 0 else nrow(extract_r2(res))
    })

    # ----------------------------------------------------------------------------
    # TABLES
    # ----------------------------------------------------------------------------
    output$model_fit_table <- renderDT({
      res <- result_rv()
      if (is.null(res)) return(empty_dt("The model has not been estimated yet."))
      df <- extract_model_fit(res)
      if (is.null(df) || nrow(df) == 0) return(empty_dt("Model fit indices not available."))

      df_fmt <- format_df(df, 3)
      dt <- datatable(df_fmt, selection = "single", rownames = FALSE,
                      options = list(dom = "tip", pageLength = 100, scrollX = TRUE))


      if ("Value" %in% names(df) && "Index" %in% names(df)) {
        for (i in seq_len(nrow(df))) {
          idx <- df$Index[i]
          val <- df$Value[i]
          if (is.na(val)) next
          color <- if (idx == "SRMR") {
            if (val < 0.08) "#d4edda" else if (val < 0.10) "#fff3cd" else "#f8d7da"
          } else if (idx == "NFI") {
            if (val > 0.90) "#d4edda" else if (val > 0.80) "#fff3cd" else "#f8d7da"
          } else if (idx == "GoF") {
            if (val > 0.36) "#d4edda" else if (val > 0.25) "#fff3cd" else "#f8d7da"
          } else {
            "white"
          }
          txt <- if (color == "#d4edda") "#155724" else if (color == "#fff3cd") "#856404" else if (color == "#f8d7da") "#721c24" else "black"
          dt <- DT::formatStyle(dt, "Value",
                                valueColumns = "Index",
                                backgroundColor = DT::styleEqual(idx, color),
                                color           = DT::styleEqual(idx, txt))
        }
      }
      dt
    })


    output$construct_type_table <- renderDT({ res <- result_rv(); if (is.null(res)) return(empty_dt("The model has not been estimated yet.")); df <- get_construct_types(res); if (is.null(df) || nrow(df) == 0) empty_dt("Construct types not available.") else standard_dt(df) })

    output$r2_table <- renderDT({
      res <- result_rv()
      if (is.null(res)) return(empty_dt("The model has not been estimated yet."))
      df <- extract_r2(res)
      if (is.null(df) || nrow(df) == 0) return(empty_dt("R-squared not available."))

      df_fmt <- format_df(df, 3)
      dt <- datatable(df_fmt, selection = "single", rownames = FALSE,
                      options = list(dom = "tip", pageLength = 100, scrollX = TRUE))


      num_cols <- names(df)[sapply(df, is.numeric)]
      for (col in num_cols) {
        dt <- DT::formatStyle(dt, col,
                              backgroundColor = DT::styleInterval(c(0.13, 0.26),
                                                                  c("#f8d7da", "#fff3cd", "#d4edda")),
                              color           = DT::styleInterval(c(0.13, 0.26),
                                                                  c("#721c24", "#856404", "#155724")))
      }
      dt
    })

    output$path_table <- renderDT({
      res <- result_rv()
      if (is.null(res)) return(empty_dt("The model has not been estimated yet."))
      df <- extract_paths(res)
      if (is.null(df) || nrow(df) == 0) return(empty_dt("Paths not available."))
      render_dt_with_pval_color(df, digits = 3)
    })

    output$loading_table <- renderDT({
      res <- result_rv()
      if (is.null(res)) return(empty_dt("The model has not been estimated yet."))
      df <- extract_loadings(res)
      if (is.null(df) || nrow(df) == 0) return(empty_dt("Loadings not available."))

      df_fmt <- format_df(df, 3)
      dt <- datatable(df_fmt, selection = "single", rownames = FALSE,
                      options = list(dom = "tip", pageLength = 100, scrollX = TRUE))


      est_col <- intersect(names(df), c("Estimate","Loading","loading"))
      if (length(est_col) > 0)
        dt <- DT::formatStyle(dt, est_col[1],
                              backgroundColor = DT::styleInterval(c(0.60, 0.70),
                                                                  c("#f8d7da", "#fff3cd", "#d4edda")),
                              color           = DT::styleInterval(c(0.60, 0.70),
                                                                  c("#721c24", "#856404", "#155724")))


      p_cols <- names(df)[grepl("p.*value|^p$|P_value", names(df), ignore.case = TRUE)]
      for (col in p_cols)
        dt <- DT::formatStyle(dt, col,
                              backgroundColor = DT::styleInterval(0.05, c("#d4edda", "white")),
                              color           = DT::styleInterval(0.05, c("#155724", "black")))
      dt
    })

    output$mm_quality_table <- renderDT({
      res <- result_rv()
      if (is.null(res)) return(empty_dt("The model has not been estimated yet."))
      df <- extract_mm_quality(res)$quality
      if (is.null(df) || nrow(df) == 0) return(empty_dt("Quality not available."))
      dt <- datatable(format_df(df, 3), selection = "single", rownames = FALSE,
                      options = list(dom = "tip", pageLength = 100, scrollX = TRUE))

      if ("AVE" %in% names(df))
        dt <- DT::formatStyle(dt, "AVE",
                              backgroundColor = DT::styleInterval(0.499, c("#fff3cd", "#d4edda")),
                              color           = DT::styleInterval(0.499, c("#856404", "#155724")))

      for (col in intersect(names(df), c("Cronbach_alpha","Omega_McDonald","rhoC_CR","rhoA")))
        dt <- DT::formatStyle(dt, col,
                              backgroundColor = DT::styleInterval(0.699, c("#fff3cd", "#d4edda")),
                              color           = DT::styleInterval(0.699, c("#856404", "#155724")))
      dt
    })

    output$htmt_table <- renderDT({
      res <- result_rv()
      if (is.null(res)) return(empty_dt("The model has not been estimated yet."))
      df <- extract_mm_quality(res)$htmt
      if (is.null(df) || nrow(df) == 0) return(empty_dt("HTMT not available."))

      df_fmt <- format_df(df, 3)
      dt <- datatable(df_fmt, selection = "single", rownames = FALSE,
                      options = list(dom = "tip", pageLength = 100, scrollX = TRUE))

      num_cols <- names(df)[sapply(df, is.numeric)]
      for (col in num_cols) {
        dt <- DT::formatStyle(dt, col,
                              backgroundColor = DT::styleInterval(c(0.85, 0.90),
                                                                  c("#d4edda", "#fff3cd", "#f8d7da")),
                              color           = DT::styleInterval(c(0.85, 0.90),
                                                                  c("#155724", "#856404", "#721c24")))
      }
      dt
    })

    output$fornell_table <- renderDT({
      res <- result_rv()
      if (is.null(res)) return(empty_dt("The model has not been estimated yet."))
      df <- extract_mm_quality(res)$fornell
      if (is.null(df) || nrow(df) == 0) return(empty_dt("Fornell-Larcker not available."))

      df_fmt <- format_df(df, 3)
      dt <- datatable(df_fmt, selection = "single", rownames = FALSE,
                      options = list(dom = "tip", pageLength = 100, scrollX = TRUE))

      num_cols <- names(df)[sapply(df, is.numeric)]
      for (col in num_cols) {
        dt <- DT::formatStyle(dt, col,
                              backgroundColor = DT::styleInterval(c(0.70, 0.90),
                                                                  c("#d4edda", "#fff3cd", "#f8d7da")),
                              color           = DT::styleInterval(c(0.70, 0.90),
                                                                  c("#155724", "#856404", "#721c24")))
      }
      dt
    })

    output$hypothesis_table <- renderDT({
      res <- result_rv()
      if (is.null(res)) return(empty_dt("The model has not been estimated yet."))
      df <- build_hypothesis_table(res, relations_rv())
      if (is.null(df) || nrow(df) == 0)
        return(empty_dt("Hypothesis table not available. Run the model with bootstrapping to see p-values."))
      dt <- datatable(format_df(df, 3), selection = "none", rownames = FALSE,
                      options = list(dom = "tip", pageLength = 100, scrollX = TRUE))
      if ("P_value" %in% names(df))
        dt <- DT::formatStyle(dt, "P_value",
                              backgroundColor = DT::styleInterval(0.05, c("#d4edda", "white")),
                              color           = DT::styleInterval(0.05, c("#155724", "black")))
      if ("Result" %in% names(df))
        dt <- DT::formatStyle(dt, "Result",
                              backgroundColor = DT::styleEqual(c("Supported","Not supported","No resampling"),
                                                               c("#d4edda",  "#fff3cd",      "#e2e3e5")),
                              color           = DT::styleEqual(c("Supported","Not supported","No resampling"),
                                                               c("#155724",  "#856404",      "#383d41")))
      dt
    })

    output$indirect_table <- renderDT({
      res <- result_rv()
      if (is.null(res)) return(empty_dt("The model has not been estimated yet."))
      df <- extract_indirect_effects(res)
      if (is.null(df) || nrow(df) == 0)
        return(empty_dt("Indirect/total effects not available. Requires bootstrapping and at least one mediated path."))
      render_dt_with_pval_color(df, digits = 3)
    })


    # ----------------------------------------------------------------------------
    # PLOTS & TEXT SUMMARIES
    # ----------------------------------------------------------------------------
    output$estimated_model_plot <- DiagrammeR::renderGrViz({
      req(result_rv())
      plot(result_rv(), .plot_structural_model_only = !isTRUE(input$show_items_estimated), .plot_labels = TRUE, .plot_correlations = "none")
    })

    output$summary_text <- renderText({
      res <- result_rv()
      if (is.null(res)) return("The model has not been estimated yet.")

      fit_df <- extract_model_fit(res)
      srmr_txt <- if (!is.null(fit_df) && any(fit_df$Index == "SRMR")) paste0("SRMR: ", sprintf("%.3f", fit_df$Value[fit_df$Index == "SRMR"][1]), ". ") else ""

      paste0("Model estimated successfully. ",
             "Resample method: ", resample_method(), "; R = ", ifelse(resample_method() == "none", 0, n_boot()), ". ",
             "Data rows used: ", nrow(analysis_data()), ". ", srmr_txt)
    })

    output$raw_sem_summary <- renderPrint({
      req(result_rv())
      cat("================ MODEL SUMMARY ================\n")
      print(cSEM::summarize(result_rv()))
      cat("\n================ MODEL ASSESSMENT ================\n")
      print(cSEM::assess(result_rv()))
    })

    # ----------------------------------------------------------------------------
    # EXCEL EXPORT
    # ----------------------------------------------------------------------------
    build_results_list <- reactive({
      res <- result_rv()
      if (is.null(res)) return(NULL)
      bun <- analysis_bundle()
      list(
        raw_data = raw_data_rv(),
        cleaned_data = bun$data,
        construct_definitions = constructs_to_df(constructs_rv()),
        structural_relations = relations_rv(),
        construct_types = get_construct_types(res),
        model_fit = extract_model_fit(res),
        r2 = extract_r2(res),
        paths = extract_paths(res),
        loadings = extract_loadings(res),
        quality = extract_mm_quality(res)$quality,
        htmt = extract_mm_quality(res)$htmt,
        fornell = extract_mm_quality(res)$fornell,
        hypothesis_table       = build_hypothesis_table(res, relations_rv()),
        indirect_total_effects = extract_indirect_effects(res)
      )
    })

    output$download_excel <- downloadHandler(
      filename = function() paste0(sanitize_name(project_name() %||% "pls_sem_results"), ".xlsx"),
      content = function(file) {
        lst <- build_results_list()
        req(lst)
        wb <- createWorkbook()
        for (nm in names(lst)) {
          if (!is.null(lst[[nm]])) {
            addWorksheet(wb, nm)
            writeData(wb, nm, lst[[nm]])
          }
        }
        saveWorkbook(wb, file, overwrite = TRUE)
      }
    )

  })
}
