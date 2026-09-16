# ==============================================================================
# MODULE: MICOM + HENSELER MGA BY USER-DEFINED SEGMENTATION VARIABLE
# ==============================================================================
# Dependencies imported by run_app(): shiny, bslib, DT, cSEM, openxlsx
# The segmentation variable only partitions observations. The model is estimated
# with the same processed dataset used by the main PLS-SEM/PATHMOX workflow.
# ==============================================================================

mod_micom_mga_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::nav_panel(
    "MICOM & MGA",

    bslib::card(
      bslib::card_header("User-defined segmentation"),
      bslib::card_body(
        bslib::layout_columns(
          col_widths = c(5, 3, 2, 2),
          shiny::selectInput(ns("segment_var"), "Segmentation variable", choices = NULL),
          shiny::selectInput(
            ns("segment_processing"),
            "Processing",
            choices = c(
              "Use observed categories" = "factor",
              "Use text categories" = "character",
              "Create quantile groups" = "quantile",
              "Create equal-width groups" = "equal"
            ),
            selected = "factor"
          ),
          shiny::conditionalPanel(
            condition = sprintf(
              "input['%s'] == 'quantile' || input['%s'] == 'equal'",
              ns("segment_processing"), ns("segment_processing")
            ),
            shiny::numericInput(ns("segment_bins"), "Number of groups", value = 2, min = 2, max = 10, step = 1)
          ),
          shiny::numericInput(ns("min_group_n"), "Minimum n per group", value = 30, min = 10, step = 5)
        ),

        bslib::layout_columns(
          col_widths = c(4, 4, 4),
          shiny::numericInput(ns("n_boot"), "MICOM/MGA bootstrap resamples", value = 500, min = 100, step = 100),
          shiny::numericInput(ns("seed"), "Random seed", value = 12345, min = 1, step = 1),
          shiny::actionButton(ns("btn_run_micom_mga"), "Run MICOM & MGA", class = "btn-success btn-compact w-100")
        ),

        shiny::tags$hr(),
        shiny::tags$div(
          class = "small-muted",
          "The selected variable defines comparison groups only. The PLS-SEM model uses the same processed dataset and model definition as the main analysis."
        ),
        shiny::uiOutput(ns("group_summary"))
      )
    ),

    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::card(
        full_screen = TRUE,
        bslib::card_header("MICOM: compositional invariance (Step 2)"),
        bslib::card_body(DT::DTOutput(ns("micom_table")))
      ),
      bslib::card(
        full_screen = TRUE,
        bslib::card_header("Henseler MGA: structural paths"),
        bslib::card_body(DT::DTOutput(ns("mga_table")))
      )
    ),

    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::card(
        full_screen = TRUE,
        bslib::card_header("Group-specific path coefficients"),
        bslib::card_body(DT::DTOutput(ns("paths_table")))
      ),
      bslib::card(
        full_screen = TRUE,
        bslib::card_header("Analysis log"),
        bslib::card_body(DT::DTOutput(ns("log_table")))
      )
    ),

    bslib::card(
      bslib::card_header("Export"),
      bslib::card_body(
        shiny::downloadButton(ns("download_results"), "Download results (.xlsx)", class = "btn-success btn-compact")
      )
    )
  )
}

mod_micom_mga_server <- function(id, analysis_data_aug_rv, model_lavaan, result_rv,
                                 handle_inadmissibles, approach_weights, approach_paths,
                                 pls_inner_scheme, plsc_disattenuate) {
  shiny::moduleServer(id, function(input, output, session) {

    results_rv <- shiny::reactiveVal(NULL)

    empty_dt <- function(message) {
      DT::datatable(
        data.frame(Message = message, check.names = FALSE),
        rownames = FALSE,
        options = list(dom = "t", paging = FALSE),
        escape = TRUE
      )
    }

    make_group_factor <- function(data, variable, processing, bins) {
      x <- data[[variable]]

      if (processing == "factor") return(droplevels(as.factor(x)))

      if (processing == "character") {
        x <- trimws(as.character(x))
        x[!nzchar(x)] <- NA_character_
        return(droplevels(as.factor(x)))
      }

      x <- suppressWarnings(as.numeric(x))
      if (sum(!is.na(x)) < 2L || length(unique(stats::na.omit(x))) < 2L) {
        stop("The selected segmentation variable does not have enough valid numeric values.")
      }

      if (processing == "quantile") {
        probabilities <- seq(0, 1, length.out = as.integer(bins) + 1L)
        breaks <- unique(stats::quantile(x, probs = probabilities, na.rm = TRUE, type = 7))
        if (length(breaks) < 3L) {
          stop("The selected variable does not have enough distinct values for the requested quantile groups.")
        }
        return(droplevels(cut(x, breaks = breaks, include.lowest = TRUE, ordered_result = TRUE)))
      }

      droplevels(cut(x, breaks = as.integer(bins), include.lowest = TRUE, ordered_result = TRUE))
    }

    extract_structural_paths <- function(model) {
      lines <- trimws(unlist(strsplit(as.character(model), "\n", fixed = TRUE), use.names = FALSE))
      lines <- lines[nzchar(lines)]
      lines <- lines[
        grepl("~", lines, fixed = TRUE) &
          !grepl("=~", lines, fixed = TRUE) &
          !grepl("<~", lines, fixed = TRUE)
      ]

      paths <- character(0)
      for (line in lines) {
        pieces <- strsplit(line, "~", fixed = TRUE)[[1]]
        if (length(pieces) != 2L) next

        target <- trimws(pieces[1])
        predictors <- trimws(unlist(strsplit(pieces[2], "+", fixed = TRUE), use.names = FALSE))
        predictors <- predictors[nzchar(predictors)]
        if (nzchar(target) && length(predictors)) paths <- c(paths, paste0(target, " ~ ", predictors))
      }
      unique(paths)
    }

    normalize_comparisons <- function(x) {
      if (is.null(x)) return(list())
      if (!is.list(x)) return(list(Comparison = x))
      if ("none" %in% names(x)) x <- x[["none"]]
      if (!is.list(x)) return(list(Comparison = x))
      x
    }

    comparison_label <- function(x) gsub("_", " versus ", as.character(x), fixed = TRUE)

    flatten_micom <- function(micom, alpha = 0.05) {
      if (is.null(micom) || !is.null(micom$error) || is.null(micom$Step2$P_value)) return(data.frame())

      comparisons <- normalize_comparisons(micom$Step2$P_value)
      rows <- list()
      for (comparison in names(comparisons)) {
        values <- unlist(comparisons[[comparison]], recursive = TRUE, use.names = TRUE)
        if (!length(values)) next
        rows[[length(rows) + 1L]] <- data.frame(
          Comparison = comparison_label(comparison),
          Construct = names(values),
          `P-value` = as.numeric(values),
          `Compositional invariance` = ifelse(
            is.na(values), NA_character_,
            ifelse(as.numeric(values) >= alpha, "Supported", "Not supported")
          ),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      }
      if (!length(rows)) return(data.frame())
      do.call(rbind, rows)
    }

    extract_group_paths <- function(fit_mg) {
      rows <- list()
      for (group_name in names(fit_mg)) {
        path_matrix <- tryCatch(as.matrix(fit_mg[[group_name]]$Estimates$Path_estimates), error = function(e) NULL)
        if (is.null(path_matrix)) next
        paths <- as.data.frame(as.table(path_matrix), stringsAsFactors = FALSE)
        names(paths) <- c("Target", "Source", "Estimate")
        paths <- paths[!is.na(paths$Estimate) & paths$Estimate != 0, , drop = FALSE]
        if (!nrow(paths)) next
        paths$Group <- group_name
        paths$Path <- paste0(paths$Source, " -> ", paths$Target)
        rows[[length(rows) + 1L]] <- paths[, c("Group", "Path", "Estimate"), drop = FALSE]
      }
      if (!length(rows)) return(data.frame())
      do.call(rbind, rows)
    }

    extract_henseler_pvalues <- function(mga_object) {
      if (is.null(mga_object) || !is.null(mga_object$error) || is.null(mga_object$Henseler$P_value)) {
        return(data.frame())
      }

      comparisons <- normalize_comparisons(mga_object$Henseler$P_value)
      rows <- list()
      for (comparison in names(comparisons)) {
        values <- unlist(comparisons[[comparison]], recursive = TRUE, use.names = TRUE)
        if (!length(values)) next

        parameter_names <- names(values)
        display_paths <- vapply(parameter_names, function(parameter) {
          pieces <- strsplit(parameter, "~", fixed = TRUE)[[1]]
          if (length(pieces) == 2L) paste0(trimws(pieces[2]), " -> ", trimws(pieces[1])) else parameter
        }, character(1))

        rows[[length(rows) + 1L]] <- data.frame(
          Comparison = comparison_label(comparison),
          Path = display_paths,
          `P-value` = as.numeric(values),
          `Significant difference` = ifelse(
            is.na(values), NA_character_,
            ifelse(as.numeric(values) < 0.05 | as.numeric(values) > 0.95, "Yes", "No")
          ),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      }
      if (!length(rows)) return(data.frame())
      do.call(rbind, rows)
    }

    run_henseler_mga <- function(fit_mg, structural_paths, R, seed) {
      if (!length(structural_paths)) {
        return(list(results = data.frame(), error = "No structural paths were identified in the model."))
      }

      test_mgd <- get("testMGD", envir = asNamespace("cSEM"))
      result_rows <- list()
      errors <- character(0)

      for (path_syntax in structural_paths) {
        test_one <- tryCatch(
          test_mgd(
            fit_mg,
            .parameters_to_compare = path_syntax,
            .approach_mgd = "Henseler",
            .R_bootstrap = R,
            .seed = seed,
            .eval_plan = "multisession",
            .verbose = FALSE
          ),
          error = function(e) list(error = conditionMessage(e))
        )

        if (!is.null(test_one$error)) {
          errors <- c(errors, paste0(path_syntax, ": ", test_one$error))
          next
        }

        table_one <- extract_henseler_pvalues(test_one)
        if (nrow(table_one)) result_rows[[length(result_rows) + 1L]] <- table_one
      }

      list(
        results = if (length(result_rows)) do.call(rbind, result_rows) else data.frame(),
        error = if (length(errors)) paste(errors, collapse = " | ") else NULL
      )
    }

    shiny::observe({
      data <- analysis_data_aug_rv()
      choices <- if (is.null(data)) character(0) else colnames(data)
      current <- isolate(input$segment_var)
      shiny::updateSelectInput(
        session,
        "segment_var",
        choices = choices,
        selected = if (!is.null(current) && current %in% choices) current else if (length(choices)) choices[1] else character(0)
      )
    })

    shiny::observeEvent(input$btn_run_micom_mga, {
      shiny::req(analysis_data_aug_rv(), model_lavaan())

      data <- as.data.frame(analysis_data_aug_rv())
      model <- model_lavaan()

      if (identical(model, "No constructs defined yet.")) {
        shiny::showNotification("Define constructs and structural relations first.", type = "error")
        return()
      }

      if (is.null(input$segment_var) || !nzchar(input$segment_var)) {
        shiny::showNotification("Select a segmentation variable.", type = "error")
        return()
      }

      analysis_out <- tryCatch({
        grouping <- make_group_factor(
          data = data,
          variable = input$segment_var,
          processing = input$segment_processing,
          bins = input$segment_bins
        )

        keep_nonmissing <- !is.na(grouping)
        analysis_data <- data[keep_nonmissing, , drop = FALSE]
        grouping <- droplevels(grouping[keep_nonmissing])

        observed_counts <- table(grouping)
        included_levels <- names(observed_counts)[observed_counts >= input$min_group_n]
        excluded_levels <- names(observed_counts)[observed_counts < input$min_group_n]

        if (length(included_levels) < 2L) {
          stop(sprintf(
            "At least two groups with n >= %d are required. Observed group sizes: %s",
            input$min_group_n,
            paste(names(observed_counts), observed_counts, sep = "=", collapse = "; ")
          ))
        }

        keep_included <- grouping %in% included_levels
        analysis_data <- analysis_data[keep_included, , drop = FALSE]
        grouping <- droplevels(grouping[keep_included])

        # Deliberately retain the complete processed data frame. This mirrors
        # PATHMOX, preserves all model information and avoids using the internal
        # parseModel representation as a source of indicator names.
        group_data <- split(analysis_data, grouping, drop = TRUE)
        group_data <- lapply(group_data, as.data.frame)

        shiny::withProgress(message = "Running MICOM and MGA...", value = 0.05, {
          shiny::incProgress(0.20, detail = "Estimating multigroup PLS-SEM model")

          fit_mg <- cSEM::csem(
            .data = group_data,
            .model = model,
            .resample_method = "none",
            .handle_inadmissibles = handle_inadmissibles(),
            .approach_weights = approach_weights(),
            .approach_paths = approach_paths(),
            .PLS_weight_scheme_inner = pls_inner_scheme(),
            .disattenuate = plsc_disattenuate(),
            .eval_plan = "multisession"
          )

          shiny::incProgress(0.25, detail = "Testing compositional invariance (MICOM)")
          micom <- tryCatch(
            cSEM::testMICOM(
              fit_mg,
              .R = max(50L, as.integer(input$n_boot)),
              .seed = as.integer(input$seed),
              .approach_p_adjust = "none",
              .verbose = FALSE
            ),
            error = function(e) list(error = conditionMessage(e))
          )

          shiny::incProgress(0.35, detail = "Testing structural path differences (Henseler MGA)")
          structural_paths <- extract_structural_paths(model)
          mga_run <- run_henseler_mga(
            fit_mg = fit_mg,
            structural_paths = structural_paths,
            R = max(50L, as.integer(input$n_boot)),
            seed = as.integer(input$seed)
          )

          shiny::incProgress(0.10, detail = "Preparing output tables")
          group_counts <- data.frame(
            Group = names(group_data),
            N = vapply(group_data, nrow, integer(1)),
            stringsAsFactors = FALSE
          )

          group_text <- paste0(group_counts$Group, " (n=", group_counts$N, ")")
          log_table <- data.frame(
            Item = c(
              "Segmentation variable",
              "Processing",
              "Groups analysed",
              "Observations analysed",
              "Excluded groups",
              "MICOM bootstrap resamples",
              "MGA bootstrap resamples",
              "Structural paths tested"
            ),
            Value = c(
              input$segment_var,
              input$segment_processing,
              paste(group_text, collapse = "; "),
              as.character(nrow(analysis_data)),
              if (length(excluded_levels)) paste(excluded_levels, collapse = "; ") else "None",
              as.character(max(50L, as.integer(input$n_boot))),
              as.character(max(50L, as.integer(input$n_boot))),
              paste(structural_paths, collapse = "; ")
            ),
            stringsAsFactors = FALSE
          )

          list(
            group_counts = group_counts,
            excluded_groups = excluded_levels,
            micom = micom,
            micom_table = flatten_micom(micom),
            mga_error = mga_run$error,
            mga_table = mga_run$results,
            paths_table = extract_group_paths(fit_mg),
            log_table = log_table
          )
        })
      }, error = function(e) list(error = conditionMessage(e)))

      if (!is.null(analysis_out$error)) {
        results_rv(NULL)
        shiny::showNotification(
          paste("MICOM/MGA analysis failed:", analysis_out$error),
          type = "error",
          duration = 10
        )
      } else {
        results_rv(analysis_out)
        if (!is.null(analysis_out$mga_error) && !nrow(analysis_out$mga_table)) {
          shiny::showNotification(
            "MICOM completed, but MGA could not be calculated for the selected paths.",
            type = "warning",
            duration = 10
          )
        } else if (!is.null(analysis_out$micom$error)) {
          shiny::showNotification(
            "MGA completed, but MICOM could not be calculated.",
            type = "warning",
            duration = 10
          )
        } else {
          shiny::showNotification("MICOM and Henseler MGA completed successfully.", type = "message")
        }
      }
    })

    output$group_summary <- shiny::renderUI({
      out <- results_rv()
      if (is.null(out)) {
        return(shiny::tags$div(
          class = "small-muted",
          style = "margin-top: 8px;",
          "Select a segmentation variable and run the analysis."
        ))
      }

      labels <- paste0(out$group_counts$Group, " (n=", out$group_counts$N, ")")
      shiny::tags$div(
        class = "small-muted",
        style = "margin-top: 8px;",
        shiny::tags$b("Groups analysed: "),
        paste(labels, collapse = " · "),
        if (length(out$excluded_groups)) {
          shiny::tags$span(" | Excluded for insufficient size: ", paste(out$excluded_groups, collapse = ", "))
        }
      )
    })

    output$micom_table <- DT::renderDT({
      out <- results_rv()
      if (is.null(out) || !nrow(out$micom_table)) {
        message <- if (!is.null(out) && !is.null(out$micom$error)) {
          paste("MICOM failed:", out$micom$error)
        } else {
          "No MICOM Step 2 results available."
        }
        return(empty_dt(message))
      }

      DT::datatable(
        out$micom_table,
        rownames = FALSE,
        options = list(pageLength = 15, scrollX = TRUE, dom = "tip")
      ) |>
        DT::formatRound("P-value", 3) |>
        DT::formatStyle(
          "Compositional invariance",
          color = DT::styleEqual(c("Supported", "Not supported"), c("#198754", "#DC3545"))
        )
    })

    output$mga_table <- DT::renderDT({
      out <- results_rv()
      if (is.null(out) || !nrow(out$mga_table)) {
        message <- if (!is.null(out) && !is.null(out$mga_error)) {
          paste("MGA failed:", out$mga_error)
        } else {
          "No Henseler MGA results available."
        }
        return(empty_dt(message))
      }

      DT::datatable(
        out$mga_table,
        rownames = FALSE,
        options = list(pageLength = 15, scrollX = TRUE, dom = "tip")
      ) |>
        DT::formatRound("P-value", 3) |>
        DT::formatStyle(
          "Significant difference",
          color = DT::styleEqual("Yes", "#DC3545")
        )
    })

    output$paths_table <- DT::renderDT({
      out <- results_rv()
      if (is.null(out) || is.null(out$paths_table) || !nrow(out$paths_table)) {
        return(empty_dt("No group-specific structural path estimates are available."))
      }

      DT::datatable(
        out$paths_table,
        rownames = FALSE,
        options = list(pageLength = 15, scrollX = TRUE, dom = "tip")
      ) |>
        DT::formatRound("Estimate", 3)
    })

    output$log_table <- DT::renderDT({
      out <- results_rv()
      if (is.null(out)) return(empty_dt("No analysis has been run."))

      DT::datatable(
        out$log_table,
        rownames = FALSE,
        options = list(dom = "t", paging = FALSE),
        escape = TRUE
      )
    })

    output$download_results <- shiny::downloadHandler(
      filename = function() paste0("MICOM_MGA_", Sys.Date(), ".xlsx"),
      content = function(file) {
        out <- results_rv()
        shiny::req(out)

        workbook <- openxlsx::createWorkbook()
        sheets <- list(
          Groups = out$group_counts,
          MICOM_Step2 = out$micom_table,
          Henseler_MGA = out$mga_table,
          Group_paths = out$paths_table,
          Log = out$log_table
        )

        for (sheet_name in names(sheets)) {
          sheet_data <- sheets[[sheet_name]]
          if (is.null(sheet_data) || !ncol(sheet_data)) {
            sheet_data <- data.frame(Message = "No results available.")
          }
          openxlsx::addWorksheet(workbook, sheet_name)
          openxlsx::writeData(workbook, sheet_name, sheet_data)
          openxlsx::setColWidths(workbook, sheet_name, cols = seq_len(ncol(sheet_data)), widths = "auto")
        }

        openxlsx::saveWorkbook(workbook, file, overwrite = TRUE)
      }
    )

    invisible(list(results = shiny::reactive(results_rv())))
  })
}
