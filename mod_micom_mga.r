# ==============================================================================
# MODULE: MICOM + HENSELER MGA BY USER-DEFINED SEGMENTATION VARIABLE
# Requires: shiny, bslib, DT, cSEM
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
            ns("segment_processing"), "Processing",
            choices = c(
              "As factor" = "factor",
              "As character" = "character",
              "Quantile groups" = "quantile",
              "Equal-width groups" = "equal"
            ),
            selected = "factor"
          ),
          shiny::numericInput(ns("segment_bins"), "Groups / bins", value = 2, min = 2, max = 10, step = 1),
          shiny::numericInput(ns("min_group_n"), "Minimum n", value = 30, min = 10, step = 5)
        ),
        bslib::layout_columns(
          col_widths = c(4, 4, 4),
          shiny::numericInput(ns("n_boot"), "MICOM/MGA bootstrap resamples", value = 500, min = 100, step = 100),
          shiny::numericInput(ns("seed"), "Random seed", value = 12345, min = 1, step = 1),
          shiny::actionButton(ns("btn_run_micom_mga"), "Run MICOM & MGA", class = "btn-success btn-compact w-100")
        ),
        shiny::tags$hr(),
        shiny::uiOutput(ns("group_summary"))
      )
    ),
    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::card(full_screen = TRUE, bslib::card_header("MICOM: compositional invariance (Step 2)"), bslib::card_body(DT::DTOutput(ns("micom_table")))),
      bslib::card(full_screen = TRUE, bslib::card_header("Henseler MGA: structural paths"), bslib::card_body(DT::DTOutput(ns("mga_table"))))
    ),
    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::card(full_screen = TRUE, bslib::card_header("Group-specific path coefficients"), bslib::card_body(DT::DTOutput(ns("paths_table")))),
      bslib::card(full_screen = TRUE, bslib::card_header("Analysis log"), bslib::card_body(DT::DTOutput(ns("log_table"))))
    ),
    bslib::card(
      bslib::card_header("Export"),
      bslib::card_body(shiny::downloadButton(ns("download_results"), "Download results (.xlsx)", class = "btn-success btn-compact"))
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
        data.frame(Message = message), rownames = FALSE,
        options = list(dom = "t", paging = FALSE), escape = TRUE
      )
    }

    segment_data <- function(data, variable, processing, bins) {
      x <- data[[variable]]
      if (processing == "factor") {
        grp <- as.factor(x)
      } else if (processing == "character") {
        grp <- as.factor(trimws(as.character(x)))
      } else if (processing == "quantile") {
        x <- suppressWarnings(as.numeric(x))
        probs <- seq(0, 1, length.out = bins + 1)
        cuts <- unique(stats::quantile(x, probs = probs, na.rm = TRUE, type = 7))
        if (length(cuts) < 3) stop("The selected variable does not have enough distinct values for the requested quantile groups.")
        grp <- cut(x, breaks = cuts, include.lowest = TRUE, ordered_result = TRUE)
      } else {
        x <- suppressWarnings(as.numeric(x))
        if (length(unique(stats::na.omit(x))) < 2) stop("The selected variable does not have enough distinct numeric values.")
        grp <- cut(x, breaks = bins, include.lowest = TRUE, ordered_result = TRUE)
      }
      as.factor(grp)
    }

    extract_pvalues <- function(x) {
      if (is.null(x)) return(list())
      if (!is.list(x)) return(list(Comparison = x))
      if ("none" %in% names(x)) x <- x[["none"]]
      if (!is.list(x)) return(list(Comparison = x))
      x
    }

    comparison_label <- function(x) gsub("_", " versus ", as.character(x), fixed = TRUE)

    flatten_micom <- function(micom, alpha = 0.05) {
      if (is.null(micom) || !is.null(micom$error)) return(data.frame())
      pvals <- tryCatch(extract_pvalues(micom$Step2$P_value), error = function(e) list())
      rows <- list()
      for (cmp in names(pvals)) {
        vals <- pvals[[cmp]]
        vals <- unlist(vals, recursive = TRUE, use.names = TRUE)
        if (!length(vals)) next
        rows[[length(rows) + 1]] <- data.frame(
          Comparison = comparison_label(cmp),
          Construct = names(vals),
          `P-value` = as.numeric(vals),
          `Compositional invariance` = ifelse(is.na(vals), NA, ifelse(vals >= alpha, "Supported", "Not supported")),
          check.names = FALSE, stringsAsFactors = FALSE
        )
      }
      if (!length(rows)) return(data.frame())
      do.call(rbind, rows)
    }

    flatten_mga <- function(mga) {
      if (is.null(mga) || !is.null(mga$error) || is.null(mga$Henseler$P_value)) return(data.frame())
      pvals <- extract_pvalues(mga$Henseler$P_value)
      rows <- list()
      for (cmp in names(pvals)) {
        vals <- pvals[[cmp]]
        vals <- unlist(vals, recursive = TRUE, use.names = TRUE)
        if (!length(vals)) next
        param <- names(vals)
        path <- vapply(param, function(z) {
          s <- strsplit(z, "~", fixed = TRUE)[[1]]
          if (length(s) == 2) paste0(trimws(s[2]), " -> ", trimws(s[1])) else z
        }, character(1))
        rows[[length(rows) + 1]] <- data.frame(
          Comparison = comparison_label(cmp),
          Path = path,
          `P-value` = as.numeric(vals),
          `Significant difference` = ifelse(is.na(vals), NA, ifelse(vals < 0.05 | vals > 0.95, "Yes", "No")),
          check.names = FALSE, stringsAsFactors = FALSE
        )
      }
      if (!length(rows)) return(data.frame())
      do.call(rbind, rows)
    }

    extract_group_paths <- function(fit) {
      est <- tryCatch(fit$Estimates$Path_estimates, error = function(e) NULL)
      if (is.null(est)) return(data.frame())
      out <- as.data.frame(as.table(as.matrix(est)), stringsAsFactors = FALSE)
      names(out) <- c("Target", "Source", "Estimate")
      out <- out[!is.na(out$Estimate) & out$Estimate != 0, , drop = FALSE]
      if (!nrow(out)) return(data.frame())
      out$Path <- paste0(out$Source, " -> ", out$Target)
      out[, c("Path", "Estimate"), drop = FALSE]
    }

    shiny::observe({
      data <- analysis_data_aug_rv()
      choices <- if (is.null(data)) character(0) else colnames(data)
      selected <- isolate(input$segment_var)
      shiny::updateSelectInput(
        session, "segment_var", choices = choices,
        selected = if (!is.null(selected) && selected %in% choices) selected else if (length(choices)) choices[1] else character(0)
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

      out <- tryCatch({
        grp <- segment_data(data, input$segment_var, input$segment_processing, input$segment_bins)
        keep <- !is.na(grp)
        data2 <- data[keep, , drop = FALSE]
        grp <- droplevels(grp[keep])
        counts <- table(grp)
        valid_levels <- names(counts)[counts >= input$min_group_n]
        excluded <- names(counts)[counts < input$min_group_n]

        if (length(valid_levels) < 2) {
          stop(sprintf("At least two groups with n >= %d are required. Current group sizes: %s",
                       input$min_group_n, paste(names(counts), counts, sep = "=", collapse = "; ")))
        }

        idx <- grp %in% valid_levels
        data2 <- data2[idx, , drop = FALSE]
        grp <- droplevels(grp[idx])
        group_data <- split(data2, grp, drop = TRUE)
        group_data <- lapply(group_data, as.data.frame)

        shiny::withProgress(message = "Running MICOM and MGA...", value = 0.15, {
          shiny::incProgress(0.25, detail = "Estimating the multigroup PLS-SEM model")
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

          shiny::incProgress(0.35, detail = "Testing compositional invariance (MICOM)")
          micom <- tryCatch(
            cSEM::testMICOM(
              fit_mg, .R = max(50, input$n_boot), .seed = input$seed,
              .approach_p_adjust = "none", .verbose = FALSE
            ),
            error = function(e) list(error = conditionMessage(e))
          )

          shiny::incProgress(0.20, detail = "Testing path differences (Henseler MGA)")
          structural_lines <- trimws(strsplit(model, "\n", fixed = TRUE)[[1]])
          structural_lines <- structural_lines[
            grepl("~", structural_lines, fixed = TRUE) &
              !grepl("=~", structural_lines, fixed = TRUE) &
              !grepl("<~", structural_lines, fixed = TRUE)
          ]

          parameters <- paste(structural_lines, collapse = "\n")


          mga <- if (nzchar(parameters)) {
            tryCatch(
              get("testMGD", envir = asNamespace("cSEM"))(
                fit_mg,
                .parameters_to_compare = parameters,
                .approach_mgd = "Henseler",
                .R_bootstrap = max(50, input$n_boot),
                .eval_plan = "multisession",
                .verbose = FALSE
              ),
              error = function(e) list(error = conditionMessage(e))
            )
          } else list(error = "No structural paths were found in the model.")

          shiny::incProgress(0.05, detail = "Preparing results")
          paths <- do.call(rbind, lapply(names(fit_mg), function(g) {
            p <- extract_group_paths(fit_mg[[g]])
            if (!nrow(p)) return(NULL)
            p$Group <- g
            p[, c("Group", "Path", "Estimate"), drop = FALSE]
          }))

          log <- data.frame(
            Item = c("Segmentation variable", "Processing", "Groups analysed", "Observations analysed", "Excluded groups", "MICOM bootstrap resamples", "MGA bootstrap resamples"),
            Value = c(
              input$segment_var,
              input$segment_processing,
              paste(names(group_data), vapply(group_data, nrow, integer(1)), sep = " (n=", collapse = "; "),
              nrow(data2),
              if (length(excluded)) paste(excluded, collapse = "; ") else "None",
              max(50, input$n_boot),
              max(50, input$n_boot)
            ),
            stringsAsFactors = FALSE
          )
          log$Value[3] <- paste0(gsub("$", ")", log$Value[3]))

          list(
            group_counts = data.frame(Group = names(group_data), N = vapply(group_data, nrow, integer(1)), stringsAsFactors = FALSE),
            excluded = excluded,
            micom = micom,
            mga = mga,
            micom_table = flatten_micom(micom),
            mga_table = flatten_mga(mga),
            paths_table = paths,
            log_table = log
          )
        })
      }, error = function(e) list(error = conditionMessage(e)))

      if (!is.null(out$error)) {
        results_rv(NULL)
        shiny::showNotification(paste("MICOM/MGA analysis failed:", out$error), type = "error", duration = 10)
      } else {
        results_rv(out)
        shiny::showNotification("MICOM and Henseler MGA completed successfully.", type = "message")
      }
    })

    output$group_summary <- shiny::renderUI({
      x <- results_rv()
      if (is.null(x)) return(shiny::tags$div(class = "small-muted", "Select a segmentation variable and run the analysis."))
      labels <- paste0(x$group_counts$Group, " (n=", x$group_counts$N, ")")
      shiny::tags$div(
        class = "small-muted",
        shiny::tags$b("Groups analysed: "), paste(labels, collapse = " · "),
        if (length(x$excluded)) shiny::tags$span(" | Excluded for insufficient size: ", paste(x$excluded, collapse = ", "))
      )
    })

    output$micom_table <- DT::renderDT({
      x <- results_rv()
      if (is.null(x) || !nrow(x$micom_table)) return(empty_dt(if (!is.null(x$micom$error)) paste("MICOM failed:", x$micom$error) else "No MICOM Step 2 results available."))
      DT::datatable(x$micom_table, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE, dom = "tip")) |>
        DT::formatRound("P-value", 3) |>
        DT::formatStyle("Compositional invariance", color = DT::styleEqual(c("Supported", "Not supported"), c("#198754", "#DC3545")))
    })

    output$mga_table <- DT::renderDT({
      x <- results_rv()
      if (is.null(x) || !nrow(x$mga_table)) return(empty_dt(if (!is.null(x$mga$error)) paste("MGA failed:", x$mga$error) else "No Henseler MGA results available."))
      DT::datatable(x$mga_table, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE, dom = "tip")) |>
        DT::formatRound("P-value", 3) |>
        DT::formatStyle("Significant difference", color = DT::styleEqual("Yes", "#DC3545"))
    })

    output$paths_table <- DT::renderDT({
      x <- results_rv()
      if (is.null(x) || is.null(x$paths_table) || !nrow(x$paths_table)) return(empty_dt("No group-specific path estimates available."))
      DT::datatable(x$paths_table, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE, dom = "tip")) |>
        DT::formatRound("Estimate", 3)
    })

    output$log_table <- DT::renderDT({
      x <- results_rv()
      if (is.null(x)) return(empty_dt("No analysis has been run."))
      DT::datatable(x$log_table, rownames = FALSE, options = list(dom = "t", paging = FALSE), escape = TRUE)
    })

    output$download_results <- shiny::downloadHandler(
      filename = function() paste0("MICOM_MGA_", Sys.Date(), ".xlsx"),
      content = function(file) {
        x <- results_rv()
        shiny::req(x)
        wb <- openxlsx::createWorkbook()
        openxlsx::addWorksheet(wb, "Groups")
        openxlsx::writeData(wb, "Groups", x$group_counts)
        openxlsx::addWorksheet(wb, "MICOM_Step2")
        openxlsx::writeData(wb, "MICOM_Step2", x$micom_table)
        openxlsx::addWorksheet(wb, "Henseler_MGA")
        openxlsx::writeData(wb, "Henseler_MGA", x$mga_table)
        openxlsx::addWorksheet(wb, "Group_paths")
        openxlsx::writeData(wb, "Group_paths", x$paths_table)
        openxlsx::addWorksheet(wb, "Log")
        openxlsx::writeData(wb, "Log", x$log_table)
        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      }
    )

    invisible(list(results = shiny::reactive(results_rv())))
  })
}
