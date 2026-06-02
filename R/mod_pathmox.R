# ==============================================================================
# MODULE: PATHMOX ANALYSIS
# ==============================================================================

# ------------------------------------------------------------------------------
# CUSTOM ALGORITHM: PATHMOX + MICOM INTEGRADO
# ------------------------------------------------------------------------------
pathmoxMICOM <- function(.model, .data, .catvar, .scheme = "path", .consistent = FALSE, 
                         .alpha = 0.05, .deep = 2, .size = 0.1, .size_candidate = 50, 
                         .tree = TRUE, .test_invariance = TRUE, .alpha_invariance = 0.05,
                         .require_critical_invariance = TRUE, .boot_micom = 200, 
                         .verbose = TRUE) {
  
  # --- Dependencies ---
  requireNamespace("cSEM", quietly = TRUE)
  requireNamespace("genpathmox", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  
  if(.verbose) {
    cat("\n*** STARTING PATHMOX-MICOM ALGORITHM ***")
    cat(sprintf("\n - Alpha: %.2f | Max Depth: %d | Min Size: %.2f", .alpha, .deep, .size))
    cat(sprintf("\n - Invariance (MICOM): %s | Require Critical Inv: %s | Bootstraps: %d\n", .test_invariance, .require_critical_invariance, .boot_micom))
  }
  
  # --- Initialization & Validation ---
  .data <- as.data.frame(.data)
  
  suppressWarnings({
    check <- genpathmox:::check_arg_mox(.model = .model, .data = .data, .catvar = .catvar, 
                                        .scheme = .scheme, .consistent = .consistent, .alpha = .alpha, 
                                        .deep = .deep, .size = .size, .size_candidate = .size_candidate, 
                                        .tree = .tree)
  })
  
  .model <- check$model; .data <- check$data; .catvar <- check$catvar; .scheme <- check$scheme
  .consistent <- check$consistent; .alpha <- check$alpha; .deep <- check$deep
  .size <- check$size; .size_candidate <- check$size_candidate; .tree <- check$tree
  
  .mod_parsed <- cSEM::parseModel(.model)
  .inner <- as.matrix(.mod_parsed$structural)
  .indicators <- unique(unlist(.mod_parsed$measurement))
  
  .par_mode <- list(alpha = .alpha, size = .size, deep = .deep, data = .data, catvar = .catvar, 
                    inner = .inner, scheme = .scheme, consistent = .consistent)
  
  .min.ind.node <- genpathmox:::percent.node(.data, .size)
  .hybrid <- list(); .id <- 1
  .dim_row <- dim(.data)[1]
  .t <- methods::new("moxtree", id = 0)
  .root <- new("node", id = 1, elements = seq(1:.dim_row), father = 0, childs = 0)
  .new_nodes <- list(.root)
  .rejected_splits <- data.frame(variable=character(), level=character(), reason=character(), info=character())
  
  # --- Main Partitioning Loop ---
  while (length(.new_nodes) > 0) {
    .n <- .new_nodes[[1]]
    
    if(.verbose) cat(sprintf("\n[NODE %s] N=%d (Depth: %d)", .n@id, length(.n@elements), genpathmox:::showDeepth(.n)))
    
    if (length(.n@elements) >= .min.ind.node$min.n.ind && genpathmox:::showDeepth(.n) < .deep) {
      
      .d <- .data[.n@elements, ]
      .s_candidates <- .catvar[.n@elements, , drop=FALSE]
      
      .vars_to_remove <- sapply(.s_candidates, function(x) length(unique(x)) < 2)
      .s_candidates <- .s_candidates[, !.vars_to_remove, drop=FALSE]
      
      .cat_cols <- sapply(.s_candidates, is.factor)
      if(sum(.cat_cols) > 0) .s_candidates[.cat_cols] <- lapply(.s_candidates[.cat_cols], factor)
      
      .split_success <- FALSE
      
      while(ncol(.s_candidates) > 0 && !.split_success) {
        set.seed(12345) 
        .tmp <- genpathmox:::partopt(.d, .s_candidates, .inner, .model, .scheme, .consistent, .size_candidate)
        
        if (is.null(.tmp$pvl.opt) || is.na(.tmp$pvl.opt) || .tmp$pvl.opt > .alpha) {
          if(.verbose) cat(sprintf("\n[STOP] No significant variable (p > %.2f).", .alpha))
          break 
        }
        
        .variable <- .tmp$variable.opt
        .level <- .tmp$level.opt
        .candidates <- .tmp$candidates
        .modtwo <- .tmp$modtwo.opt
        
        if(.verbose) cat(sprintf("\n>>> Testing candidate: '%s' (p=%.4f)", .variable, .tmp$pvl.opt))
        
        .mod <- genpathmox:::test.partition(.d, .inner, .model, .scheme, .consistent, .modtwo, .alpha)
        
        if (any(!is.na(.mod$pvc) & .mod$pvc <= .alpha)) {
          .invariance_passed <- TRUE 
          
          if (.test_invariance) {
            if(.verbose) cat("\n... Checking MICOM Invariance ...") 
            
            .idx_g1 <- .n@elements[which(.modtwo == 1)]
            .idx_g2 <- .n@elements[which(.modtwo == 2)]
            
            if(length(.idx_g1) > 20 && length(.idx_g2) > 20) {
              .dt1 <- .data[.idx_g1, ]; .dt2 <- .data[.idx_g2, ]
              
              .sd_g1 <- sapply(.dt1[, .indicators[which(.indicators %in% names(.dt1))], drop=FALSE], sd, na.rm=TRUE)
              .sd_g2 <- sapply(.dt2[, .indicators[which(.indicators %in% names(.dt2))], drop=FALSE], sd, na.rm=TRUE)
              
              if(any(.sd_g1 == 0) || any(.sd_g2 == 0)) {
                .invariance_passed <- FALSE
                .rejected_splits <- rbind(.rejected_splits, data.frame(variable=.variable, level=as.character(.level), reason="Zero Variance", info="Item variance is 0"))
              } else {
                .micom_res <- tryCatch({
                  suppressMessages(suppressWarnings({
                    .fit_mg <- cSEM::csem(.data = list(G1=.dt1, G2=.dt2), .model = .model, .PLS_weight_scheme_inner = .scheme, .disattenuate = .consistent)
                  }))
                  
                  if(inherits(.fit_mg, "cSEMResults_multi")) {
                    suppressWarnings({
                      .res <- cSEM::testMICOM(.fit_mg, .R = .boot_micom, .seed = 12345, .verbose = FALSE, .approach_p_adjust = "none")
                    })
                    .pvals_obj <- if(is.list(.res$Step2$P_value)) .res$Step2$P_value[[1]] else .res$Step2$P_value
                    .pvals_vec <- if(is.list(.pvals_obj)) as.numeric(.pvals_obj[[1]]) else as.numeric(.pvals_obj)
                    if(!is.null(.res$Step2$P_value[[1]])) names(.pvals_vec) <- names(.res$Step2$P_value[[1]][[1]]) 
                    
                    list(ok=TRUE, fails=names(.pvals_vec)[!is.na(.pvals_vec) & .pvals_vec < .alpha_invariance])
                  } else { list(ok=FALSE, msg="Model fitting failed") }
                }, error = function(e) list(ok=FALSE, msg=conditionMessage(e)))
                
                if(.micom_res$ok) {
                  if(length(.micom_res$fails) > 0) {
                    
                    # --- CONDITIONAL LOGIC ---
                    significant_paths <- names(.mod$pvc)[!is.na(.mod$pvc) & .mod$pvc <= .alpha]
                    critical_conflict <- FALSE
                    if(length(significant_paths) > 0) {
                      critical_conflict <- any(sapply(.micom_res$fails, function(f) {
                        any(grepl(paste0("\\b", f, "\\b"), significant_paths))
                      }))
                    }
                    
                    if (isTRUE(critical_conflict) && .require_critical_invariance) {
                      .invariance_passed <- FALSE
                      if(.verbose) cat(sprintf("\n[FAIL] Critical Invariance failed on paths: %s", paste(.micom_res$fails, collapse=",")))
                      .rejected_splits <- rbind(.rejected_splits, data.frame(variable=.variable, level=as.character(.level), reason="Invariance Fail (Critical)", info=paste(.micom_res$fails, collapse=",")))
                    } else {
                      if(.verbose) cat(sprintf(" [OK - Issues found in non-critical constructs: %s]", paste(.micom_res$fails, collapse=",")))
                    }
                  } else { if(.verbose) cat(" [OK]") }
                } else {
                  .invariance_passed <- FALSE
                  if(.verbose) cat(sprintf("\n[ERROR] %s", .micom_res$msg))
                  .rejected_splits <- rbind(.rejected_splits, data.frame(variable=.variable, level=as.character(.level), reason="Technical Error", info=.micom_res$msg))
                }
              }
            } else {
              .invariance_passed <- FALSE
              if(.verbose) cat("\n[FAIL] Resulting groups too small (< 20).")
            }
          } 
          
          if (.invariance_passed) {
            if(.verbose) cat(sprintf("\n[SUCCESS] Splitting node by '%s'.", .variable))
            for (i in 1:2) {
              .child_id <- (.n@id * 2) + i - 1
              .n@childs[i] <- .child_id
              .new_nodes[[length(.new_nodes) + 1]] <- new("node", id = .child_id, elements = .n@elements[which(.modtwo == i)], father = .n@id)
            }
            .n@info <- new("info", variable = .variable, level = .level, fgstatistic = .mod$Fg, fpvalg = .mod$pvg, 
                           fcstatistic = .mod$Fc, fpvalc = .mod$pvc, candidates = .candidates)
            .split_success <- TRUE 
          } else { .s_candidates[[.variable]] <- NULL }
        } else {
          if(.verbose) cat("\n[FAIL] Coefficients not significant.")
          .s_candidates[[.variable]] <- NULL
        }
      } 
    }
    
    .t@nodes[[length(.t@nodes) + 1]] <- .n
    .new_nodes[1] <- NULL
  }
  
  if(.verbose) cat("\n\n*** ALGORITHM FINISHED ***\n")
  
  if (length(.t@nodes) == 1) return(list(root = genpathmox:::root.tree(.t), par_mode = .par_mode, rejected_splits = .rejected_splits))
  
  .MOX <- genpathmox:::mox.tree(.t)
  .terminal <- genpathmox:::terminal.tree(.t)
  .terminal_paths <- NULL
  
  for (i in 2:length(.terminal)) {
    .data.node <- cbind(.data[.terminal[[i]], ], hybrid_var = rep(paste(names(.terminal)[i]), length(.terminal[[i]])))
    .hybrid[[length(.hybrid) + 1]] <- .data.node
    suppressMessages({
      .fit_temp <- cSEM::csem(.data[.terminal[[i]], ], .model, .PLS_weight_scheme_inner = .scheme, .disattenuate = .consistent)
    })
    .pls_path <- as.matrix(genpathmox:::.path(.fit_temp$Estimates$Path_estimates))
    rownames(.pls_path) <- genpathmox:::element(.inner)
    .pls_R2 <- as.matrix(.fit_temp$Estimates$R2); rownames(.pls_R2) <- paste("R^2", rownames(.pls_R2))
    .terminal_paths <- cbind(.terminal_paths, rbind(.pls_path, .pls_R2))
  }
  
  names(.hybrid) <- names(.terminal)[2:length(.terminal)]
  colnames(.terminal_paths) <- names(.terminal)[2:length(.terminal)]
  
  .res <- list(MOX = .MOX, terminal_paths = round(.terminal_paths, 4), 
               var_imp = genpathmox:::var_imp_mox(genpathmox:::candidates.tree(.t), .catvar), 
               Fg.r = genpathmox:::fglobal.tree(.t), Fc.r = genpathmox:::fcoef.tree(.t), hybrid = .hybrid, 
               other = list(candidates = genpathmox:::candidates.tree(.t), root = genpathmox:::root.tree(.t), 
                            nodes = genpathmox:::nodes.tree(.t), terminal = .terminal, par_mode = .par_mode), 
               rejected_splits = .rejected_splits)
  class(.res) <- "plstree"
  
  if (.tree) try(plot(.res, .main = "PATHMOX-MICOM SMART TREE") , silent=TRUE)
  
  return(.res)
}

# ------------------------------------------------------------------------------
# UI Y SERVER DEL MÓDULO
# ------------------------------------------------------------------------------

mod_pathmox_ui <- function(id) {
  ns <- NS(id)
  
  nav_panel("PATHMOX analysis",
            card(
              card_header("PATHMOX Configuration"),
              card_body(
                card(
                  card_header("Segmentation variables"),
                  card_body(
                    layout_columns(
                      col_widths = c(4, 4, 4),
                      selectInput(ns("pathmox_segvar_selected"), "Variable", choices = NULL),
                      selectInput(ns("pathmox_segvar_treatment"), "Processing", choices = c("As factor" = "factor", "Discretize into quantiles" = "quantile", "Discretize into equal width bins" = "equal", "Keep as character" = "character"), selected = "factor"),
                      conditionalPanel(
                        condition = sprintf("input['%s'] == 'quantile' || input['%s'] == 'equal'", ns("pathmox_segvar_treatment"), ns("pathmox_segvar_treatment")), 
                        numericInput(ns("pathmox_segvar_bins"), "Number of bins", value = 3, min = 2, max = 10, step = 1)
                      )
                    ),
                    layout_columns(
                      col_widths = c(8, 4),
                      conditionalPanel(
                        condition = sprintf("input['%s'] == 'factor' || input['%s'] == 'character'", ns("pathmox_segvar_treatment"), ns("pathmox_segvar_treatment")), 
                        checkboxInput(ns("pathmox_segvar_ordered"), "Treat as ordered", value = FALSE)
                      ), 
                      div()
                    ),
                    conditionalPanel(
                      condition = sprintf("(input['%s'] == 'factor' || input['%s'] == 'character') && input['%s'] == true", ns("pathmox_segvar_treatment"), ns("pathmox_segvar_treatment"), ns("pathmox_segvar_ordered")),
                      layout_columns(
                        col_widths = c(5, 2, 5),
                        div(tags$b("Available levels"), selectInput(ns("pathmox_levels_available"), NULL, choices = NULL, multiple = TRUE, selectize = FALSE, size = 8, width = "100%")),
                        div(style = "padding-top: 42px;", actionButton(ns("btn_pathmox_level_right"), ">", class = "btn-outline-primary btn-compact"), tags$br(), tags$br(), actionButton(ns("btn_pathmox_level_left"), "<", class = "btn-outline-secondary btn-compact"), tags$br(), tags$br(), actionButton(ns("btn_pathmox_level_up"), "Up", class = "btn-outline-dark btn-compact"), tags$br(), tags$br(), actionButton(ns("btn_pathmox_level_down"), "Down", class = "btn-outline-dark btn-compact")),
                        div(tags$b("Ordered levels"), selectInput(ns("pathmox_levels_selected"), NULL, choices = NULL, multiple = TRUE, selectize = FALSE, size = 8, width = "100%"))
                      )
                    ),
                    div(class = "action-bar", actionButton(ns("btn_add_pathmox_segvar"), "Add variable", class = "btn-success btn-compact"), actionButton(ns("btn_remove_pathmox_segvar"), "Remove selected", class = "btn-outline-danger btn-compact")),
                    tags$hr(), DTOutput(ns("pathmox_segvar_config_table"))
                  )
                ),
                tags$br(),
                layout_columns(
                  col_widths = c(3, 3, 3, 3),
                  numericInput(ns("pathmox_alpha"), "Alpha", value = 0.05, min = 0.001, max = 0.20, step = 0.001), 
                  numericInput(ns("pathmox_deep"), "Maximum depth", value = 3, min = 1, max = 10, step = 1), 
                  numericInput(ns("pathmox_size"), "Minimum node size", value = 30, min = 5, step = 1), 
                  numericInput(ns("pathmox_size_candidate"), "Minimum candidate split size", value = 15, min = 3, step = 1)
                ),
                layout_columns(
                  col_widths = c(4, 4, 4),
                  checkboxInput(ns("pathmox_test_invariance"), "Run MICOM during splitting", value = TRUE),
                  checkboxInput(ns("pathmox_require_critical_invariance"), "Stop split if critical invariance fails", value = FALSE),
                  numericInput(ns("pathmox_boot_micom"), "Bootstrap resamples (MGA/MICOM)", value = 200, min = 100, step = 100)
                ),
                actionButton(ns("btn_run_pathmox"), "Run PATHMOX", class = "btn-success btn-compact")
              )
            ),
            card(full_screen = TRUE, card_header("PATHMOX segmentation tree"), card_body(DiagrammeR::grVizOutput(ns("pathmox_tree_plot"), height = "760px"))),
            layout_columns(
              col_widths = c(6, 6),
              card(full_screen = TRUE, card_header("Terminal node results"), card_body(DTOutput(ns("pathmox_terminal_paths_table")))),
              card(full_screen = TRUE, card_header("Variable importance"), card_body(DTOutput(ns("pathmox_varimp_table"))))
            ),
            layout_columns(
              col_widths = c(6, 6),
              card(full_screen = TRUE, card_header("Global Compositional invariance (MICOM)"), card_body(DTOutput(ns("pathmox_micom_table")))),
              card(full_screen = TRUE, card_header("Global Henseler MGA results"), card_body(DTOutput(ns("pathmox_mga_table"))))
            ),
            
            card(
              full_screen = TRUE, card_header("Detailed SEM Results by Terminal Node"),
              card_body(
                layout_columns(col_widths = c(4, 8),
                               selectInput(ns("pathmox_detail_group"), "Select terminal node to view details:", choices = NULL),
                               div(class="small-muted", style="padding-top: 35px;", "These results are calculated for the selected terminal group using the global model configuration.")
                ),
                navset_card_pill(
                  nav_panel("Model fit", DTOutput(ns("pmx_dtl_fit"))),
                  nav_panel("R-squared", DTOutput(ns("pmx_dtl_r2"))),
                  nav_panel("Paths", DTOutput(ns("pmx_dtl_paths"))),
                  nav_panel("Loadings", DTOutput(ns("pmx_dtl_loadings"))),
                  nav_panel("Measurement quality", DTOutput(ns("pmx_dtl_quality"))),
                  nav_panel("HTMT", DTOutput(ns("pmx_dtl_htmt"))),
                  nav_panel("Fornell-Larcker", DTOutput(ns("pmx_dtl_fornell")))
                )
              )
            ),
            
            card(card_header("PATHMOX terminal assignment"), card_body(downloadButton(ns("download_pathmox_data"), "Download PATHMOX Data (CSV)", class = "btn-success btn-compact mb-2"), DTOutput(ns("pathmox_individual_table"))))
  )
}

mod_pathmox_server <- function(id, analysis_data_aug_rv, model_lavaan, result_rv, pathmox_segvars_rv, pathmox_results_rv, pathmox_micom_rv, pathmox_mga_rv, pmx_detail_cache, pathmox_levels_available_rv, pathmox_levels_selected_rv, handle_inadmissibles, approach_weights, approach_paths, pls_inner_scheme, plsc_disattenuate, resample_method, n_boot) {
  moduleServer(id, function(input, output, session) {
    
    # ----------------------------------------------------------------------------
    # LOCAL HELPERS PARA TABLAS Y NOMBRES DE NODOS
    # ----------------------------------------------------------------------------
    pathmox_get <- function(x, candidates) {
      if (is.null(x)) return(NULL)
      for (nm in candidates) if (!is.null(x[[nm]])) return(x[[nm]])
      NULL
    }
    
    normalize_node_label <- function(x) {
      paste0("Node ", gsub("^node\\s*", "", trimws(as.character(x)), ignore.case = TRUE))
    }
    
    format_pathmox_comp_micom <- function(x, group_names = NULL) {
      parts <- trimws(strsplit(as.character(x), "_", fixed = TRUE)[[1]])
      parts <- vapply(parts, function(p) {
        p_clean <- gsub("^node\\s*", "", p, ignore.case = TRUE)
        if (!is.null(group_names) && grepl("^[0-9]+$", p_clean)) {
          idx <- as.integer(p_clean)
          if (!is.na(idx) && idx >= 1 && idx <= length(group_names)) return(normalize_node_label(group_names[idx]))
        }
        normalize_node_label(p)
      }, character(1))
      paste(parts, collapse = " versus ")
    }
    
    format_pathmox_comp <- function(x, group_names = NULL) {
      parts <- trimws(strsplit(as.character(x), "_", fixed = TRUE)[[1]])
      map_group <- function(p) {
        p_chr   <- trimws(as.character(p))
        p_clean <- gsub("^node\\s*", "", p_chr, ignore.case = TRUE)
        
        if (!is.null(group_names) && length(group_names) > 0) {
          gn_clean <- gsub("^node\\s*", "", trimws(group_names), ignore.case = TRUE)
          hit_exact <- match(tolower(p_chr), tolower(group_names))
          if (!is.na(hit_exact)) return(normalize_node_label(group_names[hit_exact]))
          hit_clean <- match(tolower(p_clean), tolower(gn_clean))
          if (!is.na(hit_clean)) return(normalize_node_label(group_names[hit_clean]))
        }
        normalize_node_label(p_chr)
      }
      parts <- vapply(parts, map_group, character(1))
      paste(parts, collapse = " versus ")
    }
    
    # ----------------------------------------------------------------------------
    # UI Y ESTADOS DE CONFIGURACIÓN
    # ----------------------------------------------------------------------------
    observe({
      dta <- analysis_data_aug_rv()
      choices <- if (is.null(dta) || ncol(dta) == 0) character(0) else colnames(dta)
      cfg <- isolate(pathmox_segvars_rv())
      if (!is.null(cfg) && nrow(cfg) > 0) pathmox_segvars_rv(cfg[cfg$Variable %in% choices, , drop = FALSE])
      selected <- isolate(input$pathmox_segvar_selected)
      updateSelectInput(session, "pathmox_segvar_selected", choices = choices, selected = if (is.null(selected) || !selected %in% choices) (if (length(choices) > 0) choices[1] else character(0)) else selected)
    })
    
    observeEvent(list(input$pathmox_segvar_selected, input$pathmox_segvar_treatment, analysis_data_aug_rv()), {
      req(analysis_data_aug_rv(), input$pathmox_segvar_selected)
      p <- input$pathmox_segvar_treatment
      v <- input$pathmox_segvar_selected
      if (!(p %in% c("factor", "character"))) {
        pathmox_levels_available_rv(character(0))
        pathmox_levels_selected_rv(character(0))
        updateSelectInput(session, "pathmox_levels_available", choices = character(0))
        updateSelectInput(session, "pathmox_levels_selected", choices = character(0))
        return(NULL)
      }
      vals <- sort(unique(trimws(as.character(analysis_data_aug_rv()[[v]]))[!is.na(trimws(as.character(analysis_data_aug_rv()[[v]]))) & nzchar(trimws(as.character(analysis_data_aug_rv()[[v]]))) ]))
      pathmox_levels_available_rv(vals)
      pathmox_levels_selected_rv(character(0))
      updateSelectInput(session, "pathmox_levels_available", choices = vals)
      updateSelectInput(session, "pathmox_levels_selected", choices = character(0))
    }, ignoreInit = FALSE)
    
    observeEvent(input$btn_pathmox_level_right, {
      sel <- input$pathmox_levels_available; if (length(sel) == 0) return(NULL)
      pathmox_levels_selected_rv(c(pathmox_levels_selected_rv(), sel))
      pathmox_levels_available_rv(setdiff(pathmox_levels_available_rv(), sel))
      updateSelectInput(session, "pathmox_levels_available", choices = pathmox_levels_available_rv())
      updateSelectInput(session, "pathmox_levels_selected", choices = pathmox_levels_selected_rv())
    })
    
    observeEvent(input$btn_pathmox_level_left, {
      sel <- input$pathmox_levels_selected; if (length(sel) == 0) return(NULL)
      pathmox_levels_selected_rv(setdiff(pathmox_levels_selected_rv(), sel))
      pathmox_levels_available_rv(sort(unique(c(pathmox_levels_available_rv(), sel))))
      updateSelectInput(session, "pathmox_levels_available", choices = pathmox_levels_available_rv())
      updateSelectInput(session, "pathmox_levels_selected", choices = pathmox_levels_selected_rv())
    })
    
    observeEvent(input$btn_pathmox_level_up, {
      sel <- input$pathmox_levels_selected; cur <- pathmox_levels_selected_rv(); if (length(sel) != 1 || length(cur) <= 1) return(NULL)
      idx <- match(sel, cur); if (is.na(idx) || idx <= 1) return(NULL)
      tmp <- cur[idx - 1]; cur[idx - 1] <- cur[idx]; cur[idx] <- tmp
      pathmox_levels_selected_rv(cur); updateSelectInput(session, "pathmox_levels_selected", choices = cur, selected = sel)
    })
    
    observeEvent(input$btn_pathmox_level_down, {
      sel <- input$pathmox_levels_selected; cur <- pathmox_levels_selected_rv(); if (length(sel) != 1 || length(cur) <= 1) return(NULL)
      idx <- match(sel, cur); if (is.na(idx) || idx >= length(cur)) return(NULL)
      tmp <- cur[idx + 1]; cur[idx + 1] <- cur[idx]; cur[idx] <- tmp
      pathmox_levels_selected_rv(cur); updateSelectInput(session, "pathmox_levels_selected", choices = cur, selected = sel)
    })
    
    output$pathmox_segvar_config_table <- renderDT({
      cfg <- pathmox_segvars_rv()
      if (is.null(cfg) || nrow(cfg) == 0) return(empty_dt("No segmentation variables configured yet."))
      df <- cfg
      if ("Ordered" %in% names(df)) df$Ordered <- ifelse(df$Ordered, "Yes", "No")
      names(df) <- c("Variable", "Processing", "Bins", "Levels", "Ordered")
      DT::datatable(df, selection = "single", rownames = FALSE, options = list(dom = "tip", pageLength = 10, scrollX = TRUE))
    })
    
    observeEvent(input$btn_add_pathmox_segvar, {
      req(input$pathmox_segvar_selected)
      cfg <- pathmox_segvars_rv()
      v <- input$pathmox_segvar_selected
      p <- input$pathmox_segvar_treatment
      b <- if (p %in% c("quantile", "equal")) input$pathmox_segvar_bins else NA_real_
      ord <- if (p %in% c("quantile", "equal")) TRUE else if (p %in% c("factor", "character")) isTRUE(input$pathmox_segvar_ordered) else FALSE
      levs <- if (p %in% c("factor", "character") && ord) paste(pathmox_levels_selected_rv(), collapse = ", ") else ""
      
      if (p %in% c("factor", "character") && ord && length(pathmox_levels_selected_rv()) == 0) { showNotification("Select and order at least one level.", type = "error"); return(NULL) }
      
      new_row <- data.frame(Variable = v, Processing = p, Bins = b, Levels = levs, Ordered = ord, stringsAsFactors = FALSE)
      if (is.null(cfg) || nrow(cfg) == 0) cfg <- new_row else if (v %in% cfg$Variable) cfg[cfg$Variable == v, c("Processing", "Bins", "Levels", "Ordered")] <- new_row[1, c("Processing", "Bins", "Levels", "Ordered")] else cfg <- rbind(cfg, new_row)
      pathmox_segvars_rv(cfg)
    })
    
    observeEvent(input$btn_remove_pathmox_segvar, {
      idx <- input$pathmox_segvar_config_table_rows_selected
      cfg <- pathmox_segvars_rv()
      if (is.null(idx) || length(idx) == 0 || is.null(cfg) || nrow(cfg) == 0) { showNotification("Select a segmentation variable to remove.", type = "error"); return(NULL) }
      pathmox_segvars_rv(cfg[-idx, , drop = FALSE])
    })
    
    # ----------------------------------------------------------------------------
    # EJECUCIÓN PRINCIPAL DE PATHMOX (CON ALGORITMO CUSTOM Y CAPTURA DE ERRORES)
    # ----------------------------------------------------------------------------
    observeEvent(input$btn_run_pathmox, {
      req(analysis_data_aug_rv(), model_lavaan(), result_rv())
      dta <- analysis_data_aug_rv()
      model <- model_lavaan()
      res <- result_rv()
      cfg <- pathmox_segvars_rv()
      
      if (identical(model, "No constructs defined yet.")) { showNotification("Define constructs first.", type = "error"); return(NULL) }
      if (is.null(cfg) || nrow(cfg) == 0) { showNotification("Configure at least one segmentation variable.", type = "error"); return(NULL) }
      
      cat_data <- try(prepare_pathmox_segvars(dta, cfg), silent = TRUE)
      if (inherits(cat_data, "try-error") || is.null(cat_data)) { 
        showNotification(paste0("Segmentation preprocessing failed. ", gsub("Error in .*?:\\s*", "", as.character(cat_data))), type = "error", duration = 8)
        return(NULL) 
      }
      if (nrow(cat_data) == 0 || anyNA(cat_data)) { showNotification("Segmentation data invalid or missing values present.", type = "error"); return(NULL) }
      
      pathmox_input <- as.data.frame(dta)
      pathmox_input$id_case_pathmox <- seq_len(nrow(pathmox_input))
      valid_seg <- vapply(cat_data, function(x) length(unique(stats::na.omit(x))) >= 2, logical(1))
      cat_data <- cat_data[, valid_seg, drop = FALSE]
      if (ncol(cat_data) == 0) { showNotification("No segmentation variable has at least two valid categories.", type = "error"); return(NULL) }
      
      pathmox_results_rv(NULL); pathmox_micom_rv(NULL); pathmox_mga_rv(NULL); pmx_detail_cache(list())
      
      withProgress(message = "Running Custom PATHMOX+MICOM...", value = 0, {
        tryCatch({
          incProgress(0.20, detail = "Splitting nodes recursively")
          
          # LLAMADA A LA FUNCIÓN CUSTOM
          pmx <- try(pathmoxMICOM(
            .model = model,
            .data = pathmox_input,
            .catvar = cat_data,
            .scheme = pls_inner_scheme(),
            .consistent = plsc_disattenuate(),
            .alpha = input$pathmox_alpha,
            .deep = input$pathmox_deep,
            .size = input$pathmox_size,
            .size_candidate = input$pathmox_size_candidate,
            .tree = FALSE,
            .test_invariance = input$pathmox_test_invariance,
            .alpha_invariance = input$pathmox_alpha,
            .require_critical_invariance = input$pathmox_require_critical_invariance,
            .boot_micom = input$pathmox_boot_micom,
            .verbose = FALSE
          ), silent = TRUE)
          
          if (inherits(pmx, "try-error")) stop(gsub("Error in .*?:\\s*", "", as.character(pmx)))
          
          incProgress(0.75, detail = "Running Global Bootstrapping for MGA/MICOM")
          mga_res <- NULL; csem_by_group <- NULL; micom_res <- NULL
          
          # Evalución global de nodos terminales
          if (!is.null(pmx$hybrid) && is.list(pmx$hybrid) && length(pmx$hybrid) >= 2) {
            
            # MODELO BASE SIN RESAMPLE (Para evitar doble bootstrap que rompe la matriz)
            csem_by_group <- tryCatch({
              cSEM::csem(
                .data = pmx$hybrid, 
                .model = model, 
                .resample_method = "none", # <-- CAMBIO CLAVE
                .handle_inadmissibles = handle_inadmissibles(), 
                .approach_weights = approach_weights(), 
                .approach_paths = approach_paths(), 
                .PLS_weight_scheme_inner = pls_inner_scheme(), 
                .disattenuate = plsc_disattenuate()
              )
            }, error = function(e) e)
            
            if (inherits(csem_by_group, "error") || inherits(csem_by_group, "try-error")) {
              micom_res <- list(error = paste("Base model failed:", conditionMessage(csem_by_group)))
              mga_res <- list(error = paste("Base model failed:", conditionMessage(csem_by_group)))
            } else {
              # MICOM GLOBAL
              micom_res <- tryCatch({
                cSEM::testMICOM(csem_by_group, .R = max(50, input$pathmox_boot_micom), .seed = 12345,.approach_p_adjust = "none", .verbose = FALSE)
              }, error = function(e) list(error = conditionMessage(e)))
              
              # MGA GLOBAL
              mga_res <- tryCatch({
                model_lines <- trimws(strsplit(model, "\n")[[1]])
                structural_model <- paste(model_lines[grepl("~", model_lines, fixed = TRUE) & !grepl("=~", model_lines, fixed = TRUE) & !grepl("<~", model_lines, fixed = TRUE)], collapse = "\n")
                if (nzchar(structural_model)) {
                  get("testMGD", envir = asNamespace("cSEM"))(csem_by_group, .parameters_to_compare = structural_model, .approach_mgd = "Henseler", .R_bootstrap = max(50, input$pathmox_boot_micom), .verbose = FALSE)
                } else {
                  list(error = "No structural paths to compare.")
                }
              }, error = function(e) list(error = conditionMessage(e)))
            }
          } else {
            micom_res <- list(error = "Not enough terminal nodes (>1) to run global tests.")
            mga_res <- list(error = "Not enough terminal nodes (>1) to run global tests.")
          }
          
          incProgress(0.85, detail = "Calculating detailed SEM results for all nodes")
          detail_cache_all <- list()
          
          if (!is.null(pmx$hybrid) && is.list(pmx$hybrid) && length(pmx$hybrid) > 0) {
            res_glob <- isolate(try(result_rv(), silent = TRUE))
            
            if (!inherits(res_glob, "try-error") && !is.null(res_glob)) {
              base_args <- res_glob$Information$Arguments
              
              if (!is.null(base_args)) {
                for (node in names(pmx$hybrid)) {
                  dat_node <- pmx$hybrid[[node]]
                  if (is.null(dat_node)) next
                  
                  res_node <- tryCatch({
                    args <- base_args
                    args$.data <- dat_node
                    args$.resample_method <- resample_method()
                    args$.R <- if (resample_method() == "none") 0 else n_boot()
                    args$.handle_inadmissibles <- handle_inadmissibles()
                    do.call(cSEM::csem, args)
                  }, error = function(e) NULL)
                  
                  if (!is.null(res_node)) {
                    mm <- extract_mm_quality(res_node)
                    detail_cache_all[[node]] <- list(
                      fit = extract_model_fit(res_node),
                      types = get_construct_types(res_node),
                      r2 = extract_r2(res_node),
                      paths = extract_paths(res_node),
                      loadings = extract_loadings(res_node),
                      quality = mm$quality,
                      htmt = mm$htmt,
                      fornell = mm$fornell
                    )
                  }
                }
              }
            }
          }
          pmx_detail_cache(detail_cache_all)
      
          pathmox_results_rv(pmx)
          pathmox_micom_rv(micom_res)
          pathmox_mga_rv(mga_res)
          incProgress(1.0, detail = "Complete")
          showNotification("PATHMOX+MICOM and Global MGA completed successfully.", type = "message")
          
        }, error = function(e) {
          pathmox_results_rv(NULL); pathmox_micom_rv(NULL); pathmox_mga_rv(NULL); pmx_detail_cache(list())
          showNotification(paste0("PATHMOX failed: ", e$message), type = "error", duration = 8) 
        })
      })
    })
    
    # ----------------------------------------------------------------------------
    # OUTPUTS: ÁRBOL Y TABLAS GENERALES PATHMOX
    # ----------------------------------------------------------------------------
    output$pathmox_tree_plot <- DiagrammeR::renderGrViz({
      pmx <- pathmox_results_rv()
      req(pmx)
      if (is.null(pmx$MOX)) return(DiagrammeR::grViz("digraph { a [label='No PATHMOX results available', shape=box] }"))
      mox <- as.data.frame(pmx$MOX, stringsAsFactors = FALSE)
      if (nrow(mox) == 0 || !all(c("node", "parent") %in% names(mox))) return(DiagrammeR::grViz("digraph { a [label='Invalid tree structure', shape=box] }"))
      
      cfg <- pathmox_segvars_rv()
      catvar_map <- if (!is.null(cfg) && nrow(cfg) > 0) setNames(cfg$Variable, rep(".catvar", nrow(cfg))) else NULL
      fix_catvar_label <- function(x) { x <- as.character(x); if (!is.null(catvar_map)) for (i in seq_along(x)) if (x[i] %in% names(catvar_map)) x[i] <- catvar_map[[x[i]]]; x }
      
      esc <- function(x) { 
        x <- ifelse(is.na(x), "", as.character(x))
        gsub("\"", "&quot;", gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE), fixed = TRUE), fixed = TRUE) 
      }
      
      pct_col <- intersect(c("%", "percent", "pct"), names(mox))
      pct_col <- if (length(pct_col) > 0) pct_col[1] else NA_character_
      
      node_label <- vapply(seq_len(nrow(mox)), function(i) {
        if (tolower(as.character(mox$terminal[i])) %in% c("yes", "true", "terminal")) {
          pct_txt <- if (!is.na(pct_col)) paste0(" (", round(suppressWarnings(as.numeric(mox[[pct_col]][i])), 1), "%)") else ""
          paste0("<B>Node ", esc(mox$node[i]), "</B><BR/>", esc(mox$size[i]), pct_txt)
        } else {
          child_rows <- mox[mox$parent == mox$node[i], , drop = FALSE]
          split_var <- if (nrow(child_rows) > 0 && "variable" %in% names(child_rows)) as.character(unique(fix_catvar_label(as.character(child_rows$variable)[nzchar(as.character(child_rows$variable))]))[1]) else "Split"
          paste0("<B>", esc(split_var), "</B>")
        }
      }, character(1))
      
      node_style <- vapply(seq_len(nrow(mox)), function(i) {
        if (mox$parent[i] == "0" || mox$parent[i] == 0) 'shape=box style="rounded,filled" fillcolor="#DCEEFF" color="#4A90C2" penwidth=2.6 fontname="Helvetica" fontsize=13 margin="0.18,0.12"'
        else if (tolower(as.character(mox$terminal[i])) %in% c("yes", "true", "terminal")) 'shape=box style="rounded,filled" fillcolor="#DFF3E4" color="#4C9A67" penwidth=2.0 fontname="Helvetica" fontsize=12 margin="0.14,0.10"'
        else 'shape=box style="rounded,filled" fillcolor="#FFF4CC" color="#D4A017" penwidth=2.0 fontname="Helvetica" fontsize=12 margin="0.16,0.11"'
      }, character(1))
      
      node_lines <- paste0("n", mox$node, ' [label=<', node_label, '>, ', node_style, ']')
      edge_rows <- mox[mox$parent != 0 & !is.na(mox$parent), , drop = FALSE]
      edge_lines <- if (nrow(edge_rows) > 0) vapply(seq_len(nrow(edge_rows)), function(i) {
        paste0("n", edge_rows$parent[i], " -> n", edge_rows$node[i], ' [label="', esc(gsub("/", " / ", ifelse(is.na(edge_rows$category[i]), "", as.character(edge_rows$category[i])), fixed = TRUE)), '", color="', ifelse(tolower(as.character(edge_rows$terminal[i])) %in% c("yes", "true", "terminal"), "#4C9A67", "#C49A00"), '", fontcolor="#4F4F4F", fontsize=11, penwidth=1.8, arrowsize=0.9]')
      }, character(1)) else character(0)
      
      DiagrammeR::grViz(paste0("digraph pathmox_tree { graph [layout=dot, rankdir=TB, bgcolor='white', nodesep=0.55, ranksep=0.95, splines=curved, pad=0.25]; node [fontname='Helvetica']; edge [fontname='Helvetica', labelfloat=false, labeldistance=1.8, labelangle=0]; ", paste(node_lines, collapse = " "), " ", paste(edge_lines, collapse = " "), "}"))
    })
    
    output$pathmox_terminal_paths_table <- renderDT({
      pmx <- pathmox_results_rv()
      df <- pathmox_get(pmx, c("terminal_paths", "terminalpaths"))
      if (is.null(df)) return(empty_dt("Run PATHMOX analysis to see terminal model results."))
      df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
      if (!is.null(rownames(df)) && any(nzchar(rownames(df)))) {
        df <- data.frame(Parameter = rownames(df), df, check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)
      } else {
        names(df)[1] <- "Parameter"
      }
      df$Parameter <- gsub("->", " -> ", df$Parameter)
      df$Parameter <- gsub("\\s+", " ", df$Parameter)
      df$Parameter <- trimws(df$Parameter)
      names(df) <- gsub("^node\\s*", "Node ", names(df), ignore.case = TRUE)
      tbl <- DT::datatable(df, rownames = FALSE, options = list(pageLength = 20, scrollX = TRUE, dom = "tip"))
      num_cols <- names(df)[-1]
      if (length(num_cols) > 0) tbl <- DT::formatRound(tbl, columns = num_cols, digits = 3)
      tbl
    })
    
    output$pathmox_varimp_table <- renderDT({
      pmx <- pathmox_results_rv()
      df <- pathmox_get(pmx, c("varimp", "var_imp", "variable_importance"))
      if (is.null(df)) return(empty_dt("Variable importance is not available."))
      
      df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
      cfg <- pathmox_segvars_rv()
      if (!is.null(cfg) && nrow(cfg) > 0) {
        catvar_map <- setNames(cfg$Variable, rep(".catvar", nrow(cfg)))
        for (j in seq_along(df)) if (is.character(df[[j]]) || is.factor(df[[j]])) { df[[j]] <- as.character(df[[j]]); for (old in names(catvar_map)) df[[j]][df[[j]] == old] <- catvar_map[[old]] }
        if (!is.null(rownames(df))) for (old in names(catvar_map)) rownames(df)[rownames(df) == old] <- catvar_map[[old]]
        for (old in names(catvar_map)) names(df)[names(df) == old] <- catvar_map[[old]]
      }
      names(df) <- tools::toTitleCase(gsub("_", " ", names(df)))
      standard_dt(df, digits = 3)
    })
    
    output$pathmox_micom_table <- renderDT({
      micom <- pathmox_micom_rv()
      pmx <- pathmox_results_rv()
      if (is.null(micom)) return(empty_dt("MICOM results are not available."))
      if (!is.null(micom$error)) return(empty_dt(paste("MICOM Failed:", micom$error)))
      if (is.null(micom$Step2) || is.null(micom$Step2$P_value)) return(empty_dt("MICOM Step 2 p-values not found."))
      
      p_val_obj <- micom$Step2$P_value
      pvals <- if (is.list(p_val_obj) && "none" %in% names(p_val_obj)) p_val_obj$none else p_val_obj
      
      # Forzar a formato de lista para comparaciones de 2 grupos
      if (!is.list(pvals)) {
        comp_name <- "Comparison"
        if (!is.null(pmx$hybrid) && length(pmx$hybrid)>=2) comp_name <- paste(names(pmx$hybrid)[1:2], collapse="_")
        pvals <- stats::setNames(list(pvals), comp_name)
      }
      
      alpha <- input$pathmox_alpha %||% 0.05
      group_names <- if (!is.null(pmx) && !is.null(pmx$hybrid)) names(pmx$hybrid) else NULL
      rows_list <- list()
      for (comp in names(pvals)) {
        cpvals <- pvals[[comp]]
        cnames <- names(cpvals) %||% paste0("Construct_", seq_along(cpvals))
        for (i in seq_along(cpvals)) {
          pval <- as.numeric(cpvals[i])
          rows_list[[length(rows_list) + 1]] <- data.frame(Comparison = format_pathmox_comp_micom(comp, group_names), Construct = cnames[i], `p-value` = pval, Sig = ifelse(!is.na(pval) && pval >= alpha, "Y", "N"), stringsAsFactors = FALSE, check.names = FALSE)
        }
      }
      if (length(rows_list) == 0) return(empty_dt("No MICOM comparisons found."))
      tbl <- DT::datatable(do.call(rbind, rows_list), rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE, autoWidth = TRUE, dom = "tip", columnDefs = list(list(width = "280px", targets = 0), list(width = "180px", targets = 1), list(visible = FALSE, targets = 3))))
      DT::formatStyle(DT::formatRound(tbl, columns = "p-value", digits = 3), columns = "p-value", valueColumns = "Sig", color = DT::styleEqual(c("Y", "N"), c("#198754", "#DC3545")))
    })
    
    output$pathmox_mga_table <- renderDT({
      mga <- pathmox_mga_rv()
      pmx <- pathmox_results_rv()
      if (is.null(mga)) return(empty_dt("Henseler MGA results are not available."))
      if (!is.null(mga$error)) return(empty_dt(paste("MGA Failed:", mga$error)))
      if (is.null(mga$Henseler) || is.null(mga$Henseler$P_value)) return(empty_dt("Henseler MGA P-values not found."))
      
      p_val_obj <- mga$Henseler$P_value
      pvals <- if (is.list(p_val_obj) && "none" %in% names(p_val_obj)) p_val_obj$none else p_val_obj
      
      if (!is.list(pvals)) {
        comp_name <- "Comparison"
        if (!is.null(pmx$hybrid) && length(pmx$hybrid)>=2) comp_name <- paste(names(pmx$hybrid)[1:2], collapse="_")
        pvals <- stats::setNames(list(pvals), comp_name)
      }
      
      group_names <- if (!is.null(pmx) && !is.null(pmx$hybrid)) names(pmx$hybrid) else NULL
      fmt_param <- function(x) { if (grepl("~", x, fixed = TRUE)) { p <- trimws(strsplit(x, "~", fixed = TRUE)[[1]]); if (length(p) == 2) return(paste0(p[2], " -> ", p[1])) }; x }
      rows_list <- list()
      for (comp in names(pvals)) {
        for (i in seq_along(pvals[[comp]])) {
          pval <- as.numeric(pvals[[comp]][i])
          rows_list[[length(rows_list) + 1]] <- data.frame(Comparison = format_pathmox_comp(comp, group_names), Parameter = fmt_param(names(pvals[[comp]])[i]), `p-value` = pval, Sig = ifelse(!is.na(pval) && (pval < 0.05 || pval > 0.95), "Y", "N"), stringsAsFactors = FALSE, check.names = FALSE)
        }
      }
      if (length(rows_list) == 0) return(empty_dt("No Henseler MGA comparisons found."))
      tbl <- DT::datatable(do.call(rbind, rows_list), rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE, autoWidth = TRUE, dom = "tip", columnDefs = list(list(width = "280px", targets = 0), list(width = "180px", targets = 1), list(width = "90px", targets = 2), list(visible = FALSE, targets = 3))))
      DT::formatStyle(DT::formatRound(tbl, columns = "p-value", digits = 3), columns = "p-value", valueColumns = "Sig", color = DT::styleEqual(c("Y", "N"), c("#198754", "#DC3545")))
    })
    
    # ----------------------------------------------------------------------------
    # DETALLES CSEM POR NODO (Inner Evaluation)
    # ----------------------------------------------------------------------------
    observe({
      pmx <- pathmox_results_rv()
      if (!is.null(pmx) && !is.null(pmx$hybrid)) {
        nms <- names(pmx$hybrid)
        opciones <- setNames(nms, gsub("^node\\s*", "Node ", nms, ignore.case = TRUE))
        updateSelectInput(session, "pathmox_detail_group", choices = opciones)
      } else {
        updateSelectInput(session, "pathmox_detail_group", choices = character(0))
      }
    })
    
    pmx_selected_data <- reactive({
      req(input$pathmox_detail_group)
      cache <- pmx_detail_cache()
      # Simplemente retorna lo que ya calculamos en el bloque general de PATHMOX
      return(cache[[input$pathmox_detail_group]])
    })
    
    output$pmx_dtl_fit <- renderDT({ d <- pmx_selected_data(); if (is.null(d) || is.null(d$fit) || nrow(d$fit) == 0) empty_dt("Not available.") else standard_dt(d$fit, digits = 3) })
    output$pmx_dtl_r2 <- renderDT({ d <- pmx_selected_data(); if (is.null(d) || is.null(d$r2) || nrow(d$r2) == 0) empty_dt("Not available.") else standard_dt(d$r2, digits = 3) })
    output$pmx_dtl_paths <- renderDT({ d <- pmx_selected_data(); if (is.null(d) || is.null(d$paths) || nrow(d$paths) == 0) empty_dt("Not available.") else standard_dt(d$paths, digits = 3) })
    output$pmx_dtl_loadings <- renderDT({ d <- pmx_selected_data(); if (is.null(d) || is.null(d$loadings) || nrow(d$loadings) == 0) empty_dt("Not available.") else standard_dt(d$loadings, digits = 3) })
    output$pmx_dtl_quality <- renderDT({ d <- pmx_selected_data(); if (is.null(d) || is.null(d$quality) || nrow(d$quality) == 0) empty_dt("Not available.") else standard_dt(d$quality, digits = 3) })
    output$pmx_dtl_htmt <- renderDT({ d <- pmx_selected_data(); if (is.null(d) || is.null(d$htmt) || nrow(d$htmt) == 0) empty_dt("Not available.") else standard_dt(d$htmt, digits = 3) })
    output$pmx_dtl_fornell <- renderDT({ d <- pmx_selected_data(); if (is.null(d) || is.null(d$fornell) || nrow(d$fornell) == 0) empty_dt("Not available.") else standard_dt(d$fornell, digits = 3) })
    
    output$pathmox_individual_table <- renderDT({
      pmx <- pathmox_results_rv()
      if (is.null(pmx) || is.null(pmx$hybrid)) return(empty_dt("Run PATHMOX analysis to see terminal assignment."))
      pieces <- Filter(Negate(is.null), lapply(names(pmx$hybrid), function(nm) { df <- pmx$hybrid[[nm]]; if (is.null(df) || nrow(df) == 0) return(NULL); df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE); if (!"hybrid_var" %in% names(df)) df$hybrid_var <- nm; df }))
      if (length(pieces) == 0) return(empty_dt("No terminal node data available."))
      df <- do.call(rbind, pieces)
      names(df)[names(df) == "hybrid_var"] <- "Terminal Node"
      standard_dt(df[, c(which(names(df) == "Terminal Node"), setdiff(seq_along(df), which(names(df) == "Terminal Node"))), drop = FALSE], digits = 3)
    })
    
    output$download_pathmox_data <- downloadHandler(
      filename = function() paste0("pathmox_terminal_data_", Sys.Date(), ".csv"),
      content = function(file) {
        pmx <- pathmox_results_rv(); req(pmx, pmx$hybrid)
        pieces <- Filter(Negate(is.null), lapply(names(pmx$hybrid), function(nm) { df <- as.data.frame(pmx$hybrid[[nm]], stringsAsFactors = FALSE); df$Terminal_Node <- nm; df }))
        write.csv(do.call(rbind, pieces), file, row.names = FALSE, fileEncoding = "UTF-8")
      }
    )
    
  })
}