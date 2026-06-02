# ==============================================================================
# MODULE: WPI ANALYSIS
# ==============================================================================

mod_wpi_ui <- function(id) {
  ns <- NS(id)
  
  nav_panel("WPI analysis",
            card(
              card_header("Weighted Performance Index (WPI) Configuration"),
              card_body(
                selectInput(ns("wpi_target_construct"), "Select target construct:", choices = NULL),
                div(
                  class = "small-muted mb-3",
                  "WPI measures individual performance weighted by total effects from PLS-SEM. Higher WPI = better performance on important predictors."
                ),
                actionButton(ns("btn_run_wpi"), "Calculate WPI & Run t-SNE", class = "btn-success btn-compact mb-2")
              )
            ),
            
            card(
              full_screen = TRUE,
              card_header("t-SNE Projection with WPI"),
              card_body(
                plotOutput(ns("wpi_tsne_plot"), height = "600px")
              )
            ),
            
            card(
              card_header("WPI Segment Table"),
              card_body(
                DTOutput(ns("wpi_segment_table"))
              )
            ),
            
            card(
              full_screen = TRUE,
              card_header("WPI Distribution & Validation"),
              card_body(
                layout_columns(
                  col_widths = c(6, 6),
                  plotOutput(ns("wpi_hist_plot"), height = "300px"),
                  plotOutput(ns("wpi_scatter_plot"), height = "300px")
                )
              )
            ),
            
            card(
              card_header("Individual WPI Data"),
              card_body(
                downloadButton(ns("download_wpi_data"), "Download WPI Data (CSV)", class = "btn-success btn-compact mb-2"),
                DTOutput(ns("wpi_individual_table"))
              )
            )
  )
}

mod_wpi_server <- function(id, result_rv, analysis_data, constructs_rv, wpi_results_rv, analysis_data_aug_rv, loading_project_rv) {
  moduleServer(id, function(input, output, session) {
    
    # ----------------------------------------------------------------------------
    # UI DINÁMICA: ACTUALIZAR CONSTRUCTO OBJETIVO
    # ----------------------------------------------------------------------------
    observe({
      clist <- constructs_rv()
      cnames <- names(clist)
      if (is.null(cnames) || length(cnames) == 0) {
        updateSelectInput(session, "wpi_target_construct", choices = character(0), selected = character(0))
      } else {
        # Si existe input$depvar (de mod_relations), lo ideal sería enlazarlo, 
        # pero como fallback tomamos el primer constructo.
        default_target <- cnames[1] 
        updateSelectInput(session, "wpi_target_construct", choices = cnames, selected = default_target)
      }
    })
    
    # ----------------------------------------------------------------------------
    # EJECUCIÓN DEL ANÁLISIS WPI
    # ----------------------------------------------------------------------------
    observeEvent(input$btn_run_wpi, {
      req(input$wpi_target_construct)
      if (is.null(result_rv())) {
        showNotification("Please estimate the SEM model before running WPI analysis.", type = "error")
        return(NULL)
      }
      
      withProgress(message = "Calculating WPI...", detail = "Extracting effects", value = 0, {
        tryCatch({
          res <- result_rv()
          dta <- analysis_data()
          clist <- constructs_rv()
          target <- input$wpi_target_construct
          
          incProgress(0.2, detail = "Extracting total effects")
          assessment <- cSEM::assess(res)
          total_effects <- assessment$Effects$Total_effect
          effects_to_target <- total_effects[grep(paste0(target, " ~"), total_effects$Name), ]
          
          if (nrow(effects_to_target) == 0) {
            showNotification("No predictors found for the selected target construct.", type = "error")
            return(NULL)
          }
          
          predictors <- sub(paste0(target, " ~ "), "", effects_to_target$Name)
          effect_values <- effects_to_target$Estimate
          names(effect_values) <- predictors
          
          header_labels <- ifelse(effect_values < 0, paste0("(-) ", predictors), predictors)
          names(header_labels) <- predictors
          weights <- abs(effect_values) / sum(abs(effect_values))
          
          incProgress(0.3, detail = "Calculating performance scores")
          scores <- as.data.frame(res$Estimates$Construct_scores)
          scores_rescaled <- as.data.frame(lapply(scores, function(x) {
            (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * 100
          }))
          performance <- as.matrix(scores_rescaled[, predictors, drop = FALSE])
          
          calc_performance <- performance
          for (pred in predictors) {
            if (effect_values[pred] < 0) calc_performance[, pred] <- 100 - calc_performance[, pred]
          }
          wpi <- as.vector(calc_performance %*% weights)
          
          incProgress(0.5, detail = "Running t-SNE & K-Means")
          n <- nrow(performance)
          perplexity_auto <- min(50, max(5, floor((n - 1) / 3)))
          set.seed(42)
          tsne_result <- Rtsne(
            performance,
            dims = 2,
            perplexity = perplexity_auto,
            max_iter = 1000,
            check_duplicates = FALSE,
            verbose = FALSE
          )
          
          k_max <- min(10, floor(sqrt(n / 2)))
          wss <- sapply(1:k_max, function(k) kmeans(performance, centers = k, nstart = 10)$tot.withinss)
          
          if (k_max > 2) {
            x1 <- 1; y1 <- wss[1]; x2 <- k_max; y2 <- wss[k_max]
            distances <- sapply(1:k_max, function(k) {
              x0 <- k; y0 <- wss[k]
              abs((x2 - x1) * (y1 - y0) - (x1 - x0) * (y2 - y1)) / sqrt((x2 - x1)^2 + (y2 - y1)^2)
            })
            k_opt <- which.max(distances)
          } else {
            k_opt <- max(1, k_max)
          }
          
          set.seed(42)
          km_res <- kmeans(performance, centers = k_opt, nstart = 25)
          
          incProgress(0.8, detail = "Preparing results")
          wpi_data <- data.frame(
            ID = 1:n,
            tsne1 = tsne_result$Y[, 1],
            tsne2 = tsne_result$Y[, 2],
            WPI = wpi,
            Target_Score = scores_rescaled[[target]],
            Segment = as.factor(km_res$cluster)
          )
          wpi_data <- cbind(wpi_data, as.data.frame(performance))
          
          wpi_summary <- data.frame(
            Predictor = predictors,
            Total_Effect = round(effect_values, 3),
            Weight = round(weights, 3),
            Mean_Performance = round(colMeans(performance), 3),
            SD_Performance = round(apply(performance, 2, sd), 3),
            stringsAsFactors = FALSE
          )
          
          wpi_results_rv(list(
            data = wpi_data,
            summary = wpi_summary,
            target = target,
            n = n,
            perplexity = perplexity_auto,
            segments = k_opt,
            effect_values = effect_values,
            header_labels = header_labels
          ))
          
          if (!is.null(dta) && nrow(dta) == nrow(wpi_data)) {
            dta_aug <- as.data.frame(dta, stringsAsFactors = FALSE)
            dta_aug$WPI_Segment <- factor(wpi_data$Segment)
            analysis_data_aug_rv(dta_aug)
          }
          
          incProgress(1.0, detail = "Complete")
          showNotification("WPI analysis completed successfully!", type = "message")
          
        }, error = function(e) {
          showNotification(paste0("WPI calculation failed: ", e$message), type = "error", duration = 8)
          wpi_results_rv(NULL)
          
          # Rescate: Solo si no se está cargando un proyecto
          if (!isTRUE(loading_project_rv())) {
            analysis_data_aug_rv(analysis_data())
          }
        })
      })
    })
    
    # ----------------------------------------------------------------------------
    # OUTPUTS WPI: GRÁFICOS Y TABLAS
    # ----------------------------------------------------------------------------
    output$wpi_tsne_plot <- renderPlot({
      req(wpi_results_rv())
      wpi_res <- wpi_results_rv(); plot_data <- wpi_res$data
      ggplot(plot_data, aes(x = tsne1, y = tsne2)) +
        geom_mark_hull(aes(group = Segment, label = paste("Seg.", Segment)), fill = NA, color = "gray30", linetype = "dashed", concavity = 5, show.legend = FALSE) +
        geom_point(aes(color = WPI), size = 2.5, alpha = 0.9) +
        scale_color_viridis_c(option = "viridis", name = "WPI Score") +
        labs(title = "t-SNE Projection & Segment Profiling", subtitle = paste0("Target: ", wpi_res$target, " | n = ", wpi_res$n, " | Auto-segments = ", wpi_res$segments), caption = "Grouped by K-Means", x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
        theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 13), legend.position = "right", panel.grid.minor = element_blank(), panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5))
    })
    
    output$wpi_hist_plot <- renderPlot({
      req(wpi_results_rv())
      plot_data <- wpi_results_rv()$data
      ggplot(plot_data, aes(x = WPI, fill = after_stat(x))) +
        geom_histogram(bins = 30, color = "white", linewidth = 0.3) +
        scale_fill_viridis_c(option = "viridis", guide = "none") +
        geom_vline(xintercept = mean(plot_data$WPI), linetype = "dashed", color = "black", linewidth = 1) +
        annotate("text", x = mean(plot_data$WPI), y = Inf, label = paste0("M = ", round(mean(plot_data$WPI), 1), " | SD = ", round(sd(plot_data$WPI), 1)), vjust = 1.8, hjust = -0.05, color = "black", size = 3.8, fontface = "bold") +
        labs(title = "WPI Distribution", x = "Weighted Performance Index", y = "Frequency") + theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
    })
    
    output$wpi_scatter_plot <- renderPlot({
      req(wpi_results_rv())
      wpi_res <- wpi_results_rv(); plot_data <- wpi_res$data
      cor_test <- cor.test(plot_data$WPI, plot_data$Target_Score)
      label_text <- paste0("r = ", round(cor_test$estimate, 3), " | ", ifelse(cor_test$p.value < 0.001, "p < 0.001", paste0("p = ", round(cor_test$p.value, 3))))
      
      ggplot(plot_data, aes(x = WPI, y = Target_Score)) +
        geom_point(aes(color = WPI), size = 2, alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "gray20", fill = "gray80", alpha = 0.3, linewidth = 1.2) +
        scale_color_viridis_c(option = "viridis", guide = "none") +
        annotate("text", x = min(plot_data$WPI) + 3, y = max(plot_data$Target_Score) - 3, label = label_text, hjust = 0, vjust = 1, size = 4, fontface = "bold", color = "gray20") +
        labs(title = paste0("WPI vs ", wpi_res$target), x = "Weighted Performance Index", y = paste0(wpi_res$target, " Score")) + theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
    })
    
    # Originalmente wpi_component_table no tenía UI asociada visible, pero lo mantenemos por consistencia interna
    output$wpi_component_table <- renderDT({
      wpi_res <- wpi_results_rv()
      if (is.null(wpi_res)) return(empty_dt("Run WPI analysis to see component summary"))
      standard_dt(wpi_res$summary)
    })
    
    output$wpi_segment_table <- renderDT({
      wpi_res <- wpi_results_rv()
      if (is.null(wpi_res)) return(empty_dt("Run WPI analysis to see segment data"))
      df <- wpi_res$data
      base_cols <- c("ID", "tsne1", "tsne2", "WPI", "Target_Score", "Segment")
      pred_cols <- setdiff(names(df), base_cols)
      
      seg_profile <- aggregate(df[, pred_cols, drop = FALSE], by = list(Segment = df$Segment), FUN = mean)
      seg_counts <- as.data.frame(table(df$Segment))
      names(seg_counts) <- c("Segment", "N_Individuals")
      wpi_mean <- aggregate(df[, "WPI", drop = FALSE], by = list(Segment = df$Segment), FUN = mean)
      
      seg_profile <- merge(merge(seg_counts, seg_profile, by = "Segment"), wpi_mean, by = "Segment")
      seg_profile <- seg_profile[, c("Segment", "N_Individuals", "WPI", pred_cols), drop = FALSE]
      
      if (!is.null(wpi_res$header_labels)) names(seg_profile)[match(pred_cols, names(seg_profile))] <- unname(wpi_res$header_labels[pred_cols])
      standard_dt(round_df(seg_profile, 3), selection = "single")
    })
    
    output$wpi_individual_table <- renderDT({ 
      wpi_res <- wpi_results_rv()
      if (is.null(wpi_res)) empty_dt("Run WPI analysis to see individual data") 
      else standard_dt(round_df(wpi_res$data, 2), selection = "single") 
    })
    
    output$download_wpi_data <- downloadHandler(
      filename = function() paste0("wpi_individual_data_", Sys.Date(), ".csv"),
      content = function(file) { req(wpi_results_rv()); write.csv(wpi_results_rv()$data, file, row.names = FALSE) }
    )
    
  })
}