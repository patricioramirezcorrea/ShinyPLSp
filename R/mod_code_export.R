# ==============================================================================
# MÓDULO: EXPORTACIÓN DE SCRIPT R REPRODUCIBLE (cSEM)
# ==============================================================================

mod_code_export_ui <- function(id) {
  ns <- NS(id)
  nav_panel("Exportar Código R",
            layout_columns(
              col_widths = c(4, 8),
              card(
                card_header("Generador de Script PLS-SEM"),
                card_body(
                  p("Esta herramienta compila un script de R reproducible utilizando ", tags$b("cSEM"), " con los parámetros y datos exactos seleccionados en el menú."),
                  tags$hr(),
                  p(tags$b("Archivo vinculado:")),
                  div(style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; font-size: 13px;",
                      textOutput(ns("status_ruta"))
                  ),
                  div(class = "action-bar", style = "margin-top: 30px;",
                      actionButton(ns("btn_generate"), "Generar Código R", class = "btn-primary w-100 mb-2"),
                      downloadButton(ns("download_script"), "Descargar archivo .R", class = "btn-success w-100")
                  )
                )
              ),
              card(
                full_screen = TRUE,
                card_header("Vista Previa del Script Reproducible"),
                card_body(
                  div(style = "position: relative; height: 100%;",
                      verbatimTextOutput(ns("code_preview"), placeholder = TRUE)
                  )
                )
              )
            )
  )
}

mod_code_export_server <- function(id, constructs_rv, relations_rv, omission_code, 
                                   missing_treatment, dataset_name_rv,
                                   approach_weights_rv, 
                                   approach_paths_rv,
                                   pls_weight_scheme_inner_rv,
                                   disattenuate_rv,
                                   resample_method_rv, 
                                   n_boot_rv, 
                                   handle_inadmissibles_rv) { 
  
  moduleServer(id, function(input, output, session) {
    
    generated_code_rv <- reactiveVal("# Haz clic en 'Generar Código R' para compilar el script con la configuración actual.")
    
    output$status_ruta <- renderText({
      ruta <- dataset_name_rv()
      if (is.null(ruta) || ruta == "") "Esperando archivo de datos..." else ruta
    })
    
    observeEvent(input$btn_generate, {
      clist        <- constructs_rv()
      df_rel       <- relations_rv()
      nombre_arch  <- dataset_name_rv()
      omit_val     <- omission_code()
      mt_val       <- missing_treatment()
      a_weights    <- approach_weights_rv()
      a_paths      <- approach_paths_rv()
      inner_wt     <- pls_weight_scheme_inner_rv()
      disatt       <- disattenuate_rv()
      r_method     <- resample_method_rv()
      n_boot_val   <- n_boot_rv()
      h_inadmiss   <- handle_inadmissibles_rv()
      
      if (length(clist) == 0 || nrow(df_rel) == 0) {
        showNotification("No hay modelo definido. Configura constructos y relaciones primero.", type = "warning")
        return()
      }
      
      if (is.null(nombre_arch) || nombre_arch == "") nombre_arch <- "datos.xlsx"
      
      code <- c()
      code <- c(code, "# ==============================================================================")
      code <- c(code, "# SCRIPT R REPRODUCIBLE - ESTIMACIÓN PLS-SEM (cSEM)")
      code <- c(code, paste("# Generado automáticamente el:", Sys.time()))
      code <- c(code, "# ==============================================================================\n")
      code <- c(code, "# --- 1. CARGA DE LIBRERÍAS ---")
      code <- c(code, "library(dplyr)")
      code <- c(code, "library(readxl)")
      code <- c(code, "library(cSEM)\n")
      code <- c(code, "# --- 2. CARGA Y PREPROCESAMIENTO DE DATOS ---")
      code <- c(code, sprintf("raw_data <- as.data.frame(read_excel('%s'))", nombre_arch))
      code <- c(code, "analysis_data <- raw_data\n")
      
      if (!is.null(omit_val) && nzchar(trimws(omit_val))) {
        code <- c(code, sprintf("# Aplicando código de omisión: '%s'", trimws(omit_val)))
        code <- c(code, sprintf("omit_val_str <- '%s'", trimws(omit_val)))
        code <- c(code, "analysis_data[] <- lapply(analysis_data, function(col) {")
        code <- c(code, "  chr_col <- as.character(col)")
        code <- c(code, "  chr_col[chr_col == omit_val_str] <- NA")
        code <- c(code, "  if(!any(is.na(suppressWarnings(as.numeric(na.omit(chr_col)))))) return(as.numeric(chr_col)) else return(chr_col)")
        code <- c(code, "})\n")
      }
      
      if (is.null(mt_val) || mt_val == "") mt_val <- "none"
      code <- c(code, sprintf("# Tratamiento de valores perdidos seleccionado: %s", mt_val))
      if (mt_val == "listwise") {
        code <- c(code, "analysis_data <- analysis_data[complete.cases(analysis_data), , drop = FALSE]\n")
      } else if (mt_val == "mean") {
        code <- c(code, "for(j in seq_along(analysis_data)) {")
        code <- c(code, "  if(is.numeric(analysis_data[[j]]) && any(is.na(analysis_data[[j]]))) {")
        code <- c(code, "    analysis_data[[j]][is.na(analysis_data[[j]])] <- mean(analysis_data[[j]], na.rm = TRUE)")
        code <- c(code, "  }")
        code <- c(code, "}\n")
      } else if (mt_val == "median") {
        code <- c(code, "for(j in seq_along(analysis_data)) {")
        code <- c(code, "  if(is.numeric(analysis_data[[j]]) && any(is.na(analysis_data[[j]]))) {")
        code <- c(code, "    analysis_data[[j]][is.na(analysis_data[[j]])] <- median(analysis_data[[j]], na.rm = TRUE)")
        code <- c(code, "  }")
        code <- c(code, "}\n")
      }
      
      # --- 3. SINTAXIS ---
      code <- c(code, "# --- 3. ESPECIFICACIÓN DEL MODELO PLS-SEM ---")
      code <- c(code, "model_syntax <- \"")
      code <- c(code, "# Modelo Estructural")
      
      targets <- unique(df_rel$target)
      for (tg in targets) {
        origs <- df_rel$origin[df_rel$target == tg]
        code <- c(code, sprintf("  %s ~ %s", tg, paste(origs, collapse = " + ")))
      }
      
      code <- c(code, "\n# Modelo de Medida (El outer mode es automático según el operador =~ o <~)")
      for (nm in names(clist)) {
        items <- clist[[nm]]$items
        typ <- clist[[nm]]$type
        op <- if (identical(typ, "Composite")) "<~" else "=~"
        code <- c(code, sprintf("  %s %s %s", nm, op, paste(items, collapse = " + ")))
      }
      code <- c(code, "\"\n")
      if (is.null(a_weights) || a_weights == "") a_weights <- "PLS-PM"
      if (is.null(a_paths) || a_paths == "") a_paths <- "OLS"
      if (is.null(inner_wt) || inner_wt == "") inner_wt <- "path"
      disatt_str <- if (isTRUE(as.logical(disatt))) "TRUE" else "FALSE"
      
      if (is.null(r_method) || r_method == "") r_method <- "none"
      if (is.null(n_boot_val) || is.na(n_boot_val)) n_boot_val <- 200
      if (is.null(h_inadmiss) || h_inadmiss == "") h_inadmiss <- "replace"
      
      code <- c(code, "# --- 4. ESTIMACIÓN Y BOOTSTRAPPING INTEGRADO ---")
      code <- c(code, "res <- csem(")
      code <- c(code, "  .data = analysis_data,")
      code <- c(code, "  .model = model_syntax,")
      code <- c(code, sprintf("  .approach_weights = '%s',", a_weights))
      code <- c(code, sprintf("  .approach_paths = '%s',", a_paths))
      code <- c(code, sprintf("  .PLS_weight_scheme_inner = '%s',", inner_wt))
      code <- c(code, sprintf("  .disattenuate = %s,", disatt_str))
      
      if (tolower(r_method) != "none") {
        code <- c(code, sprintf("  .resample_method = '%s',", tolower(r_method)))       
        code <- c(code, sprintf("  .R = %d,", as.integer(n_boot_val))) # AQUI ESTABA EL ERROR: Es .R, no .n_boot
        code <- c(code, sprintf("  .handle_inadmissibles = '%s'", tolower(h_inadmiss))) 
      } else {
        code <- c(code, "  .resample_method = 'none'")
      }
      code <- c(code, ")\n")
      code <- c(code, "# --- 5. RESULTADOS Y VISUALIZACIÓN ---")
      
      code <- c(code, "## 1. Gráfico del modelo estimado (Paths)")
      code <- c(code, "plot(res)\n")
      
      code <- c(code, "## 2. Resumen completo de la estimación")
      code <- c(code, "summary_res <- summarize(res)")
      code <- c(code, "print(summary_res)\n")
      
      code <- c(code, "## 3. Criterios de calidad (AVE, R2, HTMT, VIF, etc.)")
      code <- c(code, "assessment_res <- assess(res)")
      code <- c(code, "print(assessment_res)\n")
      
      code <- c(code, "# ==============================================================================")
      code <- c(code, "# FIN DEL SCRIPT")
      code <- c(code, "# ==============================================================================")
      
      final_script_string <- paste(code, collapse = "\n")
      generated_code_rv(final_script_string)
      showNotification("¡Script R generado!", type = "message")
    })
    
    output$code_preview <- renderText({ generated_code_rv() })
    
    output$download_script <- downloadHandler(
      filename = function() paste0("cSEM_analysis_", format(Sys.Date(), "%Y%m%d"), ".R"),
      content = function(file) writeLines(generated_code_rv(), file)
    )
  })
}