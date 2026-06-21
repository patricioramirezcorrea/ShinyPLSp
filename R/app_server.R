app_server <- function(input, output, session) {

  session$onSessionEnded(function() {
    stopApp()
  })

  # ----------------------------------------------------------------------------
  # DECLARACIÓN DE ESTADOS GLOBALES (STATE LIFTING)
  # ----------------------------------------------------------------------------
  loading_project_rv <- reactiveVal(FALSE)
  raw_data_rv <- reactiveVal(NULL)
  analysis_data_aug_rv <- reactiveVal(NULL)

  csem_cache <- new.env(parent = emptyenv())
  make_cache_key <- function(...) {
    digest::digest(list(...), algo = "xxhash64")
  }

  pmx_detail_cache <- reactiveVal(list())
  constructs_rv <- reactiveVal(list())
  relations_rv <- reactiveVal(data.frame(origin = character(), target = character(), stringsAsFactors = FALSE))
  result_rv <- reactiveVal(NULL)

  items_available_rv <- reactiveVal(character(0))
  items_selected_rv <- reactiveVal(character(0))

  construct_selected_idx <- reactiveVal(NULL)
  relation_selected_idx <- reactiveVal(NULL)

  pathmox_levels_available_rv <- reactiveVal(character(0))
  pathmox_levels_selected_rv  <- reactiveVal(character(0))

  pathmox_segvars_rv <- reactiveVal(data.frame(Variable = character(), Processing = character(), Bins = numeric(), Levels = character(), Ordered = logical(), stringsAsFactors = FALSE))
  pathmox_results_rv <- reactiveVal(NULL)
  pathmox_micom_rv   <- reactiveVal(NULL)
  pathmox_mga_rv     <- reactiveVal(NULL)

  # ----------------------------------------------------------------------------
  # CARGA DE DATOS Y BUNDLE REACTIVO
  # ----------------------------------------------------------------------------
  observeEvent(input$data_file, {
    req(input$data_file)
    df <- as.data.frame(read_excel(input$data_file$datapath))
    raw_data_rv(df)

    # Resets al cargar nueva base de datos
    constructs_rv(list())
    relations_rv(data.frame(origin = character(), target = character(), stringsAsFactors = FALSE))
    result_rv(NULL)
    construct_selected_idx(NULL)
    relation_selected_idx(NULL)

    all_cols <- colnames(df)
    items_available_rv(all_cols)
    items_selected_rv(character(0))

    # Notificación
    showNotification("Data loaded successfully.", type = "message")
  })

  analysis_bundle <- reactive({
    df <- raw_data_rv()
    if (is.null(df)) return(NULL)
    df1 <- apply_omission_code(as.data.frame(df), input$omission_code)
    apply_missing_treatment(df1, input$missing_treatment)
  })

  analysis_data <- reactive({
    bun <- analysis_bundle()
    if (is.null(bun)) return(NULL)
    bun$data
  })

  observe({
    analysis_data_aug_rv(analysis_data())
  })

  observeEvent(c(input$omission_code, input$missing_treatment), {
    if (!is.null(raw_data_rv())) result_rv(NULL)
  }, ignoreInit = TRUE)

  # ----------------------------------------------------------------------------
  # INYECCIÓN DE MÓDULOS (LLAMADAS A LOS SERVIDORES DE MÓDULOS)
  # ----------------------------------------------------------------------------

  # 1. Módulo Data
  mod_data_server("data_tab", analysis_data_aug_rv, analysis_bundle,
                  reactive(input$omission_code), reactive(input$missing_treatment))

  # 2. Módulo Constructs
  mod_constructs_server("constructs_tab", analysis_data, constructs_rv, relations_rv, result_rv,
                        items_available_rv, items_selected_rv, construct_selected_idx, relation_selected_idx)

  # 3. Módulo Relations (Capturamos la sintaxis del modelo que retorna)
  relations_out <- mod_relations_server("relations_tab", constructs_rv, relations_rv, result_rv)

  # 4. Módulo SEM Results
  mod_sem_results_server("sem_results_tab", result_rv, constructs_rv, relations_rv, analysis_data,
                         analysis_bundle, raw_data_rv, reactive(input$resample_method),
                         reactive(input$n_boot), reactive(input$project_name))
  # 5. Módulo PATHMOX
  mod_pathmox_server("pathmox_tab", analysis_data_aug_rv, relations_out$model_lavaan, result_rv,
                     pathmox_segvars_rv, pathmox_results_rv, pathmox_micom_rv, pathmox_mga_rv, pmx_detail_cache,
                     pathmox_levels_available_rv, pathmox_levels_selected_rv,
                     reactive(input$handle_inadmissibles), reactive(input$approach_weights), reactive(input$approach_paths),
                     reactive(input$pls_inner_scheme), reactive(input$plsc_disattenuate), reactive(input$resample_method), reactive(input$n_boot))

  # ----------------------------------------------------------------------------
  # CSEM ESTIMATION (CORE LOGIC Y CACHÉ)
  # ----------------------------------------------------------------------------
  observeEvent(input$btn_run, {
    req(analysis_data())

    # Extraemos el modelo lavaan desde el retorno del módulo relations
    model <- relations_out$model_lavaan()

    if (identical(model, "No constructs defined yet.")) {
      showNotification("Define constructs first.", type = "error");
      return()
    }

    dta <- analysis_data()
    if (anyNA(dta)) {
      showNotification("Processed data still contains NAs.", type = "error");
      return()
    }

    n_boot   <- if (input$resample_method == "none") 0 else input$n_boot
    pls_modes <- if (input$pls_mode_default != "auto") input$pls_mode_default else NULL

    # Revisión en Caché
    cache_key <- make_cache_key(
      data = dta,
      model = model,
      approach_weights = input$approach_weights,
      approach_paths = input$approach_paths,
      pls_inner_scheme = input$pls_inner_scheme,
      pls_modes = pls_modes,
      resample_method = input$resample_method,
      n_boot = n_boot,
      handle_inadmissibles = input$handle_inadmissibles,
      disattenuate = input$plsc_disattenuate
    )

    if (exists(cache_key, envir = csem_cache, inherits = FALSE)) {
      result_rv(get(cache_key, envir = csem_cache))
      showNotification("Model loaded from cache.", type = "message")
      updateTabsetPanel(session, "main_tabs", selected = "SEM results")
      return()
    }

    withProgress(message = "Estimating model...", value = 0.5, {
      res <- try(
        cSEM::csem(
          .data = dta,
          .model = model,
          .approach_weights = input$approach_weights,
          .approach_paths = input$approach_paths,
          .PLS_weight_scheme_inner = input$pls_inner_scheme,
          .PLS_modes = pls_modes,
          .resample_method = input$resample_method,
          .R = n_boot,
          .handle_inadmissibles = input$handle_inadmissibles,
          .eval_plan            = "multisession",
          .disattenuate = input$plsc_disattenuate
        ),
        silent = TRUE
      )

      if (inherits(res, "try-error")) {
        result_rv(NULL)
        showNotification(paste("Estimation failed:", gsub("Error in .*?:", "", as.character(res))), type = "error")
      } else {
        # GUARDAR EN CACHÉ
        assign(cache_key, res, envir = csem_cache)

        result_rv(res)
        showNotification("Model estimated successfully.", type = "message")
        updateTabsetPanel(session, "main_tabs", selected = "SEM results")
      }
    })
  })

  # ----------------------------------------------------------------------------
  # GUARDAR Y CARGAR PROYECTO
  # ----------------------------------------------------------------------------
  output$btn_save_project <- downloadHandler(
    filename = function() {
      paste0(sanitize_name(trimws(input$project_name)), "_PLS-SEM_", Sys.Date(), ".rds")
    },
    content = function(file) {
      saveRDS(
        list(
          name = input$project_name,
          raw_data = raw_data_rv(),
          constructs = constructs_rv(),
          relations = relations_rv(),
          pathmox_segvars = pathmox_segvars_rv(),
          settings = list(
            omission_code = input$omission_code,
            missing_treatment = input$missing_treatment,
            approach_weights = input$approach_weights,
            approach_paths = input$approach_paths,
            pls_inner_scheme = input$pls_inner_scheme,
            pls_mode_default = input$pls_mode_default,
            plsc_disattenuate = input$plsc_disattenuate,
            resample_method = input$resample_method,
            n_boot = input$n_boot,
            handle_inadmissibles = input$handle_inadmissibles
          )
        ),
        file = file
      )
    }
  )

  observeEvent(input$project_file, {
    req(input$project_file)

    loading_project_rv(TRUE)

    project <- try(readRDS(input$project_file$datapath), silent = TRUE)
    if (inherits(project, "try-error")) {
      showNotification("Invalid project file.", type = "error")
      loading_project_rv(FALSE)
      return(NULL)
    }

    # Restaurar Inputs Visuales
    updateTextInput(session, "project_name", value = project$name %||% "")
    if (!is.null(project$settings$omission_code)) updateTextInput(session, "omission_code", value = project$settings$omission_code)
    if (!is.null(project$settings$missing_treatment)) updateSelectInput(session, "missing_treatment", selected = project$settings$missing_treatment)
    if (!is.null(project$settings$approach_weights)) updateSelectInput(session, "approach_weights", selected = project$settings$approach_weights)
    if (!is.null(project$settings$approach_paths)) updateSelectInput(session, "approach_paths", selected = project$settings$approach_paths)
    if (!is.null(project$settings$pls_inner_scheme)) updateSelectInput(session, "pls_inner_scheme", selected = project$settings$pls_inner_scheme)
    if (!is.null(project$settings$pls_mode_default)) updateSelectInput(session, "pls_mode_default", selected = project$settings$pls_mode_default)
    if (!is.null(project$settings$plsc_disattenuate)) updateCheckboxInput(session, "plsc_disattenuate", value = project$settings$plsc_disattenuate)
    if (!is.null(project$settings$resample_method)) updateSelectInput(session, "resample_method", selected = project$settings$resample_method)
    if (!is.null(project$settings$n_boot)) updateNumericInput(session, "n_boot", value = project$settings$n_boot)
    if (!is.null(project$settings$handle_inadmissibles)) updateSelectInput(session, "handle_inadmissibles", selected = project$settings$handle_inadmissibles)

    # Restaurar variables globales
    raw_data_rv(if (!is.null(project$raw_data)) as.data.frame(project$raw_data) else NULL)

    all_cols <- if (!is.null(project$raw_data)) colnames(project$raw_data) else character(0)
    used_items <- unique(unlist(lapply(project$constructs, `[[`, "items"), use.names = FALSE))
    items_available_rv(sort(setdiff(all_cols, used_items)))
    items_selected_rv(character(0))

    constructs_rv(if (!is.null(project$constructs)) project$constructs else list())
    relations_rv(if (!is.null(project$relations)) project$relations else data.frame(origin = character(), target = character(), stringsAsFactors = FALSE))
    pathmox_segvars_rv( if (!is.null(project$pathmox_segvars)) project$pathmox_segvars  else data.frame(Variable = character(), Processing = character(), Bins = numeric(), Levels = character(), Ordered = logical(), stringsAsFactors = FALSE))

    result_rv(NULL)
    pathmox_results_rv(NULL)
    pathmox_micom_rv(NULL)
    pathmox_mga_rv(NULL)
    pmx_detail_cache(list())

    construct_selected_idx(NULL)
    relation_selected_idx(NULL)
    pathmox_levels_available_rv(character(0))
    pathmox_levels_selected_rv(character(0))

    session$onFlushed(function() {
       loading_project_rv(FALSE)
    }, once = TRUE)

    showNotification("Project loaded successfully.", type = "message")
  })

}
