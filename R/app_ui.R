app_ui <- function() {
  page_sidebar(
    fillable = FALSE,
    theme = bs_theme(
      version = 5,
      bootswatch = "flatly",
      primary   = "#1F4E79",
      secondary = "#5C6770",
      success   = "#198754",
      info      = "#0D6EFD",
      warning   = "#FFC107",
      danger    = "#DC3545",
      base_font    = font_google("Inter"),
      heading_font = font_google("Inter")
    ),

    tags$head(
      tags$style(HTML("
        .metric-box { background: #f8fafc; border: 1px solid #e9ecef; border-radius: 12px; padding: 14px 16px; margin-bottom: 10px; }
        .metric-label { font-size: .85rem; color: #5C6770; margin-bottom: 4px; }
        .metric-value { font-size: 1.18rem; font-weight: 700; color: #1F2937; }
        .small-muted { color: #6c757d; font-size: .92rem; }
        .btn-compact { padding: .45rem .75rem; }
        .control-label { font-weight: 600; }
        .action-bar { display: flex; gap: 10px; flex-wrap: wrap; }
        pre { white-space: pre-wrap; }
      "))
    ),

    sidebar = sidebar(
      width = 340,
      div(class = "mb-3", actionButton("btn_run", "▶ Run PLS-SEM", class = "btn-success btn-compact w-100")),
      accordion(
        accordion_panel("Project", card(
          card_header("Project settings"),
          card_body(
            textInput("project_name", "Project name", value = ""),
            fileInput("project_file", "Load project (.rds)", accept = ".rds"),
            downloadButton("btn_save_project", "Save project", class = "btn-success btn-compact w-100")
          )
        )),
        accordion_panel("Data", card(
          card_header("Data import"),
          card_body(
            fileInput("data_file", "Upload Excel file", accept = c(".xls", ".xlsx")),
            textInput("omission_code", "Code to be treated as missing", value = "-1", placeholder = "Example: -1"),
            selectInput("missing_treatment", "Treatment of missing values", choices = c("Listwise deletion (recommended)" = "listwise", "Mean imputation (numeric variables)" = "mean", "Median imputation (numeric variables)" = "median", "None" = "none"), selected = "listwise")
          )
        )),
        accordion_panel("Estimation & Resampling", card(
          card_header("PLS-SEM estimation settings"),
          card_body(
            selectInput("approach_weights", "Weighting approach", choices = c("PLS-PM (PLS)" = "PLS-PM", "GSCA" = "GSCA", "SUMCORR" = "SUMCORR", "ML" = "ML", "OLS" = "OLS", "2SLS" = "2SLS"), selected = "PLS-PM"),
            selectInput("approach_paths", "Path estimation", choices = c("OLS" = "OLS", "2SLS" = "2SLS"), selected = "OLS"),
            selectInput("pls_inner_scheme", "PLS inner weighting", choices = c("Path" = "path", "Centroid" = "centroid", "Factorial"= "factorial"), selected = "path"),
            selectInput("pls_mode_default", "Default outer mode", choices = c("Automatic by construct type" = "auto", "Mode A" = "modeA", "Mode B" = "modeB", "PCA" = "PCA", "Unit" = "unit"), selected = "auto"),
            checkboxInput("plsc_disattenuate", "Use PLSc correction for reflective constructs", value = TRUE)
          )
        ), card(
          card_header("Resampling settings"),
          card_body(
            selectInput("resample_method", "Resample method", choices = c("Bootstrap" = "bootstrap", "Jackknife" = "jackknife", "None" = "none"), selected = "bootstrap"),
            numericInput("n_boot", "Number of resamples", value = 200, min = 0, step = 100),
            selectInput("handle_inadmissibles", "Handle inadmissibles", choices = c("Replace (repeat until admissible)" = "replace", "Drop (discard inadmissible)" = "drop", "Ignore (retain all)" = "ignore"), selected = "replace")
          )
        ))
      ),
      tags$hr(style = "margin-top: auto; margin-bottom: 10px; opacity: 0.3;"),
      div(
        class = "small-muted",
        style = "text-align: center; font-size: 0.82rem; line-height: 1.4;",
        tags$b("ShinyPLSp"), tags$br(),
        "© 2026, Patricio Ramírez Correa.", tags$br(),
        "Open Source Software (GPL-3).", tags$br(),
        tags$span(style = "font-size: 0.75rem;", "Provided 'as is' without warranty.")
      )
    ),

    navset_card_tab(
      id = "main_tabs",
      full_screen = TRUE,
      mod_data_ui("data_tab"),
      mod_constructs_ui("constructs_tab"),
      mod_relations_ui("relations_tab"),
      mod_sem_results_ui("sem_results_tab"),
      mod_wpi_ui("wpi_tab"),
      mod_pathmox_ui("pathmox_tab")
    )
  )
}
