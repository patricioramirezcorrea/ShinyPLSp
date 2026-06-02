# ==============================================================================
# HELPERS GENERALES
# ==============================================================================

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || identical(x, "")) y else x
}

sanitize_name <- function(x) {
  x <- trimws(x)
  gsub("\\s+", "_", x)
}

to_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

round_df <- function(df, k = 3) {
  if (is.null(df)) return(NULL)
  df[] <- lapply(df, function(col) if (is.numeric(col)) round(col, k) else col)
  df
}

format_df <- function(df, k = 3) {
  if (is.null(df)) return(NULL)
  df[] <- lapply(df, function(col) if (is.numeric(col)) formatC(col, digits = k, format = "f") else col)
  df
}

named_to_df <- function(x, name_col = "name", value_col = "value") {
  if (is.null(x) || length(x) == 0) return(NULL)
  df <- data.frame(
    temp_name = names(x),
    temp_value = as.numeric(x),
    stringsAsFactors = FALSE
  )
  names(df) <- c(name_col, value_col)
  df
}

matrix_to_long <- function(mat, row_name = "row", col_name = "col", value_name = "value", drop_zeros = FALSE, drop_diag = FALSE) {
  if (is.null(mat)) return(NULL)
  out <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(out) <- c(row_name, col_name, value_name)
  if (drop_diag) out <- out[out[[row_name]] != out[[col_name]], , drop = FALSE]
  if (drop_zeros) out <- out[!is.na(out[[value_name]]) & out[[value_name]] != 0, , drop = FALSE]
  rownames(out) <- NULL
  out
}

merge_by_first_col <- function(df_list) {
  df_list <- Filter(function(x) !is.null(x) && nrow(x) > 0, df_list)
  if (length(df_list) == 0) return(NULL)
  if (length(df_list) == 1) return(df_list[[1]])
  key <- names(df_list[[1]])[1]
  Reduce(function(x, y) merge(x, y, by = key, all = TRUE), df_list)
}

find_name <- function(nms, patterns) {
  if (is.null(nms) || length(nms) == 0) return(NA_character_)
  for (p in patterns) {
    hit <- nms[grepl(p, nms, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

significance_stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
}

empty_dt <- function(message) {
  datatable(data.frame(Message = message), options = list(dom = "t"), rownames = FALSE)
}

standard_dt <- function(df, digits = NULL, selection = "single") {
  if (!is.null(digits)) df <- format_df(df, digits)
  datatable(df, selection = selection, options = list(dom = "tip", pageLength = 100, scrollX = TRUE), rownames = FALSE)
}

metric_box <- function(label, output_id) {
  div(
    class = "metric-box",
    div(class = "metric-label", label),
    div(class = "metric-value", textOutput(output_id, inline = TRUE))
  )
}