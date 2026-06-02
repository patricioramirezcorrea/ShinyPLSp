# ==============================================================================
# HELPERS DE DATOS Y PREPROCESAMIENTO
# ==============================================================================

prepare_pathmox_segvars <- function(df, config_df) {
  if (is.null(df) || is.null(config_df) || nrow(config_df) == 0) return(NULL)
  out <- vector("list", nrow(config_df))
  names(out) <- config_df$Variable
  
  for (i in seq_len(nrow(config_df))) {
    v <- config_df$Variable[i]
    p <- config_df$Processing[i]
    b <- config_df$Bins[i]
    levels_txt <- if ("Levels" %in% names(config_df)) {
      lv <- config_df$Levels[i]
      if (is.na(lv)) "" else as.character(lv)
    } else ""
    
    ordered_flag <- if ("Ordered" %in% names(config_df)) isTRUE(config_df$Ordered[i]) else FALSE
    if (!v %in% names(df)) next
    x <- df[[v]]
    
    if (p == "factor" || p == "character") {
      x_chr <- trimws(as.character(x))
      x_chr[!nzchar(x_chr)] <- NA
      levs <- if (nzchar(levels_txt)) {
        trimws(unlist(strsplit(levels_txt, ",", fixed = TRUE)))
      } else {
        sort(unique(x_chr[!is.na(x_chr)]))
      }
      levs <- levs[nzchar(levs)]
      out[[v]] <- if (ordered_flag) ordered(x_chr, levels = levs) else factor(x_chr, levels = levs)
      
    } else if (p %in% c("quantile", "equal")) {
      xn <- if (is.numeric(x)) x else suppressWarnings(as.numeric(as.character(x)))
      if (all(is.na(xn))) stop(sprintf("Variable '%s' no valid numeric values.", v))
      rng <- range(xn, na.rm = TRUE)
      if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) stop(sprintf("Variable '%s' insufficient range.", v))
      
      if (p == "quantile") {
        qs <- unique(stats::quantile(xn, probs = seq(0, 1, length.out = b + 1), na.rm = TRUE, type = 7))
        if (length(qs) < 2) stop(sprintf("Variable '%s' insufficient variation.", v))
        out[[v]] <- cut(xn, breaks = qs, include.lowest = TRUE, ordered_result = TRUE)
      } else {
        br <- unique(seq(rng[1], rng[2], length.out = b + 1))
        if (length(br) < 2) stop(sprintf("Variable '%s' insufficient variation.", v))
        out[[v]] <- cut(xn, breaks = br, include.lowest = TRUE, ordered_result = TRUE)
      }
    } else {
      x_chr <- trimws(as.character(x))
      x_chr[!nzchar(x_chr)] <- NA
      out[[v]] <- if (ordered_flag) ordered(x_chr) else factor(x_chr)
    }
  }
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) return(NULL)
  
  out_df <- data.frame(row.names = seq_len(length(out[[1]])))
  for (col_name in names(out)) out_df[[col_name]] <- out[[col_name]]
  out_df
}

extract_first_df <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) return(x)
  if (is.list(x)) {
    for (nm in names(x)) {
      obj <- extract_first_df(x[[nm]])
      if (!is.null(obj)) return(obj)
    }
  }
  NULL
}

constructs_to_df <- function(clist) {
  if (length(clist) == 0) {
    return(data.frame(Construct = character(), Type = character(), Operator = character(), Indicators = character(), stringsAsFactors = FALSE))
  }
  data.frame(
    Construct = names(clist),
    Type = vapply(clist, function(x) x$type, character(1)),
    Operator = vapply(clist, function(x) ifelse(x$type == "Composite", "<~", "=~"), character(1)),
    Indicators = vapply(clist, function(x) paste(x$items, collapse = ", "), character(1)),
    stringsAsFactors = FALSE
  )
}

apply_omission_code <- function(df, omission_code) {
  if (is.null(df)) return(NULL)
  omission_code <- trimws(omission_code)
  if (!nzchar(omission_code)) return(df)
  num_code <- suppressWarnings(as.numeric(omission_code))
  
  for (j in seq_along(df)) {
    col <- df[[j]]
    if (inherits(col, "Date") || inherits(col, "POSIXct") || inherits(col, "POSIXt")) next
    if (is.numeric(col) || is.integer(col)) {
      if (!is.na(num_code)) col[col == num_code] <- NA
      df[[j]] <- col
    } else {
      chr <- trimws(as.character(col))
      idx_chr <- chr == omission_code
      idx_num <- rep(FALSE, length(chr))
      if (!is.na(num_code)) {
        suppressWarnings(num_view <- as.numeric(chr))
        idx_num <- !is.na(num_view) & num_view == num_code
      }
      chr[idx_chr | idx_num] <- NA
      df[[j]] <- chr
    }
  }
  df
}

apply_missing_treatment <- function(df, method = "listwise") {
  if (is.null(df)) return(list(data = NULL, removed = 0, imputed = 0))
  method <- match.arg(method, c("listwise", "mean", "median", "none"))
  original_n <- nrow(df)
  
  if (method == "none") return(list(data = df, removed = 0, imputed = sum(is.na(df))))
  if (method == "listwise") {
    keep <- stats::complete.cases(df)
    out <- df[keep, , drop = FALSE]
    return(list(data = out, removed = original_n - nrow(out), imputed = 0))
  }
  
  out <- df
  imputed_count <- 0
  for (j in seq_along(out)) {
    col <- out[[j]]
    if (is.numeric(col) || is.integer(col)) {
      miss <- is.na(col)
      if (any(miss)) {
        fill_value <- if (method == "mean") mean(col, na.rm = TRUE) else stats::median(col, na.rm = TRUE)
        if (is.finite(fill_value)) {
          col[miss] <- fill_value
          imputed_count <- imputed_count + sum(miss)
        }
      }
      out[[j]] <- col
    }
  }
  out <- out[stats::complete.cases(out), , drop = FALSE]
  list(data = out, removed = original_n - nrow(out), imputed = imputed_count)
}