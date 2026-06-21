# ==============================================================================
# HELPERS DE RESULTADOS CSEM
# ==============================================================================

get_assess_safe <- function(res) {
  out <- try(cSEM::assess(res), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  out
}

get_summary_safe <- function(res) {
  out <- try(cSEM::summarize(res), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  out
}

get_infer_safe <- function(res) {
  out <- try(cSEM::infer(res), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  out
}

collect_recursive_tables <- function(x, patterns) {
  out <- list()
  if (is.null(x)) return(out)
  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- rep("", length(x))
    for (i in seq_along(x)) {
      nm <- nms[i]
      obj <- x[[i]]
      if ((is.data.frame(obj) || is.matrix(obj)) && any(grepl(paste(patterns, collapse = "|"), nm, ignore.case = TRUE))) {
        out[[length(out) + 1]] <- obj
      }
      if (is.list(obj)) {
        child <- collect_recursive_tables(obj, patterns)
        if (length(child) > 0) out <- c(out, child)
      }
    }
  }
  out
}

pick_best_table <- function(tbls, relation_type = c("generic", "path", "loading")) {
  relation_type <- match.arg(relation_type)
  if (length(tbls) == 0) return(NULL)

  scores <- sapply(tbls, function(obj) {
    df <- as.data.frame(obj, stringsAsFactors = FALSE)
    nms <- names(df)
    score <- ncol(df)
    if (any(grepl("estimate|original|loading|beta", nms, ignore.case = TRUE))) score <- score + 5
    if (any(grepl("std|error|se", nms, ignore.case = TRUE))) score <- score + 5
    if (any(grepl("t.*value|t.*stat|^t$", nms, ignore.case = TRUE))) score <- score + 5
    if (any(grepl("p.*value|^p$", nms, ignore.case = TRUE))) score <- score + 5
    if (relation_type == "path" && any(grepl("dependent|target|predictor|origin|path|relation", nms, ignore.case = TRUE))) score <- score + 5
    if (relation_type == "loading" && any(grepl("construct|indicator|item|manifest|loading|relation", nms, ignore.case = TRUE))) score <- score + 5
    score
  })
  tbls[[which.max(scores)]]
}

normalize_table_source <- function(x, relation_type = c("generic", "path", "loading")) {
  relation_type <- match.arg(relation_type)
  if (is.null(x)) return(NULL)

  if (is.data.frame(x)) {
    df <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
    # Normalizar nombres: espacios y puntos → underscore, limpiar duplicados
    names(df) <- gsub("[[:space:]]+", "_", names(df))
    names(df) <- gsub("\\.+",         "_", names(df))
    names(df) <- gsub("_+",           "_", names(df))
    names(df) <- gsub("_$",           "",  names(df))
    return(df)
  }

  if (is.matrix(x)) {
    if (is.numeric(x) && !is.null(rownames(x)) && !is.null(colnames(x))) {
      if (relation_type == "path")    return(matrix_to_long(x, "Dependent",  "Predictor", "Estimate", drop_zeros = TRUE))
      if (relation_type == "loading") return(matrix_to_long(x, "Construct",  "Indicator", "Estimate", drop_zeros = TRUE))
    }
    df <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
    names(df) <- gsub("[[:space:]]+", "_", names(df))
    names(df) <- gsub("\\.+",         "_", names(df))
    names(df) <- gsub("_+",           "_", names(df))
    names(df) <- gsub("_$",           "",  names(df))
    return(df)
  }

  NULL
}

parse_path_relation <- function(rel) {
  rel <- trimws(rel)
  if (grepl("~", rel, fixed = TRUE)) {
    parts <- strsplit(rel, "~", fixed = TRUE)[[1]]
    if (length(parts) >= 2) return(c(trimws(parts[1]), trimws(parts[2])))
  }
  if (grepl("->", rel, fixed = TRUE)) {
    parts <- strsplit(rel, "->", fixed = TRUE)[[1]]
    if (length(parts) >= 2) return(c(trimws(parts[2]), trimws(parts[1])))
  }
  c(NA_character_, NA_character_)
}

parse_loading_relation <- function(rel) {
  rel <- trimws(rel)
  if (grepl("=~", rel, fixed = TRUE)) {
    parts <- strsplit(rel, "=~", fixed = TRUE)[[1]]
    if (length(parts) >= 2) return(c(trimws(parts[1]), trimws(parts[2])))
  }
  if (grepl("<~", rel, fixed = TRUE)) {
    parts <- strsplit(rel, "<~", fixed = TRUE)[[1]]
    if (length(parts) >= 2) return(c(trimws(parts[1]), trimws(parts[2])))
  }
  if (grepl("->", rel, fixed = TRUE)) {
    parts <- strsplit(rel, "->", fixed = TRUE)[[1]]
    if (length(parts) >= 2) return(c(trimws(parts[1]), trimws(parts[2])))
  }
  c(NA_character_, NA_character_)
}

parse_paths_table <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  d <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  nms <- names(d)

  dep_col <- find_name(nms, c("^Dependent$", "^Target$", "^Endogenous$", "dependent", "target", "endo"))
  pred_col <- find_name(nms, c("^Predictor$", "^Origin$", "^Independent$", "^Exogenous$", "predictor", "origin", "independent", "exogenous"))
  rel_col <- find_name(nms, c("^Path$", "^Relation$", "^Parameter$", "path", "relation", "parameter"))

  dep <- pred <- NULL
  if (!is.na(dep_col) && !is.na(pred_col)) {
    dep <- as.character(d[[dep_col]])
    pred <- as.character(d[[pred_col]])
  } else if (!is.na(rel_col)) {
    pairs <- t(vapply(as.character(d[[rel_col]]), parse_path_relation, character(2)))
    dep <- pairs[, 1]
    pred <- pairs[, 2]
  } else if (ncol(d) >= 2) {
    dep <- as.character(d[[1]])
    pred <- as.character(d[[2]])
  } else {
    return(NULL)
  }

  est_col <- find_name(nms, c("^Estimate$", "^Est$", "estimate", "original", "beta"))
  se_col  <- find_name(nms, c("^SE$", "std.*error", "standard.*error", "^Std\\. Error$", "std.*err"))
  t_col   <- find_name(nms, c("^t$", "t.*value", "t.*stat"))
  p_col   <- find_name(nms, c("^p$", "p.*value"))
  cil_col <- find_name(nms, c("ci.*low", "ci.*lower", "lower.*ci", "2\\.5", "lower"))
  ciu_col <- find_name(nms, c("ci.*up", "ci.*upper", "upper.*ci", "97\\.5", "upper"))

  if (is.na(est_col)) {
    num_candidates <- names(d)[sapply(d, is.numeric)]
    if (length(num_candidates) > 0) est_col <- num_candidates[1]
  }

  out <- data.frame(
    Dependent = dep,
    Predictor = pred,
    Estimate = if (!is.na(est_col)) to_num(d[[est_col]]) else NA_real_,
    Std_Error = if (!is.na(se_col)) to_num(d[[se_col]]) else NA_real_,
    T_value = if (!is.na(t_col)) to_num(d[[t_col]]) else NA_real_,
    P_value = if (!is.na(p_col)) to_num(d[[p_col]]) else NA_real_,
    CI_low = if (!is.na(cil_col)) to_num(d[[cil_col]]) else NA_real_,
    CI_high = if (!is.na(ciu_col)) to_num(d[[ciu_col]]) else NA_real_,
    stringsAsFactors = FALSE
  )

  keep <- !(is.na(out$Dependent) | is.na(out$Predictor) | out$Dependent == "" | out$Predictor == "")
  out <- out[keep, , drop = FALSE]
  out$key <- paste0(out$Dependent, "|||", out$Predictor)
  out$Significance <- significance_stars(out$P_value)
  rownames(out) <- NULL
  out
}

parse_loadings_table <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  d <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  nms <- names(d)

  con_col <- find_name(nms, c("^Construct$", "^Latent$", "construct", "latent", "factor"))
  ind_col <- find_name(nms, c("^Indicator$", "^Item$", "^Manifest$", "indicator", "item", "manifest"))
  rel_col <- find_name(nms, c("^Loading$", "^Relation$", "^Parameter$", "loading", "relation", "parameter"))

  cons <- inds <- NULL
  if (!is.na(con_col) && !is.na(ind_col)) {
    cons <- as.character(d[[con_col]])
    inds <- as.character(d[[ind_col]])
  } else if (!is.na(rel_col) && !any(grepl("^Estimate$|estimate|original", nms, ignore.case = TRUE))) {
    pairs <- t(vapply(as.character(d[[rel_col]]), parse_loading_relation, character(2)))
    cons <- pairs[, 1]
    inds <- pairs[, 2]
  } else if (ncol(d) >= 2) {
    cons <- as.character(d[[1]])
    inds <- as.character(d[[2]])
  } else {
    return(NULL)
  }

  est_col <- find_name(nms, c("^Estimate$", "^Loading$", "estimate", "original", "loading"))
  se_col  <- find_name(nms, c("^SE$", "std.*error", "standard.*error", "^Std\\. Error$", "std.*err"))
  t_col   <- find_name(nms, c("^t$", "t.*value", "t.*stat"))
  p_col   <- find_name(nms, c("^p$", "p.*value"))
  cil_col <- find_name(nms, c("ci.*low", "ci.*lower", "lower.*ci", "2\\.5", "lower"))
  ciu_col <- find_name(nms, c("ci.*up", "ci.*upper", "upper.*ci", "97\\.5", "upper"))

  if (is.na(est_col)) {
    num_candidates <- names(d)[sapply(d, is.numeric)]
    if (length(num_candidates) > 0) est_col <- num_candidates[1]
  }

  out <- data.frame(
    Construct = cons,
    Indicator = inds,
    Estimate = if (!is.na(est_col)) to_num(d[[est_col]]) else NA_real_,
    Std_Error = if (!is.na(se_col)) to_num(d[[se_col]]) else NA_real_,
    T_value = if (!is.na(t_col)) to_num(d[[t_col]]) else NA_real_,
    P_value = if (!is.na(p_col)) to_num(d[[p_col]]) else NA_real_,
    CI_low = if (!is.na(cil_col)) to_num(d[[cil_col]]) else NA_real_,
    CI_high = if (!is.na(ciu_col)) to_num(d[[ciu_col]]) else NA_real_,
    stringsAsFactors = FALSE
  )

  keep <- !(is.na(out$Construct) | is.na(out$Indicator) | out$Construct == "" | out$Indicator == "")
  out <- out[keep, , drop = FALSE]
  out$key <- paste0(out$Construct, "|||", out$Indicator)
  out$Significance <- significance_stars(out$P_value)
  rownames(out) <- NULL
  out
}

get_construct_types <- function(res) {
  if (is.null(res)) return(NULL)
  ct <- NULL
  if (!is.null(res$Information) && !is.null(res$Information$Model) && !is.null(res$Information$Model$construct_type)) {
    ct <- res$Information$Model$construct_type
  }
  if (is.null(ct)) return(NULL)
  data.frame(Construct = names(ct), Type = as.vector(ct), stringsAsFactors = FALSE)
}

extract_r2 <- function(res) {
  if (is.null(res)) return(NULL)
  b <- get_assess_safe(res)
  if (!is.null(b) && !is.null(b$R2)) {
    r2 <- named_to_df(b$R2, "Construct", "R_squared")
    r2a <- if (!is.null(b$R2_adj)) named_to_df(b$R2_adj, "Construct", "Adjusted_R_squared") else NULL
    return(round_df(merge_by_first_col(list(r2, r2a)), 3))
  }
  if (!is.null(res$Estimates) && !is.null(res$Estimates$R2)) {
    r2 <- named_to_df(res$Estimates$R2, "Construct", "R_squared")
    r2a <- if (!is.null(res$Estimates$R2adj)) named_to_df(res$Estimates$R2adj, "Construct", "Adjusted_R_squared") else NULL
    return(round_df(merge_by_first_col(list(r2, r2a)), 3))
  }
  NULL
}

extract_paths <- function(res) {
  if (is.null(res)) return(NULL)
  s <- get_summary_safe(res)
  inf <- get_infer_safe(res)

  tbls <- c(
    collect_recursive_tables(s,  c("^Path_estimates$", "^Paths?$", "path.*estimate", "structural.*path")),
    collect_recursive_tables(inf, c("^Path_estimates$", "^Paths?$", "path.*estimate", "structural.*path"))
  )

  best <- pick_best_table(tbls, "path")
  best <- normalize_table_source(best, "path")

  # Normalizar notación de rutas: "IU ~ PE" -> "PE -> IU"
  if (!is.null(best) && nrow(best) > 0) {
    if ("Name" %in% names(best)) {
      parts <- strsplit(as.character(best$Name), "~", fixed = TRUE)
      lhs <- vapply(parts, function(z) trimws(z[1]), character(1))
      rhs <- vapply(parts, function(z) if (length(z) >= 2) trimws(z[2]) else NA_character_, character(1))
      best$Name <- paste(rhs, "->", lhs)
    }
    return(round_df(best, 3))
  }

  if (!is.null(res$Estimates) && !is.null(res$Estimates$Path_estimates)) {
    df <- matrix_to_long(
      res$Estimates$Path_estimates,
      row_name  = "Dependent",
      col_name  = "Predictor",
      value_name= "Estimate",
      drop_zeros= TRUE
    )
    df <- parse_paths_table(df)
    df$Name <- paste(df$Predictor, "->", df$Dependent)
    return(round_df(df, 3))
  }
  NULL
}

extract_loadings <- function(res) {
  if (is.null(res)) return(NULL)
  s <- get_summary_safe(res)
  inf <- get_infer_safe(res)

  tbls <- c(
    collect_recursive_tables(s, c("^Loading_estimates$", "^Loadings?$", "loading.*estimate", "outer.*loading")),
    collect_recursive_tables(inf, c("^Loading_estimates$", "^Loadings?$", "loading.*estimate", "outer.*loading"))
  )

  best <- pick_best_table(tbls, "loading")
  best <- normalize_table_source(best, "loading")
  if (!is.null(best) && nrow(best) > 0) return(round_df(best, 3))

  if (!is.null(res$Estimates) && !is.null(res$Estimates$Loading_estimates)) {
    df <- matrix_to_long(res$Estimates$Loading_estimates, row_name = "Construct", col_name = "Indicator", value_name = "Estimate", drop_zeros = TRUE)
    df <- parse_loadings_table(df)
    return(round_df(df, 3))
  }
  NULL
}

extract_model_fit <- function(res) {
  b <- get_assess_safe(res)
  if (is.null(b)) return(NULL)

  fit_names <- c("SRMR", "d_ULS", "d_G", "dL", "RMS_theta", "GoF", "NFI", "GFI")
  out <- data.frame(Index = character(), Value = numeric(), stringsAsFactors = FALSE)

  recurse_fit <- function(x, nm = "") {
    if (is.null(x)) return(NULL)
    if (is.numeric(x) && length(x) == 1 && nzchar(nm) && nm %in% fit_names) {
      out <<- rbind(out, data.frame(Index = nm, Value = as.numeric(x), stringsAsFactors = FALSE))
      return(NULL)
    }
    if (is.numeric(x) && length(x) > 1 && !is.null(names(x))) {
      for (i in seq_along(x)) recurse_fit(x[i], names(x)[i])
      return(NULL)
    }
    if (is.list(x)) {
      nms <- names(x)
      if (is.null(nms)) nms <- rep("", length(x))
      for (i in seq_along(x)) recurse_fit(x[[i]], nms[i])
    }
  }

  recurse_fit(b)
  if (nrow(out) == 0) return(NULL)
  out <- out[!duplicated(out$Index), , drop = FALSE]
  round_df(out, 3)
}

extract_mm_quality <- function(res) {
  b <- get_assess_safe(res)
  if (is.null(b)) return(list(quality = NULL, htmt = NULL, fornell = NULL))
  omega_df <- NULL
  rel <- b$Reliability
  if (!is.null(rel)) {
      if (!is.null(rel$omega)) {
      omega_df <- named_to_df(rel$omega, "Construct", "Omega_McDonald")
    }
    if (is.null(omega_df) && !is.null(rel$rho_T)) {
      omega_df <- named_to_df(rel$rho_T, "Construct", "Omega_McDonald")
    }
    if (is.null(omega_df) && !is.null(res$Estimates$Loading_estimates)) {
      L <- res$Estimates$Loading_estimates
      if (is.matrix(L)) {
        omega_vals <- apply(L, 1, function(row) {
          lam <- row[row != 0]
          if (length(lam) < 2) return(NA_real_)
          num   <- sum(lam)^2
          denom <- num + sum(1 - lam^2)
          if (denom == 0) NA_real_ else num / denom
        })
        omega_df <- data.frame(
          Construct      = names(omega_vals),
          Omega_McDonald = as.numeric(omega_vals),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  quality <- merge_by_first_col(list(
    if (!is.null(b$AVE))                                  named_to_df(b$AVE,                                  "Construct", "AVE")           else NULL,
    if (!is.null(rel$Cronbachs_alpha))                    named_to_df(rel$Cronbachs_alpha,                    "Construct", "Cronbach_alpha") else NULL,
    omega_df,
    if (!is.null(rel$Joereskogs_rho))                     named_to_df(rel$Joereskogs_rho,                     "Construct", "rhoC_CR")        else NULL,
    if (!is.null(rel$`Dijkstra-Henselers_rho_A`))         named_to_df(rel$`Dijkstra-Henselers_rho_A`,         "Construct", "rhoA")           else NULL
  ))


  htmt <- NULL
  if (!is.null(b$HTMT) && !is.null(b$HTMT$htmts)) {
    htmt <- as.data.frame(b$HTMT$htmts, stringsAsFactors = FALSE, check.names = FALSE)
    nr <- nrow(htmt)
    nc <- ncol(htmt)
    for (i in seq_len(nr)) {
      for (j in seq_len(nc)) {
        if (j >= i) htmt[i, j] <- NA
      }
    }
    htmt <- cbind(Construct = rownames(b$HTMT$htmts), htmt, row.names = NULL)
  }

  fornell <- NULL
  if (!is.null(b$`Fornell-Larcker`)) {
    fornell <- as.data.frame(
      b$`Fornell-Larcker`,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    fornell <- cbind(
      Construct = rownames(fornell),
      fornell,
      row.names = NULL
    )
  }

  list(quality = round_df(quality, 3), htmt = round_df(htmt, 3), fornell = round_df(fornell, 3))
}

