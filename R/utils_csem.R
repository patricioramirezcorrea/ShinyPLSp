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
  if (is.data.frame(x)) return(as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE))
  if (is.matrix(x)) {
    if (is.numeric(x) && !is.null(rownames(x)) && !is.null(colnames(x))) {
      if (relation_type == "path") return(matrix_to_long(x, "Dependent", "Predictor", "Estimate", drop_zeros = TRUE))
      if (relation_type == "loading") return(matrix_to_long(x, "Construct", "Indicator", "Estimate", drop_zeros = TRUE))
    }
    return(as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE))
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
    htmt <- cbind(Construct = rownames(htmt), htmt, row.names = NULL)
  }

  fornell <- NULL
  if (!is.null(b$`Fornell-Larcker`)) {
    fornell <- as.data.frame(b$`Fornell-Larcker`, stringsAsFactors = FALSE, check.names = FALSE)
    fornell <- cbind(Construct = rownames(fornell), fornell, row.names = NULL)
  }

  list(quality = round_df(quality, 3), htmt = round_df(htmt, 3), fornell = round_df(fornell, 3))
}

extract_indirect_effects <- function(res) {
  if (is.null(res)) return(NULL)
  inf <- get_infer_safe(res)
  indirect_df <- NULL
  total_df    <- NULL
  if (!is.null(inf)) {
    for (nm in c("Indirect_effect", "Indirect_effects", "indirect_effect")) {
      if (!is.null(inf[[nm]])) {
        raw <- normalize_table_source(inf[[nm]], "path")
        if (!is.null(raw) && nrow(raw) > 0) { indirect_df <- raw; break }
      }
    }
    for (nm in c("Total_effect", "Total_effects", "total_effect")) {
      if (!is.null(inf[[nm]])) {
        raw <- normalize_table_source(inf[[nm]], "path")
        if (!is.null(raw) && nrow(raw) > 0) { total_df <- raw; break }
      }
    }
  }
  if (is.null(indirect_df) && !is.null(res$Estimates)) {
    est <- res$Estimates
    for (nm in c("Indirect_effect", "Indirect_effects", "indirect")) {
      if (!is.null(est[[nm]]) && (is.matrix(est[[nm]]) || is.data.frame(est[[nm]]))) {
        m <- est[[nm]]
        if (is.matrix(m)) {
          df_long <- matrix_to_long(m, "Dependent", "Predictor", "Indirect_Estimate", drop_zeros = TRUE)
          if (!is.null(df_long) && nrow(df_long) > 0) { indirect_df <- df_long; break }
        }
      }
    }
    for (nm in c("Total_effect", "Total_effects", "total")) {
      if (!is.null(est[[nm]]) && (is.matrix(est[[nm]]) || is.data.frame(est[[nm]]))) {
        m <- est[[nm]]
        if (is.matrix(m)) {
          df_long <- matrix_to_long(m, "Dependent", "Predictor", "Total_Estimate", drop_zeros = TRUE)
          if (!is.null(df_long) && nrow(df_long) > 0) { total_df <- df_long; break }
        }
      }
    }
  }
  parse_effects_table <- function(df, effect_label = "Indirect_Estimate") {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    nms     <- names(df)
    dep_col <- find_name(nms, c("^Dependent$","^Target$","dependent","target","endo"))
    pred_col<- find_name(nms, c("^Predictor$","^Origin$","predictor","origin","exo","independent"))
    rel_col <- find_name(nms, c("^Path$","^Relation$","^Parameter$","path","relation","parameter","Name"))
    est_col <- find_name(nms, c("^Estimate$","estimate","original","beta", effect_label))
    se_col  <- find_name(nms, c("^SE$","std.*error","standard.*error","std.*err"))
    t_col   <- find_name(nms, c("^t$","t.*value","t.*stat"))
    p_col   <- find_name(nms, c("^p$","p.*value"))
    cil_col <- find_name(nms, c("ci.*low","ci.*lower","lower.*ci","2\\.5","lower"))
    ciu_col <- find_name(nms, c("ci.*up","ci.*upper","upper.*ci","97\\.5","upper"))
    dep <- pred <- NULL
    if (!is.na(dep_col) && !is.na(pred_col)) {
      dep  <- as.character(df[[dep_col]])
      pred <- as.character(df[[pred_col]])
    } else if (!is.na(rel_col)) {
      pairs <- t(vapply(as.character(df[[rel_col]]), parse_path_relation, character(2)))
      dep <- pairs[,1]; pred <- pairs[,2]
    } else if (ncol(df) >= 2) {
      dep <- as.character(df[[1]]); pred <- as.character(df[[2]])
    } else return(NULL)
    if (is.na(est_col)) {
      num_candidates <- names(df)[sapply(df, is.numeric)]
      if (length(num_candidates) > 0) est_col <- num_candidates[1]
    }
    out <- data.frame(
      Dependent = dep,  Predictor = pred,
      Estimate  = if (!is.na(est_col))  to_num(df[[est_col]])  else NA_real_,
      Std_Error = if (!is.na(se_col))   to_num(df[[se_col]])   else NA_real_,
      T_value   = if (!is.na(t_col))    to_num(df[[t_col]])    else NA_real_,
      P_value   = if (!is.na(p_col))    to_num(df[[p_col]])    else NA_real_,
      CI_low    = if (!is.na(cil_col))  to_num(df[[cil_col]])  else NA_real_,
      CI_high   = if (!is.na(ciu_col))  to_num(df[[ciu_col]])  else NA_real_,
      stringsAsFactors = FALSE
    )
    keep <- !(is.na(out$Dependent) | is.na(out$Predictor) | out$Dependent == "" | out$Predictor == "")
    out  <- out[keep, , drop = FALSE]
    out$Significance <- significance_stars(out$P_value)
    rownames(out) <- NULL
    out
  }
  indirect_parsed <- parse_effects_table(indirect_df, "Indirect_Estimate")
  total_parsed    <- parse_effects_table(total_df,    "Total_Estimate")
  direct_df  <- extract_paths(res)
  all_pairs  <- unique(rbind(
    if (!is.null(direct_df)       && all(c("Predictor","Dependent") %in% names(direct_df)))       direct_df[,       c("Predictor","Dependent")] else NULL,
    if (!is.null(indirect_parsed) && all(c("Predictor","Dependent") %in% names(indirect_parsed))) indirect_parsed[, c("Predictor","Dependent")] else NULL,
    if (!is.null(total_parsed)    && all(c("Predictor","Dependent") %in% names(total_parsed)))    total_parsed[,    c("Predictor","Dependent")] else NULL
  ))
  if (is.null(all_pairs) || nrow(all_pairs) == 0) return(NULL)

  out <- all_pairs
  out$Relationship <- paste0(out$Predictor, " -> ", out$Dependent)
  if (!is.null(direct_df) && all(c("Predictor","Dependent","Estimate") %in% names(direct_df))) {
    d2 <- direct_df[, intersect(names(direct_df), c("Predictor","Dependent","Estimate","P_value","Significance"))]
    names(d2)[names(d2) == "Estimate"]     <- "Direct_Estimate"
    names(d2)[names(d2) == "P_value"]      <- "Direct_P_value"
    names(d2)[names(d2) == "Significance"] <- "Direct_Sig"
    out <- merge(out, d2, by = c("Predictor","Dependent"), all.x = TRUE)
  } else {
    out$Direct_Estimate <- NA_real_; out$Direct_P_value <- NA_real_; out$Direct_Sig <- NA_character_
  }
  if (!is.null(indirect_parsed) && nrow(indirect_parsed) > 0) {
    i2 <- indirect_parsed[, intersect(names(indirect_parsed), c("Predictor","Dependent","Estimate","P_value","Significance","CI_low","CI_high"))]
    names(i2)[names(i2) == "Estimate"]     <- "Indirect_Estimate"
    names(i2)[names(i2) == "P_value"]      <- "Indirect_P_value"
    names(i2)[names(i2) == "Significance"] <- "Indirect_Sig"
    names(i2)[names(i2) == "CI_low"]       <- "Indirect_CI_low"
    names(i2)[names(i2) == "CI_high"]      <- "Indirect_CI_high"
    out <- merge(out, i2, by = c("Predictor","Dependent"), all.x = TRUE)
  } else {
    out$Indirect_Estimate <- NA_real_; out$Indirect_P_value <- NA_real_
    out$Indirect_Sig      <- NA_character_
    out$Indirect_CI_low   <- NA_real_;  out$Indirect_CI_high  <- NA_real_
  }
  if (!is.null(total_parsed) && nrow(total_parsed) > 0) {
    t2 <- total_parsed[, intersect(names(total_parsed), c("Predictor","Dependent","Estimate","P_value","Significance"))]
    names(t2)[names(t2) == "Estimate"]     <- "Total_Estimate"
    names(t2)[names(t2) == "P_value"]      <- "Total_P_value"
    names(t2)[names(t2) == "Significance"] <- "Total_Sig"
    out <- merge(out, t2, by = c("Predictor","Dependent"), all.x = TRUE)
  } else {
    out$Total_Estimate <- ifelse(
      !is.na(out$Direct_Estimate) | !is.na(out$Indirect_Estimate),
      rowSums(cbind(ifelse(is.na(out$Direct_Estimate), 0, out$Direct_Estimate),
                    ifelse(is.na(out$Indirect_Estimate), 0, out$Indirect_Estimate))),
      NA_real_)
    out$Total_P_value <- NA_real_; out$Total_Sig <- NA_character_
  }

  cols_ord <- c("Relationship","Predictor","Dependent",
                "Direct_Estimate","Direct_P_value","Direct_Sig",
                "Indirect_Estimate","Indirect_P_value","Indirect_Sig",
                "Indirect_CI_low","Indirect_CI_high",
                "Total_Estimate","Total_P_value","Total_Sig")
  round_df(out[, intersect(cols_ord, names(out)), drop = FALSE], 3)
}

build_hypothesis_table <- function(res, relations_rv_df = NULL) {
  if (is.null(res)) return(NULL)
  paths_df <- extract_paths(res)
  if (is.null(paths_df) || nrow(paths_df) == 0) return(NULL)
  dep_col  <- if ("Dependent"  %in% names(paths_df)) "Dependent"  else NULL
  pred_col <- if ("Predictor"  %in% names(paths_df)) "Predictor"  else NULL
  if (is.null(dep_col) || is.null(pred_col)) return(NULL)
  out <- data.frame(
    Hypothesis   = paste0("H", seq_len(nrow(paths_df))),
    Relationship = paste0(paths_df[[pred_col]], " -> ", paths_df[[dep_col]]),
    Predictor    = paths_df[[pred_col]],
    Dependent    = paths_df[[dep_col]],
    stringsAsFactors = FALSE
  )

  est_col <- find_name(names(paths_df), c("^Estimate$","estimate","beta","original"))
  if (!is.na(est_col)) out$Beta <- to_num(paths_df[[est_col]])

  p_col <- find_name(names(paths_df), c("^P_value$","^p$","p.*value"))
  out$P_value <- if (!is.na(p_col)) to_num(paths_df[[p_col]]) else NA_real_
  t_col <- find_name(names(paths_df), c("^T_value$","^t$","t.*value","t.*stat"))
  if (!is.na(t_col)) out$T_value <- to_num(paths_df[[t_col]])
  cil_col <- find_name(names(paths_df), c("CI_low","ci.*low","lower"))
  ciu_col <- find_name(names(paths_df), c("CI_high","ci.*up","upper"))
  if (!is.na(cil_col)) out$CI_low  <- to_num(paths_df[[cil_col]])
  if (!is.na(ciu_col)) out$CI_high <- to_num(paths_df[[ciu_col]])
  out$Stars  <- significance_stars(out$P_value)
  out$Result <- ifelse(is.na(out$P_value), "No resampling",
                       ifelse(out$P_value <= 0.05, "Supported", "Not supported"))
    show_cols <- intersect(c("Hypothesis","Relationship","Beta","T_value","P_value",
                           "CI_low","CI_high","Stars","Result"), names(out))
  round_df(out[, show_cols, drop = FALSE], 3)
}

