# ==============================================================================
# HELPERS DE MODELO Y GRAFICOS
# ==============================================================================

get_model_syntax <- function(clist, df_rel) {
  if (length(clist) == 0) return("No constructs defined yet.")
  
  meas_lines <- unlist(lapply(names(clist), function(nm) {
    items <- clist[[nm]]$items
    typ <- clist[[nm]]$type
    op <- if (identical(typ, "Composite")) "<~" else "=~"
    paste0(nm, " ", op, " ", paste(items, collapse = " + "))
  }))
  
  struct_lines <- character(0)
  if (nrow(df_rel) > 0) {
    targets <- unique(df_rel$target)
    struct_lines <- sapply(targets, function(tg) {
      orig <- df_rel$origin[df_rel$target == tg]
      paste0(tg, " ~ ", paste(orig, collapse = " + "))
    })
  }
  paste(c(struct_lines, meas_lines), collapse = "\n")
}

has_cycle <- function(df_rel) {
  if (is.null(df_rel) || nrow(df_rel) == 0) return(FALSE)
  nodes <- unique(c(df_rel$origin, df_rel$target))
  adj <- split(df_rel$target, df_rel$origin)
  adj <- lapply(adj, unique)
  visited <- setNames(rep(FALSE, length(nodes)), nodes)
  stack <- setNames(rep(FALSE, length(nodes)), nodes)
  
  dfs <- function(v) {
    visited[[v]] <<- TRUE
    stack[[v]] <<- TRUE
    neigh <- adj[[v]]
    if (!is.null(neigh)) {
      for (u in neigh) {
        if (!isTRUE(visited[[u]])) {
          if (dfs(u)) return(TRUE)
        } else if (isTRUE(stack[[u]])) {
          return(TRUE)
        }
      }
    }
    stack[[v]] <<- FALSE
    FALSE
  }
  
  for (v in nodes) {
    if (!isTRUE(visited[[v]])) {
      if (dfs(v)) return(TRUE)
    }
  }
  FALSE
}

edge_exists <- function(df_rel, origin, target) {
  if (is.null(df_rel) || nrow(df_rel) == 0) return(FALSE)
  any(df_rel$origin == origin & df_rel$target == target)
}

would_create_cycle <- function(df_rel, origin, target) {
  new_df <- rbind(df_rel, data.frame(origin = origin, target = target, stringsAsFactors = FALSE))
  has_cycle(new_df)
}

plot_graph <- function(clist, df_rel, res = NULL, show_items = TRUE, estimated = FALSE) {
  if (is.null(clist) || length(clist) == 0) {
    return(DiagrammeR::grViz("
      digraph sem {
        graph [layout = dot, rankdir = LR]
        node [shape = box, style = filled, color = lightgray]
        empty [label = 'Define at least one construct']
      }
    "))
  }
  
  r2_map <- NULL
  if (isTRUE(estimated) && !is.null(res)) {
    r2_df <- try(extract_r2(res), silent = TRUE)
    if (!inherits(r2_df, "try-error") && !is.null(r2_df) &&
        all(c("Construct", "R_squared") %in% names(r2_df))) {
      vals <- suppressWarnings(as.numeric(r2_df$R_squared))
      names(vals) <- as.character(r2_df$Construct)
      r2_map <- vals
    }
  }
  
  endo <- character(0)
  if (!is.null(df_rel) && nrow(df_rel) > 0) {
    if ("target" %in% names(df_rel)) {
      endo <- unique(as.character(df_rel$target))
    }
  }
  
  cnames <- names(clist)
  if (is.null(cnames) || length(cnames) == 0) {
    return(DiagrammeR::grViz("
      digraph sem {
        graph [layout = dot, rankdir = LR]
        node [shape = box, style = filled, color = lightgray]
        empty [label = 'Define at least one construct']
      }
    "))
  }
  
  latent_nodes <- paste(
    vapply(cnames, function(nm) {
      x <- clist[[nm]]
      x_type <- if (!is.null(x$type)) x$type else "Common factor"
      node_shape <- if (identical(x_type, "Composite")) "hexagon" else "ellipse"
      fill <- "white"
      border <- "#444444"
      label <- nm
      
      if (isTRUE(estimated) && !is.null(r2_map) &&
          nm %in% names(r2_map) && nm %in% endo &&
          !is.na(r2_map[[nm]])) {
        label <- paste0(nm, "\\nR2 = ", sprintf("%.3f", r2_map[[nm]]))
      }
      
      paste0(
        '"', nm, '" [shape=', node_shape,
        ', style="filled"',
        ', color="', border, '"',
        ', fillcolor="', fill, '"',
        ', penwidth=1.2',
        ', fontsize=10',
        ', width=1.0',
        ', height=0.8',
        ', label="', label, '"];'
      )
    }, character(1)),
    collapse = "\n"
  )
  
  item_nodes <- ""
  meas_edges <- ""
  
  if (isTRUE(show_items)) {
    all_items <- unique(unlist(lapply(clist, function(x) {
      if (!is.null(x$items)) x$items else character(0)
    }), use.names = FALSE))
    
    if (length(all_items) > 0) {
      item_nodes <- paste(
        vapply(all_items, function(it) {
          paste0(
            '"', it, '" [shape=box, style="rounded,filled", fillcolor="white", color="#777777", penwidth=1, fontsize=9];'
          )
        }, character(1)),
        collapse = "\n"
      )
    }
    
    meas_list <- lapply(cnames, function(nm) {
      x <- clist[[nm]]
      its <- if (!is.null(x$items)) x$items else character(0)
      typ <- if (!is.null(x$type)) x$type else "Common factor"
      if (length(its) == 0) return(NULL)
      data.frame(from = nm, to = its, type = typ, stringsAsFactors = FALSE)
    })
    
    meas_df <- do.call(rbind, meas_list)
    
    if (!is.null(meas_df) && nrow(meas_df) > 0) {
      meas_edges <- paste(
        apply(meas_df, 1, function(r) {
          edge_style <- if (r["type"] == "Composite") {
            'style="dashed", color="#7A7A7A", penwidth=1.1'
          } else {
            'style="solid", color="#7A7A7A", penwidth=1.1'
          }
          paste0('"', r["from"], '" -> "', r["to"], '" [', edge_style, ', arrowsize=0.7];')
        }),
        collapse = "\n"
      )
    }
  }
  
  struct_edges <- ""
  if (!is.null(df_rel) && nrow(df_rel) > 0) {
    from_col <- if ("origin" %in% names(df_rel)) "origin" else if ("from" %in% names(df_rel)) "from" else NULL
    to_col   <- if ("target" %in% names(df_rel)) "target" else if ("to" %in% names(df_rel)) "to" else NULL
    
    if (!is.null(from_col) && !is.null(to_col)) {
      struct_edges <- paste(
        apply(df_rel, 1, function(r) {
          paste0(
            '"', r[[from_col]], '" -> "', r[[to_col]],
            '" [color="#222222", penwidth=1.5, arrowsize=0.8];'
          )
        }),
        collapse = "\n"
      )
    }
  }
  
  dot_code <- paste0(
    "digraph sem {\n",
    'graph [layout=dot, rankdir=LR, nodesep=0.35, ranksep=0.5, splines=true, overlap=false, bgcolor="transparent"];\n',
    'node [fontname=Helvetica, fontsize=10];\n',
    'edge [fontname=Helvetica, fontsize=9, labelfloat=true];\n\n',
    latent_nodes, "\n\n",
    item_nodes, "\n\n",
    meas_edges, "\n\n",
    struct_edges, "\n",
    "}\n"
  )
  
  DiagrammeR::grViz(dot_code)
}