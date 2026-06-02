# Ejecutar desde la raíz del proyecto

root_dir <- "."
output_file <- "proyecto_consolidado_para_IA2.txt"

exclude_dirs <- c(
  "^renv($|/)",
  "^\\.git($|/)",
  "^rsconnect($|/)",
  "^packrat($|/)",
  "^node_modules($|/)",
  "^www($|/)",
  "^dist($|/)",
  "^docs($|/)",
  "^output($|/)"
)

all_r_files <- list.files(
  path = root_dir,
  pattern = "\\.[Rr]$",
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

rel_paths <- gsub("^\\./", "", all_r_files)

keep_file <- function(x) {
  !any(sapply(exclude_dirs, function(p) grepl(p, x)))
}

# Primero filtras por directorios
r_files <- all_r_files[keep_file(rel_paths)]
rel_paths <- rel_paths[keep_file(rel_paths)]

# Ahora excluyes específicamente join-tool.R
exclude_files <- c("join-tool.R")

keep_by_name <- !basename(r_files) %in% exclude_files
r_files <- r_files[keep_by_name]
rel_paths <- rel_paths[keep_by_name]

header <- c(
  paste0("PROYECTO: ", basename(normalizePath(root_dir))),
  "",
  "ESTRUCTURA DE ARCHIVOS INCLUIDOS:",
  paste0("- ", rel_paths),
  "",
  "========================================",
  "INICIO DEL CODIGO CONSOLIDADO",
  "========================================",
  ""
)

writeLines(header, output_file, useBytes = TRUE)

for (i in seq_along(r_files)) {
  block_header <- c(
    "",
    "========================================",
    paste0("ARCHIVO: ", basename(r_files[i])),
    paste0("RUTA: ", rel_paths[i]),
    "========================================"
  )
  
  write(block_header, file = output_file, append = TRUE)
  
  lines <- readLines(r_files[i], warn = FALSE, encoding = "UTF-8")
  
  if (length(lines) == 0) {
    lines <- c("[archivo vacio]")
  }
  
  write(lines, file = output_file, append = TRUE)
}

cat("Archivo generado:", output_file, "\n")