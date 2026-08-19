project_markers <- c(
  "config/simulation.R",
  "R/simulate.R",
  "README.qmd"
)
if (!all(file.exists(project_markers))) {
  stop("Run this command from the project root.")
}

generated_directories <- c(
  "data/derived",
  "reports/_freeze"
)
generated_artifacts <- list.files(
  "artifacts",
  pattern = "[.](pt|pth|rds|csv)$",
  full.names = TRUE
)

existing_directories <- generated_directories[dir.exists(generated_directories)]
existing_artifacts <- generated_artifacts[file.exists(generated_artifacts)]

if (length(existing_directories)) {
  unlink(existing_directories, recursive = TRUE, force = TRUE)
}
if (length(existing_artifacts)) {
  unlink(existing_artifacts, force = TRUE)
}

cat(
  "Removed generated directories:",
  if (length(existing_directories)) {
    paste(existing_directories, collapse = ", ")
  } else {
    "none"
  },
  "\n"
)
cat(
  "Removed generated model artifacts:",
  if (length(existing_artifacts)) {
    paste(existing_artifacts, collapse = ", ")
  } else {
    "none"
  },
  "\n"
)
