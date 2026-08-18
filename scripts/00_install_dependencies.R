required <- c("epiworldR", "torch")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)

if (!torch::torch_is_installed()) {
  torch::install_torch()
  message("Torch runtime installed. Restart R before training if requested.")
}

