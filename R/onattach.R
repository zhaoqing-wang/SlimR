.onAttach <- function(libname, pkgname) {
  if (requireNamespace("crayon", quietly = TRUE)) {
    message <- paste0(
      "Please cite: Wang Z (2026). ",
      crayon::italic("SlimR: Adaptive Machine Learning-Powered, Context-Matching Tool for Single-Cell and Spatial Transcriptomics Annotation."),
      " R package version", crayon::bold(" 1.1.1."),
      " Available at: https://github.com/zhaoqing-wang/SlimR"
    )
  } else {
    message <- paste0(
      "Please cite: Wang Z (2026). SlimR: Adaptive Machine Learning-Powered, Context-Matching Tool for Single-Cell and Spatial Transcriptomics Annotation.",
      "R package version 1.1.1. Available at: https://github.com/zhaoqing-wang/SlimR"
    )
  }

  packageStartupMessage(message)
}
