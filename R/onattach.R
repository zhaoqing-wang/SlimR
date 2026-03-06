.onAttach <- function(libname, pkgname) {
  pkg_version <- utils::packageVersion("SlimR")

  if (requireNamespace("crayon", quietly = TRUE)) {
    msg <- paste0(
      "Please cite: Wang Z (2026). ",
      crayon::italic(
        "SlimR: Adaptive Machine Learning-Powered, Context-Matching Tool for Single-Cell and Spatial Transcriptomics Annotation."
      ),
      " R package version",
      crayon::bold(paste0(" ", pkg_version, ".")),
      " Available at: https://github.com/zhaoqing-wang/SlimR"
    )
  } else {
    msg <- paste0(
      "Please cite: Wang Z (2026). SlimR: Adaptive Machine Learning-Powered, Context-Matching Tool for Single-Cell and Spatial Transcriptomics Annotation.",
      "R package version ",
      pkg_version,
      ". Available at: https://github.com/zhaoqing-wang/SlimR"
    )
  }

  packageStartupMessage(msg)
}
