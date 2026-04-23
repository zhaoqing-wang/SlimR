.onAttach <- function(libname, pkgname) {
  pkg_version <- utils::packageVersion("SlimR")

  banner <- paste(
    c(
      "  ____  _ _           ____   ",
      " / ___|| (_)_ __ ___ |  _ \\  ",
      " \\___ \\| | | '_ ` _ \\| |_) | ",
      "  ___) | | | | | | | |  _ <  ",
      " |____/|_|_|_| |_| |_|_| \\_\\ "
    ),
    collapse = "\n"
  )

  title <- "SlimR: Adaptive Machine Learning-Powered, Context-Matching Tool for Single-Cell and Spatial Transcriptomics Annotation."
  url <- "https://github.com/zhaoqing-wang/SlimR"
  use_crayon <- requireNamespace("crayon", quietly = TRUE)

  if (use_crayon) {
    citation <- paste0(
      "Please cite: Wang Z (2026). ",
      crayon::italic(title),
      " R package version",
      crayon::bold(paste0(" ", pkg_version, ".")),
      " Available at: ",
      url
    )
    msg <- paste0(crayon::cyan$bold(banner), "\n", citation)
  } else {
    citation <- paste0(
      "Please cite: Wang Z (2026). ",
      title,
      " R package version ",
      pkg_version,
      ". Available at: ",
      url
    )
    msg <- paste0(banner, "\n", citation)
  }

  packageStartupMessage(msg)
}
