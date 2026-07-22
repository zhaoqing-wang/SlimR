#' Discrete colour palette from ArchR
#'
#' @description
#' Returns a set of colours for discrete categories, using the same
#' palettes as the **ArchR** package (e.g., *stallion*, *stallion2*, *calm*,
#' *kelly*, *bear*, *ironMan*, etc.).  When the number of categories exceeds
#' the palette length, colours are interpolated via
#' \code{\link[grDevices]{colorRampPalette}}.
#'
#' This function is a standalone reimplementation of
#' \code{ArchR::paletteDiscrete()} and is provided so that **SlimR** can
#' generate the same colour palettes without requiring the external
#' **ArchR** package.  **ArchR** is released under the MIT license.
#'
#' @param values A character or factor vector of categories to colour.
#' @param set The palette name, corresponding to one of the palettes
#'   defined in \code{ArchRPalettes}.  Default \code{"stallion"}.
#' @param reverse Logical indicating whether to reverse the order of the
#'   colours.  Default \code{FALSE}.
#'
#' @return A named character vector of hex colours, where names are the
#'   unique sorted categories from \code{values}.
#'
#' @section Palette sets:
#' Available palettes include \code{"stallion"}, \code{"stallion2"},
#' \code{"calm"}, \code{"kelly"}, \code{"bear"}, \code{"ironMan"},
#' \code{"circus"}, \code{"paired"}, \code{"grove"}, \code{"summerNight"},
#' \code{"zissou"}, \code{"darjeeling"}, \code{"rushmore"}, \code{"captain"},
#' \code{"horizon"}, \code{"horizonExtra"}, \code{"blueYellow"},
#' \code{"sambaNight"}, \code{"solarExtra"}, \code{"whitePurple"},
#' \code{"whiteBlue"}, \code{"whiteRed"}, \code{"comet"}, \code{"greenBlue"},
#' \code{"beach"}, \code{"coolwarm"}, \code{"fireworks"}, \code{"greyMagma"},
#' \code{"fireworks2"}, \code{"purpleOrange"}.
#'
#' @references
#' Granja JM, Corces MR et al. ArchR is a scalable software package for
#' integrative single-cell chromatin accessibility analysis.
#' *Nature Genetics* (2021).  \url{https://doi.org/10.1038/s41588-021-00790-6}
#'
#' \url{http://www.archrproject.com/}
#'
#' \url{https://github.com/GreenleafLab/ArchR}
#'
#' @family Section_5_Other_Functions_Provided
#' 
#' @importFrom gtools mixedsort
#' @importFrom grDevices colorRampPalette
#' @keywords internal
#' 
paletteDiscrete <- function(values = NULL, set = "stallion", reverse = FALSE) {
  values <- unique(values)
  values <- gtools::mixedsort(values)
  n <- length(values)
  pal <- ArchRPalettes[[set]]
  palOrdered <- pal[gtools::mixedsort(names(pal))]
  if (n > length(palOrdered)) {
    message("Length of unique values greater than palette, interpolating..")
    palOut <- colorRampPalette(pal)(n)
  } else {
    palOut <- palOrdered[seq_len(n)]
  }
  if (reverse) {
    palOut <- rev(palOut)
  }
  names(palOut) <- values
  palOut
}

#' ArchR colour palettes
#'
#' A list of colour vectors reproducing the palettes used by the
#' \pkg{ArchR} package.  Names correspond to the \code{set} argument of
#' \code{paletteDiscrete()}.
#'
#' @format A named list of character vectors.
#' @source \url{https://github.com/GreenleafLab/ArchR}
#' 
ArchRPalettes <- list(
  stallion = c(
    "1" = "#D51F26", "2" = "#272E6A", "3" = "#208A42", "4" = "#89288F",
    "5" = "#F47D2B", "6" = "#FEE500", "7" = "#8A9FD1", "8" = "#C06CAB",
    "19" = "#E6C2DC", "10" = "#90D5E4", "11" = "#89C75F", "12" = "#F37B7D",
    "13" = "#9983BD", "14" = "#D24B27", "15" = "#3BBCA8", "16" = "#6E4B9E",
    "17" = "#0C727C", "18" = "#7E1416", "9" = "#D8A767", "20" = "#3D3D3D"
  ),
  stallion2 = c(
    "1" = "#D51F26", "2" = "#272E6A", "3" = "#208A42", "4" = "#89288F",
    "5" = "#F47D2B", "6" = "#FEE500", "7" = "#8A9FD1", "8" = "#C06CAB",
    "19" = "#E6C2DC", "10" = "#90D5E4", "11" = "#89C75F", "12" = "#F37B7D",
    "13" = "#9983BD", "14" = "#D24B27", "15" = "#3BBCA8", "16" = "#6E4B9E",
    "17" = "#0C727C", "18" = "#7E1416", "9" = "#D8A767"
  ),
  calm = c(
    "1" = "#7DD06F", "2" = "#844081", "3" = "#688EC1", "4" = "#C17E73",
    "5" = "#484125", "6" = "#6CD3A7", "7" = "#597873", "8" = "#7B6FD0",
    "9" = "#CF4A31", "10" = "#D0CD47", "11" = "#722A2D", "12" = "#CBC594",
    "13" = "#D19EC4", "14" = "#5A7E36", "15" = "#D4477D", "16" = "#403552",
    "17" = "#76D73C", "18" = "#96CED5", "19" = "#CE54D1", "20" = "#C48736"
  ),
  kelly = c(
    "1" = "#FFB300", "2" = "#803E75", "3" = "#FF6800", "4" = "#A6BDD7",
    "5" = "#C10020", "6" = "#CEA262", "7" = "#817066", "8" = "#007D34",
    "9" = "#F6768E", "10" = "#00538A", "11" = "#FF7A5C", "12" = "#53377A",
    "13" = "#FF8E00", "14" = "#B32851", "15" = "#F4C800", "16" = "#7F180D",
    "17" = "#93AA00", "18" = "#593315", "19" = "#F13A13", "20" = "#232C16"
  ),
  bear = c(
    "1" = "#faa818", "2" = "#41a30d", "3" = "#fbdf72", "4" = "#367d7d",
    "5" = "#d33502", "6" = "#6ebcbc", "7" = "#37526d", "8" = "#916848",
    "9" = "#f5b390", "10" = "#342739", "11" = "#bed678", "12" = "#a6d9ee",
    "13" = "#0d74b6", "14" = "#60824f", "15" = "#725ca5", "16" = "#e0598b"
  ),
  ironMan = c(
    "9" = "#371377", "3" = "#7700FF", "2" = "#9E0142", "10" = "#FF0080",
    "14" = "#DC494C", "12" = "#F88D51", "1" = "#FAD510", "8" = "#FFFF5F",
    "4" = "#88CFA4", "13" = "#238B45", "5" = "#02401B", "7" = "#0AD7D3",
    "11" = "#046C9A", "6" = "#A2A475", "15" = "grey35"
  ),
  circus = c(
    "1" = "#D52126", "2" = "#88CCEE", "3" = "#FEE52C", "4" = "#117733",
    "5" = "#CC61B0", "6" = "#99C945", "7" = "#2F8AC4", "8" = "#332288",
    "9" = "#E68316", "10" = "#661101", "11" = "#F97B72", "12" = "#DDCC77",
    "13" = "#11A579", "14" = "#89288F", "15" = "#E73F74"
  ),
  paired = c(
    "9" = "#A6CDE2", "1" = "#1E78B4", "3" = "#74C476", "12" = "#34A047",
    "11" = "#F59899", "2" = "#E11E26", "10" = "#FCBF6E", "4" = "#F47E1F",
    "5" = "#CAB2D6", "8" = "#6A3E98", "6" = "#FAF39B", "7" = "#B15928"
  ),
  grove = c(
    "11" = "#1a1334", "9" = "#01545a", "1" = "#017351", "6" = "#03c383",
    "8" = "#aad962", "2" = "#fbbf45", "10" = "#ef6a32", "3" = "#ed0345",
    "7" = "#a12a5e", "5" = "#710162", "4" = "#3B9AB2"
  ),
  summerNight = c(
    "1" = "#2a7185", "2" = "#a64027", "3" = "#fbdf72", "4" = "#60824f",
    "5" = "#9cdff0", "6" = "#022336", "7" = "#725ca5"
  ),
  zissou = c(
    "1" = "#3B9AB2", "4" = "#78B7C5", "3" = "#EBCC2A", "5" = "#E1AF00", "2" = "#F21A00"
  ),
  darjeeling = c(
    "1" = "#FF0000", "2" = "#00A08A", "3" = "#F2AD00", "4" = "#F98400", "5" = "#5BBCD6"
  ),
  rushmore = c(
    "1" = "#E1BD6D", "5" = "#EABE94", "2" = "#0B775E", "4" = "#35274A", "3" = "#F2300F"
  ),
  captain = c(
    "1" = "grey", "2" = "#A1CDE1", "3" = "#12477C", "4" = "#EC9274", "5" = "#67001E"
  ),
  horizon = c(
    "1" = "#000075", "4" = "#2E00FF", "6" = "#9408F7", "10" = "#C729D6",
    "8" = "#FA4AB5", "3" = "#FF6A95", "7" = "#FF8B74", "5" = "#FFAC53",
    "9" = "#FFCD32", "2" = "#FFFF60"
  ),
  horizonExtra = c(
    "1" = "#000436", "4" = "#021EA9", "6" = "#1632FB", "8" = "#6E34FC",
    "3" = "#C732D5", "9" = "#FD619D", "7" = "#FF9965", "5" = "#FFD32B",
    "2" = "#FFFC5A"
  ),
  blueYellow = c(
    "1" = "#352A86", "2" = "#343DAE", "3" = "#0262E0", "4" = "#1389D2",
    "5" = "#2DB7A3", "6" = "#A5BE6A", "7" = "#F8BA43", "8" = "#F6DA23", "9" = "#F8FA0D"
  ),
  sambaNight = c(
    "6" = "#1873CC", "2" = "#1798E5", "8" = "#00BFFF", "5" = "#4AC596",
    "1" = "#00CC00", "4" = "#A2E700", "9" = "#FFFF00", "7" = "#FFD200", "3" = "#FFA500"
  ),
  solarExtra = c(
    "5" = "#3361A5", "7" = "#248AF3", "1" = "#14B3FF", "8" = "#88CEEF",
    "9" = "#C1D5DC", "4" = "#EAD397", "3" = "#FDB31A", "2" = "#E42A2A", "6" = "#A31D1D"
  ),
  whitePurple = c(
    "9" = "#f7fcfd", "6" = "#e0ecf4", "8" = "#bfd3e6", "5" = "#9ebcda",
    "2" = "#8c96c6", "4" = "#8c6bb1", "7" = "#88419d", "3" = "#810f7c", "1" = "#4d004b"
  ),
  whiteBlue = c(
    "9" = "#fff7fb", "6" = "#ece7f2", "8" = "#d0d1e6", "5" = "#a6bddb",
    "2" = "#74a9cf", "4" = "#3690c0", "7" = "#0570b0", "3" = "#045a8d", "1" = "#023858"
  ),
  whiteRed = c("1" = "white", "2" = "red"),
  comet = c("1" = "#E6E7E8", "2" = "#3A97FF", "3" = "#8816A7", "4" = "black"),
  greenBlue = c(
    "4" = "#e0f3db", "7" = "#ccebc5", "2" = "#a8ddb5", "5" = "#4eb3d3",
    "3" = "#2b8cbe", "6" = "#0868ac", "1" = "#084081"
  ),
  beach = c(
    "4" = "#87D2DB", "1" = "#5BB1CB", "6" = "#4F66AF", "3" = "#F15F30",
    "5" = "#F7962E", "2" = "#FCEE2B"
  ),
  coolwarm = c(
    "1" = "#4858A7", "4" = "#788FC8", "5" = "#D6DAE1", "3" = "#F49B7C", "2" = "#B51F29"
  ),
  fireworks = c(
    "5" = "white", "2" = "#2488F0", "4" = "#7F3F98", "3" = "#E22929", "1" = "#FCB31A"
  ),
  greyMagma = c(
    "2" = "grey", "4" = "#FB8861FF", "5" = "#B63679FF", "3" = "#51127CFF", "1" = "#000004FF"
  ),
  fireworks2 = c(
    "5" = "black", "2" = "#2488F0", "4" = "#7F3F98", "3" = "#E22929", "1" = "#FCB31A"
  ),
  purpleOrange = c(
    "5" = "#581845", "2" = "#900C3F", "4" = "#C70039", "3" = "#FF5744", "1" = "#FFC30F"
  )
)