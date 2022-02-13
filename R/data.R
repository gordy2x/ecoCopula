#' Spider data
#'
#' Abundance of hunting spiders and associated environmental variables 
#'
#' The matrix \code{abund} has the following species abundances (column name abbreviation in brackets)
#' \itemize{
#'   \item Alopecosa accentuata (Alopacce)
#'   \item Alopecosa cuneata (Alopcune)
#'   \item Alopecosa fabrilis (Alopfabr)
#'   \item Arctosa lutetiana (Arctlute)
#'   \item Arctosa perita(Arctperi)
#'   \item Aulonia albimana (Auloalbi)
#'   \item Pardosa lugubris (Pardlugu)
#'   \item Pardosa monticola (Pardmont)
#'   \item Pardosa nigriceps (Pardnigr)
#'   \item Pardosa pullata (Pardpull)
#'   \item Trochosa terricola (Trocterr)
#'   \item Zora spinimana (Zoraspin)
#' }
#' The data frame \code{x} has the following log(x+1)-transformed environmental variables 
#' \itemize{
#'   \item soil.dry - Soil dry mass
#'   \item bare.sand - Cover bare sand
#'   \item fallen.leaves - Cover fallen leaves / twigs
#'   \item moss - Cover moss
#'   \item herb.layer - Cover herb layer
#'   \item reflection - Reflection of the soil surface with a cloudless sky
#' }
#'  The data frame \code{trait} has the following variables
#' \itemize{
#'   \item length (numeric) - Length (log-transformed), averaged across typical lengths (in centimetres) for male and females
#'   \item (factor) - Predominant colour, "yellow" or "dark"
#'   \item (factor) - Whether the spider typically has markings on it: "none", "spots" or "stripes"
#' }
#' 
#' @format A list containing the elements
#' \describe{
#'   \item{abund}{A matrix with 28 observations of abundance of 12 hunting spider species.}
#'   \item{x}{A data frame of six (transformed) environmental variables at each of the 28 sites.}
#'   \item{trait}{A data frame of three species trait variables for each of the 12 species.}
#' }
#' @source Data attributed to van der Aart & Smeenk-Enserink (1975), 
#' obtained from the spider2 directory, CANOCO FORTRAN package, with trait data added by David Warton,
#' exported from mvabund R package.
"spider"