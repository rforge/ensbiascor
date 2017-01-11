# Documentation of data files in ensbiascoR
# -----------------------------------------

# Sebastian Sippel
# 10.10.2016

#' @name ensbiascoR.example1
#' @title A test dataset to run example1 on http://ensbiascor.r-forge.r-project.org/ensbiascoR_example1.html
#' @docType data
#' @usage data(ensbiascoR.example1)
#' @format A list with two entries:
#' \describe{
#'   \item{obs.data}{Data frame with reanalysis data from ERA-Interim, aggregated over Central Europe (JJA), for the years 1979-2014.}
#'   \item{mod.data}{Data frame with HadRM3P simulated data, aggregated over Central Europe (JJA), for the years 1985-2010 and 800 ensemble members for each year.}
#'   ...
#' }
#' @keywords datasets
#' @examples 
#' data(ensbiascoR.example1)
#' @source Sippel, S., Otto, F. E. L., Forkel, M., Allen, M. R., Guillod, B. P., Heimann, M., Reichstein, M., Seneviratne, S. I., Thonicke, K., and Mahecha, M. D. (2016) A novel bias correction methodology for climate impact simulations. Earth System Dynamics, 7, 71-88. doi:10.5194/esd-7-71-2016.
"ensbiascoR.example1"


#' @name ensbiascoR.example2
#' @title A test dataset to run example2 on http://ensbiascor.r-forge.r-project.org/ensbiascoR_example2.html
#' @docType data
#' @usage data(ensbiascoR.example2)
#' @format A list with three locations: AT-Vienna (AT_WIE), DE-Jena (DE_JEN), and NL-De Bilt (NL_DEB). Each location contains a list with two entries described below:
#' \describe{
#'   \item{obs.data}{Data frame with station data (from the ECA&D dataset), processed to seasonal maxima of indiviual heat-health related metrics}
#'   \item{mod.data}{Data frame with HadRM3P simulated data, for the years 1985-2010, ~20.000 ensemble members in total, and processed to seasonal maxima of indiviual heat-health related metrics}
#'   ...
#' }
#' @keywords datasets
#' @examples 
#' data(ensbiascoR.example2)
#' @source Sippel, S., Otto, F. E. L., Flach, M., and van Oldenborgh, G. J. (2016). The Role of  Anthropogenic Warming in 2015 Centra European Heat Waves. In Herring, S. C., Hoell, A., Hoerling, M. P., Kossin, J. P., Schreck III, C. J., and Stott, P. A. (Eds.), Explaining Extremes of 2015 from a Climate Perspective. Bull. Amer. Meteor. Soc., 97(12), S51–S56. doi:10.1175/BAMS-D-16-0149.
"ensbiascoR.example2"

