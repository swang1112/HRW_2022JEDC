#' Simulate Data from given DGP
#'
#' @param DGP Integer, indicates DGP
#' @param Config List, configuration for simulation
#' @param ... Other parameters for DGP
#'
#' @return A list containing simulated dats.
#' @export simu.datas
#'
#' @examples NULL
simu.datas <- function(Config, DGP, ...){
  Y          <- list()
  start.time <- Sys.time()
  Size       <- Config$Size

  if (DGP == 1){
    for (i in 1:Size) {
      Y[[i]] <- DGP1(Par = Config$Par, Tob = Config$Tob, Break = Config$Break,
                     Volaregim = Config$Volaregim, VolaRand = Config$VolaRand,
                     Edist = Config$Edist, v = Config$v, Combi = Config$Combi,
                     ME = Config$ME, MERand = Config$MERand)
      progress(x = i, max = Size, start.time = start.time)
    }
  } else if (DGP == 2){
    for (i in 1:Size) {
      Y[[i]] <- DGP2(Par = Config$Par, Tob = Config$Tob, Break = Config$Break,
                     Volaregim = Config$Volaregim, VolaRand = Config$VolaRand,
                     Edist = Config$Edist, v = Config$v, Weak = Config$Weak,
                     Noise = Config$Noise, NoiseRand = Config$NoiseRand,
                     Combi = Config$Combi, ME = Config$ME, MERand = Config$MERand)
      progress(x = i, max = Size, start.time = start.time)
    }
  } else {stop("Undefined DGP!")}

  return(Y)
}
