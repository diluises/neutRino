##' Oscillation Probability for oscillations mediated
##' by sterilne neutrinos
#'
#'
#'


#' @export
#' @title NuE Disppearance via Sterile Transition
#' @details
#' NuE disappearance
#' in the Model 3+1.\cr
#' Approximate expression
#' @param E,L: neutrino Energy and Baseline. Units: [GeV] and [km] or [MeV] and [m], ...
#' @param s2_ee: sin^{2}_{ee}
#' @param dm2_41: in  Delta_m^2_{41} [eV^2]
#' @return probability
Sterile.ProbNuEToNuEv2 <- function(E, L, pars = c(s2_2ee, dm2_41)) {

    p <- as.list(pars)

    s2 <- p$s2_2ee
    dm2 <- p$dm2_41

    arg <- 1.27 * dm2 * L/E

    P <- 1 - s2 * sin(arg) * sin(arg)

    invisible(P)
}

#' @export
#' @title NuE Disppearance via Sterile Transition
#' @details
#' NuE disappearance
#' in the Model 3+1.\cr
#' Approximate expression
#' @param E,L: neutrino Energy and Baseline. Units: [GeV] and [km] or [MeV] and [m], ...
#' @param s2_ee: sin^{2}_{ee}
#' @param dm2_41: in  Delta_m^2_{41} [eV^2]
#' @return probability
Sterile.ProbNuEToNuE <- function(E, L, s2_2ee, dm2_41) {

  arg <- 1.27 * dm2_41 * L/E

  P <- 1 - s2_2ee * sin(arg) * sin(arg)

  invisible(P)
}



