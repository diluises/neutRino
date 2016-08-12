# Oscillation Probabilities for the 3 Flavour Model


### Load common parameters
data(db.pars)

A0 <- db.pars[db.pars$param == "A0", "value"]
rho <- db.pars[db.pars$param == "rho", "value"]

A <- A0*rho

####

require(compiler)



#' @export
#' @title
#' NuMuDisappearance
#' @details
#' 3 Flavour approximate formulta for NuMu->NuMu
#' @param E,L neutrino Energy and Baseline. Units: [GeV] and [km] or [MeV] and [m], ...
#' @param s2_12: sin^{2}_{12}
#' @param s2_13: sin^{2}_{13}
#' @param s2_23: sin^{2}_{23}
#' @param dm2_21: Delta_m^2_{21} [eV^2]
#' @param dm2_32: Delta_m^2_{32} [eV^2]
#' @param MH: Mass Hierachy, MH=1 for Normal Hierachy NH, MH=-1 for Inverted Hierarchy IH
#' @param dcp: delta CP phase
#' @return Probability
ProbNuMuToNuMu <- function(E, L, s2_12, s2_13, s2_23, dm2_21, dm2_32, MH, dcp) {

    #
    #
    #s2_12 <- pars[1]
    #s2_13 <- pars[2]
    #s2_23 <- pars[3]
    #dm2_21 <- pars[4]
    #dm2_32 <- pars[5]
    #MH <- pars[6]
    #dcp <- pars[7]
    #
    #

    s12 <- sqrt(s2_12)
    s13 <- sqrt(s2_13)
    s23 <- sqrt(s2_23)

    c12 <- sqrt(1 - s2_12)
    c13 <- sqrt(1 - s2_13)
    c23 <- sqrt(1 - s2_23)

    c2_12 <- (1 - s2_12)
    c2_13 <- (1 - s2_13)
    c2_23 <- (1 - s2_23)

    cosd <- cos(dcp)
    sind <- sin(dcp)

    dm2_31 <- dm2_21 + dm2_32

    if (MH == -1) {
        dm2_31 <- -dm2_31
        dm2_32 <- dm2_31 - dm2_21
    }

    phi21 <- 1.27 * dm2_21 * L/E
    phi31 <- 1.27 * dm2_31 * L/E
    phi32 <- 1.27 * dm2_32 * L/E

    s_phi31 <- sin(phi31)
    s_phi32 <- sin(phi32)
    s_phi21 <- sin(phi21)

    s2_phi31 <- s_phi31 * s_phi31
    s2_phi32 <- s_phi32 * s_phi32
    s2_phi21 <- s_phi21 * s_phi21

    c_phi31 <- cos(phi31)
    c_phi32 <- cos(phi32)
    c_phi21 <- cos(phi21)

    c2_phi31 <- c_phi31 * c_phi31
    c2_phi32 <- c_phi32 * c_phi32
    c2_phi21 <- c_phi21 * c_phi21

    P1 <- (s2_12 * c2_23 + s2_13 * s2_23 * c2_12 + 2 * s12 * s13 * s23 * c12 * c23 * cosd) * s2_23 * c2_13 * s2_phi31

    P2 <- (c2_12 * c2_23 + s2_13 * s2_23 * s2_12 - 2 * s12 * s13 * s23 * c12 * c23 * cosd) * s2_23 * c2_13 * s2_phi32

    P3 <- s2_12 * c2_23 + s2_13 * s2_23 * c2_12 + 2 * s12 * s13 * s23 * c12 * c23 * cosd

    P4 <- (c2_12 * c2_23 + s2_13 * s2_23 * s2_12 - 2 * s12 * s13 * s23 * c12 * c23 * cosd) * s2_phi21

    P <- 1 - 4 * P1 - 4 * P2 - 4 * P3 * P4

    return(P)
}


#' @export
#' @title
#' NuMuDisappearance compiled
#' Check ProbNuMuToNuMu for reference
#' Temporary keeping both
ProbNuMuToNuMu_cmp<- cmpfun(ProbNuMuToNuMu)





#' @export
#' @title
#' NuEAppearance
#' @details
#' 3 Flavour approximate formulta for NuMu->NuE
#' @param E,L neutrino Energy and Baseline. Units: [GeV] and [km] or [MeV] and [m], ...
#' @param s2_12: sin^{2}_{12}
#' @param s2_13: sin^{2}_{13}
#' @param s2_23: sin^{2}_{23}
#' @param dm2_21: Delta_m^2_{21} [eV^2]
#' @param dm2_32: Delta_m^2_{32} [eV^2]
#' @param MH: Mass Hierachy, MH=1 for Normal Hierachy NH, MH=-1 for Inverted Hierarchy IH
#' @param dcp: delta CP phase
#' @return Probability
ProbNuMuToNuE <- function(E, L, s2_12, s2_13, s2_23, dm2_21, dm2_32, MH, dcp ) {

    a <- A*E

    #
    #s2_12 <- pars[1]
    #s2_13 <- pars[2]
    #s2_23 <- pars[3]
    #dm2_21 <- pars[4]
    #dm2_32 <- pars[5]
    #MH <- pars[6]
    #dcp <- pars[7]
    #

    s12 <- sqrt(s2_12)
    s13 <- sqrt(s2_13)
    s23 <- sqrt(s2_23)

    c12 <- sqrt(1 - s2_12)
    c13 <- sqrt(1 - s2_13)
    c23 <- sqrt(1 - s2_23)

    c2_12 <- (1 - s2_12)
    c2_13 <- (1 - s2_13)
    c2_23 <- (1 - s2_23)

    cosd <- cos(dcp)
    sind <- sin(dcp)

    dm2_31 <- dm2_21 + dm2_32

    if (MH == -1) {
      dm2_31 <- -dm2_31
      dm2_32 <- dm2_31 - dm2_21
    }

    phi21 <- 1.27 * dm2_21 * L/E
    phi31 <- 1.27 * dm2_31 * L/E
    phi32 <- 1.27 * dm2_32 * L/E

    s_phi31 <- sin(phi31)
    s_phi32 <- sin(phi32)
    s_phi21 <- sin(phi21)

    s2_phi31 <- s_phi31 * s_phi31
    s2_phi32 <- s_phi32 * s_phi32
    s2_phi21 <- s_phi21 * s_phi21

    c_phi31 <- cos(phi31)
    c_phi32 <- cos(phi32)
    c_phi21 <- cos(phi21)

    c_phi23 <- c_phi32

    c2_phi31 <- c_phi31 * c_phi31
    c2_phi32 <- c_phi32 * c_phi32
    c2_phi21 <- c_phi21 * c_phi21


   P1 <- c2_13*s2_13*s2_23*s2_phi31*(1+(1-2*s2_13)*2*a/dm2_31)

   P2 <- c2_13*s12*s13*s23*(c12*c23*cosd - s12*s13*s23)*c_phi23*s_phi31*s_phi21

   P3 <- c2_13*c12*c23*s12*s13*s23*sind*s_phi32*s_phi31*s_phi21

   P4 <- s2_12*c2_13*(c2_12*c2_23+s2_12*s2_23*s2_13-2*c12*c23*s12*s23*s13*cosd)*s_phi21*s_phi21

   P5 <- c2_13*s2_13*s2_23*(1-2*s2_13)*(1.27*a*L/E)*c_phi32*s_phi31

   P <- 4*P1 + 8*P2 - 8*P3 +4*P4 -8*P5

   return(P)

}


#' @export
#' @title
#' NuEAppearance compiled
#' Check ProbNuMuToNuE for reference
#' Temporary keeping both
ProbNuMuToNuE_cmp<- cmpfun(ProbNuMuToNuMu)


