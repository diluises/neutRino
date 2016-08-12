# Utility function for neutrino physics

# Retrive constants
data(db.pars)

Eb <- db.pars[db.pars$param == "Eb", "value"]
me <- db.pars[db.pars$param == "mass_ele", "value"]
mmu <- db.pars[db.pars$param == "mass_muon", "value"]
Mp <- db.pars[db.pars$param == "mass_proton", "value"]
Mn <- db.pars[db.pars$param == "mass_neutron", "value"]


Mp2 <- Mp^2
me2 <- me^2
mmu2 <- mmu^2

#

#' @export
#' @title Nu Reconstructed Energy
#' @description Reconstructed Neutrino Energy based on quasielastic kinematic
#' @param Mlepton: Mass of the outgoing lepton [MeV] \cr
#' @param Plepton: Momentum of the outgoing lepton  [MeV/c] \cr
#' @param cos.theta: theta = angle of the outgoing lepton  wrt to the incoming neutrino \cr
Nu.Ereco <- function(Mlepton, Plepton, cos.theta) {

    Elepton <- sqrt(Mlepton^2 + Plepton^2)

    Ereco <- (Mp2 - (Mn - Eb)^2 - Mlepton^2 + 2 * (Mn - Eb) * Elepton)

    Ereco <- Ereco/2/(Mn - Eb - Elepton + Plepton * cos.theta)

}



#' @export
#' @title NuE Reconstructed Energy
#' @description Reconstructed Neutrino Energy based on quasielastic kinematic
#' electron neutrino hyp.
#' @param Elepton: Energy of the outgoing lepton (electron) [MeV]\cr
#' @param Plepton: Momentum of the outgoing lepton (electron) [MeV/c]\cr
#' @param cos.theta: theta = angle of the outgoing lepton (electron) wrt to the incoming neutrino\cr
Nu.ErecoAsEle <- function(Plepton, cos.theta) {

    Elepton <- sqrt(me2 + Plepton^2)

    Ereco <- (Mp2 - (Mn - Eb)^2 - me2 + 2 * (Mn - Eb) * Elepton)

    Ereco <- Ereco/2/(Mn - Eb - Elepton + Plepton * cos.theta)

}


#' @export
#' @title NuMu Reconstructed Energy
#' @description Reconstructed Neutrino Energy based on quasielastic kinematic
#' muon neutrino hyp.
#' @param Elepton: Energy of the outgoing lepton (muon) [MeV] \cr
#' @param Plepton: Momentum of the outgoing lepton (muon) [MeV]\cr
#' @param cos.theta: theta = angle of the outgoing lepton (muon) wrt to the incoming neutrino\cr
Nu.ErecoAsMuon <- function(Plepton, cos.theta) {

    Elepton <- sqrt(mmu2 + Plepton^2)

    Ereco <- (Mp2 - (Mn - Eb)^2 - me2 + 2 * (Mn - Eb) * Elepton)

    Ereco <- Ereco/2/(Mn - Eb - Elepton + Plepton * cos.theta)

}

