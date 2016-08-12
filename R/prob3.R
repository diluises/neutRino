## standard-R implementation of the \code{prob3} C-class
## \code{prob3} implements 3-flavour neutrino oscillation probability
## as developed in the paper of
## The "C"-version of \code{prob3} can be found here

#
#' @import compiler
#' @import parallel
#' @import doParallel
#' @import foreach


#
# Define global variables
#

matrixTypes <- as.factor( c("standard","barger") )
nuTypes     <- as.factor( c("data","nue","numu","nutau","sterile","unknown") )


#' @export
#' @title environment for global variables
#' @description main functions \code{\link{setMNS}},\code{\link{computeProb}},\code{\link{computeVaccumProb}}
my.env <- new.env(parent=emptyenv())

my.env$matrixTypes <- matrixTypes
my.env$nuTypes     <- nuTypes


#
##
#
my.env$matrixType <- matrixTypes[1]
my.env$nuType     <- nuTypes[2]

#
## Init Mixing (MNS) matrix (Re and Im part) and mass matrix
#
my.env$Mix <- array(0,dim=c(3,3,2))
my.env$dm  <- array(0,dim=c(3,3))

#
## Flag for dcp=0 case
#
my.env$zero_cp <- FALSE


#
##
#
my.env$dmVacVac <- array(0,dim=c(3,3))
my.env$dmMatMat <- array(0,dim=c(3,3))
my.env$dmMatVac <- array(0,dim=c(3,3))

#
## Amplitude
#
my.env$A     <- array(0,dim=c(3,3,2))

#
## Calculated transition probabilities between the 3 falvour states, matrix is filled after a
## a call to "propgate" or "compute" functions.
#
my.env$Prob  <- array(0,dim=c(3,3))

#
## Input parameters
#
my.env$dm21   <- 0
my.env$dm32   <- 0
my.env$dm31   <- 0
my.env$dmAtm  <- 0
my.env$dmSol  <- 0
my.env$MH     <- 1 #(NH)

my.env$s12   <- 0
my.env$s23   <- 0
my.env$s31   <- 0
my.env$delta <- 0
my.env$rho   <- 0
my.env$isAnti<-FALSE
#
## Flavor Indexes
#
my.env$elec <- 1
my.env$muon <- 2
my.env$tau  <- 3

#
#
## Dump input parameters
#
my.env$dump.params <- function() {

  cat("----prob3 input parameters --- \n")
  cat(" sin_12: ",my.env$s12,"\n")
  cat(" sin_23: ",my.env$s23,"\n")
  cat(" sin_31: ",my.env$s31,"\n")
  cat("\n")
  cat(" sin2_12: ",(my.env$s12)^2,"\n")
  cat(" sin2_23: ",(my.env$s23)^2,"\n")
  cat(" sin2_31: ",(my.env$s31)^2,"\n")
  cat("\n")
  cat(" MH: ",(my.env$MH),"\n")
  cat(" delta: ",(my.env$delta),"\n")
  cat(" is anti: ",(my.env$isAnti),"\n")
  cat("\n")
  cat(" dmAtm: ",my.env$dmAtm," [ev^2/c^4]\n")
  cat(" dmSol: ",my.env$dmSol," [ev^2/c^4]\n")
  cat("\n")
  cat(" dm21: ",my.env$dm21," [eV^2/c^4]\n")
  cat(" dm31: ",my.env$dm31," [ev^2/c^4]\n")
  cat(" dm32: ",my.env$dm32," [ev^2/c^4]\n")
  cat(" density ",my.env$rho," g/cm3\n")
  cat("----\n")

}

#
#
#

#
## test function for internal debugging
#
prob3.test <- function() {

  #
  ## Test prob3
  #

  data("db.pars")

  rho <- db.pars[db.pars$param=="rho","value"]

  #
  ##
  #

  data("db.oscpars")
  row.names(db.oscpars) <- db.oscpars[,1]

  print(db.oscpars)

  sin2_12 <- db.oscpars[db.oscpars$param=="sin2_21","value"]
  sin2_13 <- db.oscpars[db.oscpars$param=="sin2_13","value"]
  sin2_23 <- db.oscpars[db.oscpars$param=="sin2_23","value"]
  squared <- TRUE

  dm2_21 <- db.oscpars[db.oscpars$param=="dm2_21","value"]
  dm2_23 <- db.oscpars[db.oscpars$param=="dm2_32","value"]

  MH <- db.oscpars[db.oscpars$param=="MH","value"]

  delta <- pi/2


  #
  ##
  #
  MH <- 1

  if(delta == 0) my.env$zero_cp <- TRUE

  is_anti <- T

#  numCores <- detectCores()
#  cl <- makeCluster(numCores)
#  registerDoParallel(cl)

  for(i in 1:1){

    setMNS(sin2_12, sin2_13, sin2_23, dm2_21, MH*dm2_23, delta, squared,is_anti)

    L <- 2500
    E <- 0.4182558

    nu_in  <- 2
    nu_out <- 1
    if(is_anti){
      nu_in <- -1*nu_in
      nu_out <- -1*nu_out
    }

    prob <- computeProb(nu_in,nu_out,L,rho,E,is_anti)

    MMH <- if(MH==1) factor("NH") else factor("IH")

    prob2 <- setAndComputeProb(sin2_12, sin2_13, sin2_23, dm2_21, dm2_23, MMH, delta, is_anti=T, L, rho, E, abs(nu_in), abs(nu_out))
    prob3 <- setAndComputeProb(sin2_12, sin2_13, sin2_23, dm2_21, dm2_23, MMH, delta, is_anti=T, L, rho, E, abs(nu_in), abs(nu_out))


    cat("Compare ",prob," ",prob2," ",prob3,"\n")
    #prob <- computeVacuumProb(1,1,L,E)
    #cat("prob: ",i," ",prob,"\n")
 }

  #return(0)

  cat("----Mix----\n")
  print(my.env$Mix)

  cat("----dm----\n")
  print(my.env$dm)

  cat("----dmVacVac----\n")
  print(my.env$dmVacVac)


  cat("----dmMatMat----\n")
  print(my.env$dmMatMat)

  cat("----dmMatVac----\n")
  print(my.env$dmMatVac)

  cat("----A----\n")
  print(my.env$A)

  cat("----Prob----\n")
  print(prob)

  cat("----Prob----\n")
  print(my.env$Prob)

  cat("sums: ",colSums(my.env$Prob),"\n")

  my.env$dump.params()

  invisible()

}# prob3



#' @export
#' @title oscillation probability
#' @description return oscillation probabilities for the transition nu_in -> nu_out
#' @details  probabilities should be pre-computed via \code{\link{computeProb}}
#' @param  flavour codes are: nu_  1:e 2:mu 3:tau  -1:e_bar -2:mu_bar -3:tau_bar
#' @examples  getProb(2,1): nu_mu -> nu_e \cr
#' getProb(-1,-1) nu_e_bar -> nu_e_bar
getProb <- function(nu_in, nu_out) {

  if( nu_in*nu_out <0 ) stop(call.=T, " asking for neutrino-antineutrino transition probability")

  nu_in <- abs(nu_in)
  nu_out <- abs(nu_out)

  invisible( my.env$Prob[nu_in,nu_out] )
}



#' @export
#' @title compute vacuum probability
#' @description  matrix of 3-flavour transition probabilities is filled, refer to \code{\link{computeProb}}
#' @param refers to \code{\link{computeProb}}
computeVacuumProb <- function(alpha, beta, L, E) {

  lovere <- 1.26693281*L/E

  dm21 <- my.env$dm21
  dm32 <- my.env$dm32
  dm31 <- my.env$dm31


  ss21 <- sin(dm21*lovere)^2
  ss32 <- sin(dm32*lovere)^2
  ss31 <- sin( (dm31)*lovere )^2
  #ss31 <- sin( (dm21+dm32)*lovere )^2

  mix <- my.env$Mix

  prob <- array(0,dim=c(3,3))

  re <- 1
  im <- 2

  for(ista in 1:3){
    for(iend in 1:3){

      prob[ista,iend]  = mix[ista,1,re]*mix[iend,1,re]*
        mix[ista,2,re]*mix[iend,2,re]*ss21;

      prob[ista,iend] = prob[ista,iend] + mix[ista,2,re]*mix[iend,2,re]*
        mix[ista,3,re]*mix[iend,3,re]*ss32;

      prob[ista,iend] = prob[ista,iend] + mix[ista,3,re]*mix[iend,3,re]*
        mix[ista,1,re]*mix[iend,1,re]*ss31;

      if ( iend == ista ) {
        prob[ista,iend]  = 1 - 4*prob[ista,iend];
      } else {
        prob[ista,iend]  =    -4*prob[ista,iend];
      }

    }#end

    prob[ista,3]=1.0-prob[ista,1]-prob[ista,2];

  }#start

  my.env$Prob <- prob

  if(missing(alpha) | missing(beta) ) return(invisible(NULL))

  nu_in <- abs(alpha)
  nu_out <- abs(beta)


  invisible( prob[nu_in,nu_out])

}#getVacuumProb


#' @export
#' @title set and compute ascill. probability nu_in -> nu_out
#' @description wrapper function to setup and compute transition probabilities, it calls in sequence \code{\link{setMNS}}) and \code{\link{computeProb}})
#' @details Assumes \code{sin} are passed as \code{sin^2} (i.e.: \code{squared=T} in \code{\link{setMNS}})\cr
#' \code{MH} is of type \code{factor}\cr
setAndComputeProb <- function(sin2_12, sin2_13, sin2_23, dm2_21, dm2_23, MH=factor(c("NH","IH")), delta, is_anti=c(T,F), L, rho, E, nu_in, nu_out)
 {

  #print(as.data.frame(as.list(environment())))

  iMH <- if(MH=="NH") 1 else -1;

  setMNS(x12=sin2_12, x13=sin2_13, x23=sin2_23, dm21=dm2_21, dmAtm=iMH*dm2_23, delta=delta, kSquared=TRUE, is_anti=is_anti)

  if(is_anti){
    nu_in  <- -1*nu_in
    nu_out <- -1*nu_out
  }

  P <- computeProb(nu_in, nu_out, L,rho,E, is_anti)

 }#setAndCompute



#' @export
#' @title compute oscillation probability
#' @description compute oscillations probabilities in the 3 flavour framework following the
#' linear propagation in constant density matter as derived in the paper or Barger at al.\cr
#' All transition probabilities are computed and can be retrieved via \code{\link{getProb}}.
#' the function itself returns the probability for the transition \code{alpha->beta}.\cr
#' This function must be called after a call to \code{\link{setMNS}} functions in order to
#' configure the oscillation parameters.
#' @param \code{alpha,beta} flavour code for the transition \code{alpha->beta}. Check \code{\link{getProb}}
#' for the index values. If not provided a \code{NULL} is returned.
computeProb <- function(alpha, beta, L, rho, E, is_anti=c(TRUE,FALSE)) {

#
##  minimally readapted from original C code
#

  useMassEigenstate <- FALSE

  #
  ## propagate Linear
  #


  if( (my.env$isAnti & !is_anti) | (!my.env$isAnti & is_anti) ){

    stop(call.=TRUE," Propagating neutrino flavor and MNS matrix definition differ ")
  }

  my.env$rho <- rho


  ## getTransitionMatrix

  antitype <- alpha
  phase_offset <- 0

  getM(E, rho/2, antitype ) #call this before getA!


  getA(L, E, rho/2, antitype, phase_offset)


  # now transition matrix A should be filled

  rawInputPsi <- array(0,dim=c(3,2))
  outputPsi   <- array(0,dim=c(3,2))

  probability <- array(0,dim=c(3,3))


  for(i in 1:3){

    for(j in 1:3){
     rawInputPsi[j,1] = 0.0; rawInputPsi[j,2] = 0.0;
    }

    if( !useMassEigenstate ){
      rawInputPsi[i,1] = 1.0;
    }else{
       stop(call.=T, " useMassEigenstate option not implemented yet")
    }

    outputPsi <- multiply_complex_matvec(rawInputPsi)

    #
    ###
    #

    probability[i,1] = probability[i,1] + outputPsi[1,1] * outputPsi[1,1] + outputPsi[1,2]*outputPsi[1,2];
    probability[i,2] = probability[i,2] + outputPsi[2,1] * outputPsi[2,1] + outputPsi[2,2]*outputPsi[2,2];
    probability[i,3] = probability[i,3] + outputPsi[3,1] * outputPsi[3,1] + outputPsi[3,2]*outputPsi[3,2];

  }#i

  #
  ## copy to global matrix
  #
  my.env$Prob <- probability

  #
  ##m Return if no transition specified (Probabiloty matrix can be accessed afterwards via global functions)
  #
  if( missing(alpha) | missing(beta) ) invisible(NULL)

  #
  ## Return the transition requested
  #
  nin  <- abs(alpha)
  nout <- abs(beta)

  prob <- probability[nin,nout]

  invisible(prob)

}#getProb



# @title multiply_complex_matvec
multiply_complex_matvec <- function( V ) {

  A <- my.env$A

  W <- array(0,dim=c(3,2))

  re <- 1
  im <- 2

  for(i in 1:3) {

    W[i,re] = A[i,1,re]*V[1,re]-A[i,1,im]*V[1,im]+
      A[i,2,re]*V[2,re]-A[i,2,im]*V[2,im]+
      A[i,3,re]*V[3,re]-A[i,3,im]*V[3,im] ;

     W[i,im] = A[i,1,re]*V[1,im]+A[i,1,im]*V[1,re]+
      A[i,2,re]*V[2,im]+A[i,2,im]*V[2,re]+
      A[i,3,re]*V[3,im]+A[i,3,im]*V[3,re] ;
  }


  invisible(W)
}


#' @export
#' @title set MNS matrix, mass splittings etc.
#' @description  determine the neutrino oscillation parameters
#' This function must be called before propagation functions
#' \code{\link{getProb}}, \code{\link{getVaccumProb}}
#' @param Specify the neutrino oscillation parameters, and energy,
#' the final boolean specifies which form  of mixing angle is input \cr
#'           \code{x12   ,  x13   ,  x23   ,  dm21  ,  dmAtm  ,  d_cp, kSquared , is_anti}\cr
#' Notation: m21=mSol \cr
#' If \code{dmAtm>0} Normal Hierachy is assumed (NH) and \code{dm32 = dmAtm (and dm31=dm21+dm32)}\cr
#' If \code{dmAtm<0} Inverted Hierarchy is assumed (IH), that is \code{dmAtm=dm31<0 and dm32 = -abs(dmAtm)-dm21 <0 }\cr
#' By default the "One dominant mass scale" approx is used, therefore only \code{dm21} and \code{dm32} are used\cr
#' \code{kSquared=T: sin^2(x) F: sin^2(2x)} \cr
#' \code{is_anti=T} if for anti-neutrino propagation \cr
#' @return Following global objects are filled: \code{mix} (MNS matrix), \code{dm,dmVacVac} (mass matrix).
#' @details Computation follows the derivation of Barger et al. \code{Phys. Rev. D22(1980) 2718.}, code is adapted form its original C-version
#' available here \code{www.phy.duke.edu/~raw22/public/Prob3++}
setMNS <- function(x12, x13, x23, dm21, dmAtm, delta, kSquared=c(TRUE,FALSE), is_anti=c(TRUE,FALSE) ) {

  ldm32 <- dmAtm

  dm32 <- dmAtm
  dm31 <- dm21+dm32

  if( dmAtm < 0 ) {

     my.env$MH <- -1

     ldm32 <- dmAtm - dm21

     dm31 <- dmAtm
     dm32 <- dm31-dm21

  }


  sin21 <- 0
  sin13 <- 0
  sin23 <- 0

  if( kSquared ){

    sin12 <- sqrt(x12)
    sin13 <- sqrt(x13)
    sin23 <- sqrt(x23)

  }else{

    sin12 <- sqrt(0.5*(1 - sqrt(1 - x12)))
    sin13 <- sqrt(0.5*(1 - sqrt(1 - x13)))
    sin23 <- sqrt(0.5*(1 - sqrt(1 - x23)))

  }


  my.env$dmAtm <- dmAtm
  my.env$dmSol <- dm21

  my.env$dm21 <- dm21
  my.env$dm32 <- dm32
  my.env$dm31 <- dm31

  my.env$isAnti <- is_anti

  if(is_anti) delta <- -1*delta

  init_mixing_matrix(dm21, ldm32, sin12, sin23, sin13, delta)

  invisible()

}



# init MNS mixing matrix and mass matrix in the vacuum
init_mixing_matrix <- function(dm21, dm23, s12, s23, s13, dcp) {


  #
  #
  ###
  #
  #

  tol <- 5e-9

  mVac <- array(0,dim=3)

  mVac[1] = 0
  mVac[2] = dm21
  mVac[3] = dm21 + dm23

  # remove any degeneracy
  if (dm21==0.) mVac[1] = mVac[1] - tol;
  if (dm23==0.) mVac[3] = mVac[3] + tol;

  dmVacVac <- array(0,dim=c(3,3))

  dmVacVac[1,1] = 0
  dmVacVac[2,2] = 0
  dmVacVac[3,3] = 0

  dmVacVac[1,2] = mVac[1]-mVac[2]
  dmVacVac[2,1] = -dmVacVac[1,2]

  dmVacVac[1,3] = mVac[1]-mVac[3]
  dmVacVac[3,1] = -dmVacVac[1,3]

  dmVacVac[2,3] = mVac[2]-mVac[3]
  dmVacVac[3,2] = -dmVacVac[2,3]


  #
  #
  ###
  #
  #

  if ( s12>1.0 ) s12=1.0;
  if ( s23>1.0 ) s23=1.0;
  if ( s13>1.0 ) s13=1.0;

  sd = sin( dcp );
  cd = cos( dcp );

  c12 = sqrt(1.0-s12*s12);
  c23 = sqrt(1.0-s23*s23);
  c13 = sqrt(1.0-s13*s13);

  re <- 1
  im <- 2

  Mix <- array(0,dim=c(3,3,2))

  if ( my.env$matrixType == matrixTypes[1] )
  {
    #cat(" standard type ",dcp,cd,sd,s13,"\n",sep=" ")
    Mix[1,1,re] =  c12*c13;
    Mix[1,1,im] =  0.0;
    Mix[1,2,re] =  s12*c13;
    Mix[1,2,im] =  0.0;
    Mix[1,3,re] =  s13*cd;
    Mix[1,3,im] = -s13*sd;
    Mix[2,1,re] = -s12*c23-c12*s23*s13*cd;
    Mix[2,1,im] =         -c12*s23*s13*sd;
    Mix[2,2,re] =  c12*c23-s12*s23*s13*cd;
    Mix[2,2,im] =         -s12*s23*s13*sd;
    Mix[2,3,re] =  s23*c13;
    Mix[2,3,im] =  0.0;
    Mix[3,1,re] =  s12*s23-c12*c23*s13*cd;
    Mix[3,1,im] =         -c12*c23*s13*sd;
    Mix[3,2,re] = -c12*s23-s12*c23*s13*cd;
    Mix[3,2,im] =         -s12*c23*s13*sd;
    Mix[3,3,re] =  c23*c13;
    Mix[3,3,im] =  0.0;
  }
  else
  {
    Mix[1,1,re] =  c12;
    Mix[1,1,im] =  0.0;
    Mix[1,2,re] =  s12*c23;
    Mix[1,2,im] =  0.0;
    Mix[1,3,re] =  s12*s23;
    Mix[1,3,im] =  0.0;
    Mix[2,1,re] = -s12*c13;
    Mix[2,1,im] =  0.0;
    Mix[2,2,re] =  c12*c13*c23+s13*s23*cd;
    Mix[2,2,im] =              s13*s23*sd;
    Mix[2,3,re] =  c12*c13*s23-s13*c23*cd;
    Mix[2,3,im] =             -s13*c23*sd;
    Mix[3,1,re] = -s12*s13;
    Mix[3,1,im] =  0.0;
    Mix[3,2,re] =  c12*s13*c23-c13*s23*cd;
    Mix[3,2,im] =             -c13*s23*sd;
    Mix[3,3,re] =  c12*s13*s23+c13*c23*cd;
    Mix[3,3,im] =              c13*c23*sd;
  }


  #
  ## save to globals
  #

  my.env$dm21 <- dm21
  my.env$dm32 <- dm23

  my.env$s12   <- s12
  my.env$s23   <- s23
  my.env$s31   <- s13
  my.env$delta <- dcp

  my.env$dmVacVac <- dmVacVac
  my.env$dm <- dmVacVac
  my.env$Mix <- Mix

  invisible()
}


#
#
#
#
#
#
#
# getM
getM <- function( Enu, rho, antitype ) {

#
# minimally readapted from the C version
#
  alpha = beta = gamma = fac = arg = tmp = 0;
  alphaV = betaV = gammaV = argV = tmpV = 0;
  theta0 = theta1 = theta2 = 0;
  theta0V = theta1V = theta2V = 0;

  mMatU <- array(0,dim=3)
  mMatV <- array(0,dim=3)
  mMat  <- array(0,dim=3)

  tworttwoGf = 1.52588e-4

  if (antitype<0) fac = tworttwoGf*Enu*rho # Anti-neutrinos
  else            fac = -tworttwoGf*Enu*rho # Real-neutrinos


  #
  # get some global vars
  #
  elec <- my.env$elec
  muon <- my.env$muon
  tau  <- my.env$tau

  dmVacVac <- my.env$dmVacVac
  Mix <- my.env$Mix

  #
  #

  alpha  = fac + dmVacVac[1,2] + dmVacVac[1,3];
  alphaV = dmVacVac[1,2] + dmVacVac[1,3];


  re <- 1
  im <- 2


  if( my.env$zero_cp ) {

    beta = dmVacVac[1,2]*dmVacVac[1,3] +
      fac*(dmVacVac[1,2]*(1.0 -
                             Mix[elec,2,re]*Mix[elec,2,re] -
                             Mix[elec,2,im]*Mix[elec,2,im]) +
             dmVacVac[1,3]*(1.0-
                               Mix[elec,3,re]*Mix[elec,3,re] -
                               Mix[elec,3,im]*Mix[elec,3,im]));

    betaV = dmVacVac[1,2]*dmVacVac[1,3];

  } else {

    beta = dmVacVac[1,2]*dmVacVac[1,3] +
      fac*(dmVacVac[1,2]*(1.0 -
                             Mix[elec,2,re]*Mix[elec,2,re]) +
             dmVacVac[1,3]*(1.0-
                               Mix[elec,3,re]*Mix[elec,3,re]));
    betaV = dmVacVac[1,2]*dmVacVac[1,3];

  }#


  if( my.env$zero_cp) {

    gamma = fac*dmVacVac[1,2]*dmVacVac[1,3]*
      (Mix[elec,1,re]*Mix[elec,1,re]+Mix[elec,1,im]*Mix[elec,1,im]);
    gammaV = 0.0;

  } else {

    gamma = fac*dmVacVac[1,2]*dmVacVac[1,3]*
      (Mix[elec,1,re]*Mix[elec,1,re]);
    gammaV = 0.0;

  }#



  #
  ##
  #


  # Compute the argument of the arc-cosine #
  tmp = alpha*alpha-3.0*beta;
  tmpV = alphaV*alphaV-3.0*betaV;

  if (tmp<0.0) {
    cat(" getM: alpha^2-3*beta < 0 !\n")
    tmp = 0.0;
  }

  #
  ## Equation (21)
  #

  arg = (2.0*alpha*alpha*alpha-9.0*alpha*beta+27.0*gamma)/
    (2.0*sqrt(tmp*tmp*tmp));
  if (abs(arg)>1.0) arg = arg/abs(arg);


  argV = (2.0*alphaV*alphaV*alphaV-9.0*alphaV*betaV+27.0*gammaV)/
    (2.0*sqrt(tmpV*tmpV*tmpV));
  if (abs(argV)>1.0) argV = argV/abs(argV);

#  cat("zero_cp ",my.env$zero_cp,"\n")
#  cat("beta ",beta," ",dmVacVac[1,2]," ",dmVacVac[1,3],"\n")
#  cat("betaV ",betaV,"\n")
#  cat("gamma ",gamma,"\n")
#  cat("gammaV ",gammaV,"\n")

#  cat("arg ",arg,"\n")
#  cat("argV ",argV,"\n")
#  cat("tmp ",tmp,"\n")
#  cat("tmpV ",tmpV,"\n")



  # These are the three roots the paper refers to #
  theta0 = acos(arg)/3.0;
  theta1 = theta0-(2.0*pi/3.0);
  theta2 = theta0+(2.0*pi/3.0);

#  cat("theta",theta0,theta1,theta2,"\n",sep=" ")

  theta0V = acos(argV)/3.0;
  theta1V = theta0V-(2.0*pi/3.0);
  theta2V = theta0V+(2.0*pi/3.0);


  mMatU[1] = mMatU[2] = mMatU[3] = -(2.0/3.0)*sqrt(tmp);
  mMatU[1] = mMatU[1]*cos(theta0);
  mMatU[2] = mMatU[2]*cos(theta1);
  mMatU[3] = mMatU[3]*cos(theta2);

  tmp = dmVacVac[1,1] - alpha/3.0;
  mMatU[1] = mMatU[1]+tmp; mMatU[2] = mMatU[2]+tmp; mMatU[3] = mMatU[3]+tmp;

  mMatV[1] = mMatV[2] = mMatV[3] = -(2.0/3.0)*sqrt(tmpV);
  mMatV[1] = mMatV[1]*cos(theta0V);
  mMatV[2] = mMatV[2]*cos(theta1V);
  mMatV[3] = mMatV[3]*cos(theta2V);

  tmpV = dmVacVac[1,1] - alphaV/3.0;
  mMatV[1] = mMatV[1]+tmpV; mMatV[2] = mMatV[2]+tmpV; mMatV[3] = mMatV[3]+tmpV;


  # Sort according to which reproduce the vaccum eigenstates
  for(i in 1:3) {

    tmpV = abs(dmVacVac[i,1]-mMatV[1])

    k = 1
    for(j in 2:3) {

      tmp =   abs(dmVacVac[i,1] - mMatV[j])

      if(tmp<tmpV) {
        k=j; tmpV=tmp;
      }

    }#j

    mMat[i] = mMatU[k]

  }#i


  for (i in 1:3) {
    for (j in 1:3) {
      my.env$dmMatMat[i,j] = mMat[i] - mMat[j];
      my.env$dmMatVac[i,j] = mMat[i] - dmVacVac[j,1];
    }
  }

#  cat("---mMatU--- \n")
#  print(mMatU)

#  cat("---mMatV--- \n")
#  print(mMatV)

#  cat("---mMat -- \n")
#  print(mMat)


#  cat("---dmMatMat -- \n")
#  print(my.env$dmMatMat)

#  cat("---dmMatVac -- \n")
#  print(my.env$dmMatVac)


#  cat("getM ",Enu," ",rho," fac ", fac,"\n")

# stopifnot(F)

  invisible()

}# getM


#
#
#
#
# @title getA
getA <- function(L, E, rho, antitype, phase_offset) {

#
# minimally readapted from the C version
#
  arc = c = s = 0

  X <- array(0, dim=c(3,3,2))
  product <- array(0, dim=c(3,3,3,2))

  # (1/2)*(1/(h_bar*c)) in units of GeV/(eV^2-km)
  LoEfac <- 2.534

  if( phase_offset == 0 ){
    product <- get_product(L,E,rho,antitype)
  }



  #
  re <- 1
  im <- 2
  #

  # Make the sum with the exponential factor
  for( k in 1:3 ) {

    arg = -LoEfac * my.env$dmMatVac[k][1]*L/E

    if(k==2) arg = arg + phase_offset

    c = cos(arg)
    s = sin(arg)

    for( i in 1:3 ) {
      for( j in 1:3) {

        if(my.env$zero_cp){

          X[i,j,re] = X[i,j,re] + c*product[i,j,k,re] - s*product[i,j,k,im];
          X[i,j,im] = X[i,j,im] + c*product[i,j,k,im] + s*product[i,j,k,re];

        }else{

          X[i,j,re] = X[i,j,re] + c*product[i,j,k,re];
          X[i,j,im] = X[i,j,im] + s*product[i,j,k,re];

        }#zero_cp

      }#j
    }#i

  }#k


  # Compute the product with the mixing matrices

  A <- array(0,dim=c(3,3,2))

  Mix <- my.env$Mix

  for(n in 1:3 ){
    for(m in 1:3 ){
      for(i in 1:3 ){
        for(j in 1:3 ){

          if(my.env$zero_cp){

            A[n,m,re] = A[n,m,re] +
              Mix[n,i,re]*X[i,j,re]*Mix[m,j,re] +
              Mix[n,i,re]*X[i,j,im]*Mix[m,j,im] +
              Mix[n,i,im]*X[i,j,re]*Mix[m,j,im] -
              Mix[n,i,im]*X[i,j,im]*Mix[m,j,re];

            A[n,m,im] = A[n,m,im] +
              Mix[n,i,im]*X[i,j,im]*Mix[m,j,im] +
              Mix[n,i,im]*X[i,j,re]*Mix[m,j,re] +
              Mix[n,i,re]*X[i,j,im]*Mix[m,j,re] -
              Mix[n,i,re]*X[i,j,re]*Mix[m,j,im];

          }else{

            A[n,m,re] = A[n,m,re]+
              Mix[n,i,re]*X[i,j,re]*Mix[m,j,re];

            A[n,m,im] = A[n,m,im]+
              Mix[n,i,re]*X[i,j,im]*Mix[m,j,re];

          }#zero_cp

        }#n
      }#m
    }#i
  }#j

  my.env$A <- A

  invisible()

}# getA



#
#
#
#
# @get_product
get_product <- function(L, E, rho, antitype) {

  #
  # minimally readapted from the C version
  #

#  cat("get_product\n")


  twoEHmM <- array(0,dim=c(3,3,3,2))

  # (1/2)*(1/(h_bar*c)) in units of GeV/(eV^2-km)
  # Reverse the sign of the potential depending on neutrino type

  tworttwoGf = 1.52588e-4

  fac = 0

  # If we're doing matter effects for electron neutrinos
  if (antitype<0) fac =  tworttwoGf*E*rho  # Anti-neutrinos
  else            fac = -tworttwoGf*E*rho      # Real-neutrinos

#  cat("fac= ",fac,"\n")

  #
  re <- 1
  im <- 2
  #

  Mix <- my.env$Mix

  for(n in 1:3) {
    for( m in 1:3 ) {


     if(my.env$zero_cp){

       #cat("zero_cp ",my.env$zero_cp,"\n")

       twoEHmM[n,m,1,re] =
         -fac*(Mix[1,n,re]*Mix[1,m,re]+Mix[1,n,im]*Mix[1,m,im]);

       twoEHmM[n,m,1,im] =
         -fac*(Mix[1,n,re]*Mix[1,m,im]-Mix[1,n,im]*Mix[1,m,re]);

       twoEHmM[n,m,2,re] = twoEHmM[n,m,3,re] = twoEHmM[n,m,1,re];

       twoEHmM[n,m,2,im] = twoEHmM[n,m,3,im] = twoEHmM[n,m,1,im];

     }else{

       twoEHmM[n,m,1,re] =
         -fac*(Mix[1,n,re]*Mix[1,m,re]);

       twoEHmM[n,m,1,im] = 0;
       twoEHmM[n,m,2,re] = twoEHmM[n,m,3,re] = twoEHmM[n,m,1,re];
       twoEHmM[n,m,2,im] = twoEHmM[n,m,3,im] = twoEHmM[n,m,1,im];
     }#zero_cp

      if(n==m) {
        for(j in 1:3) {twoEHmM[n,m,j,re] = twoEHmM[n,m,j,re] - my.env$dmMatVac[j,n]}
      }

    }#m
  }#n

#  cat("twoEHnM\n")
#  print(twoEHmM)

  #
  ##
  ###
  ##
  #

  product <- array(0,dim=c(3,3,3,2))

  dmMatMat <- my.env$dmMatMat


  for(i in 1:3) {
    for(j in 1:3) {
      for(k in 1:3){

        if( my.env$zero_cp) {

          #cat("zero_cp ",my.env$zero_cp,"\n")

          product[i,j,1,re] = product[i,j,1,re] +
            twoEHmM[i,k,2,re]*twoEHmM[k,j,3,re] -
            twoEHmM[i,k,2,im]*twoEHmM[k,j,3,im];

          product[i,j,1,im] = product[i,j,1,im] +
            twoEHmM[i,k,2,re]*twoEHmM[k,j,3,im] +
            twoEHmM[i,k,2,im]*twoEHmM[k,j,3,re];

          product[i,j,2,re] = product[i,j,2,re] +
            twoEHmM[i,k,3,re]*twoEHmM[k,j,1,re] -
            twoEHmM[i,k,3,im]*twoEHmM[k,j,1,im];

          product[i,j,2,im] = product[i,j,2,im] +
            twoEHmM[i,k,3,re]*twoEHmM[k,j,1,im] +
            twoEHmM[i,k,3,im]*twoEHmM[k,j,1,re];

          product[i,j,3,re] = product[i,j,3,re] +
            twoEHmM[i,k,1,re]*twoEHmM[k,j,2,re] -
            twoEHmM[i,k,1,im]*twoEHmM[k,j,2,im];

          product[i,j,3,im] = product[i,j,3,im] +
            twoEHmM[i,k,1,re]*twoEHmM[k,j,2,im] +
            twoEHmM[i,k,1,im]*twoEHmM[k,j,2,re];

          #cat(product[i,j,1,re],product[i,j,1,im],product[i,j,2,re],product[i,j,2,re],product[i,j,3,re],product[i,j,3,re],"\n",sep=" ")

        }else{

          product[i,j,1,re] = product[i,j,1,re] +
            twoEHmM[i,k,2,re]*twoEHmM[k,j,3,re];

          product[i,j,2,re] =  product[i,j,2,re] +
            twoEHmM[i,k,3,re]*twoEHmM[k,j,1,re];

          product[i,j,3,re] = product[i,j,3,re] +
            twoEHmM[i,k,1,re]*twoEHmM[k,j,2,re];
        }#zero_cp
      }#k

      if( my.env$zero_cp) {

        product[i,j,1,re] = product[i,j,1,re]/(dmMatMat[1,2]*dmMatMat[1,3]);
        product[i,j,1,im] = product[i,j,1,im]/(dmMatMat[1,2]*dmMatMat[1,3]);
        product[i,j,2,re] = product[i,j,2,re]/(dmMatMat[2,3]*dmMatMat[2,1]);
        product[i,j,2,im] = product[i,j,2,im]/(dmMatMat[2,3]*dmMatMat[2,1]);
        product[i,j,3,re] = product[i,j,3,re]/(dmMatMat[3,1]*dmMatMat[3,2]);
        product[i,j,3,im] = product[i,j,3,im]/(dmMatMat[3,1]*dmMatMat[3,2]);


        #cat(product[i,j,1,re],product[i,j,1,im],product[i,j,2,re],product[i,j,2,re],product[i,j,3,re],product[i,j,3,re],"\n",sep=" ")

      }else{

        product[i,j,1,re] = product[i,j,1,re]/(dmMatMat[1,2]*dmMatMat[1,3]);
        product[i,j,2,re] = product[i,j,2,re]/(dmMatMat[2,3]*dmMatMat[2,1]);
        product[i,j,3,re] = product[i,j,3,re]/(dmMatMat[3,1]*dmMatMat[3,2]);
      }#zero_cp

    }#j
  }#i


  invisible(product)

}# get_product






