######################################
#
#
#    Studyes with prob3-R-version functions
#    3-flavour neutrino oscillation probability in
#    constant density matter
#
#
#####################################


library(plyr)
library(compiler)
library(parallel)
library(doParallel)
library(reshape2)
library(fields)
library(ggplot2)
library(neutRino)


#
##
#

enableJIT(3)

#
##
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

delta <- db.oscpars[db.oscpars$param=="dcp","value"]

# flavour indexes
nu_in   <- my.env$muon
nu_out  <- my.env$elec

#
## Vectorize main function for computation prob.
#

vec.computeProb <- Vectorize(computeProb)

vec.setAndComputeProb <- Vectorize(setAndComputeProb)

#
##
#

base.pars <- list(sin2_12=sin2_12, sin2_13=sin2_13, sin2_23=sin2_23, dm2_21=dm2_21, dm2_23=dm2_23, rho=rho)

#
##
#
print("--- Scan grid --- ")

scan.MH    <- factor(c("NH","IH"))
scan.delta <- c(0,pi/4,pi/2,pi)
scan.anti  <- c(F,T)



scan.grid <- grid.list(MH=scan.MH, delta = scan.delta, is_anti=scan.anti)

print( list.to.df(scan.grid) )

##
#params.grid <- grid.list(MH=scan.MH, delta = scan.delta, is_anti=scan.anti, pre=base.pars, post = list(L=2500, nu_in=nu_mu, nu_out=nu_e) )
#print( list.to.df(params.grid) )

L <- 2300
Emin <- 0.1 #GeV
Emax <- 8

rho <- 2.8
sin2_23 <- 0.49
sin2_12 <- 0.32
sin2_13 <- 0.026
dm2_21 <- 7.62e-5
dm2_23 <- 2.53e-3

E <- logspaced(1000,Emin,Emax)


do.parallel <- F
if( do.parallel ) registerDoParallel( makeCluster(detectCores()) )


scan.list <- list(MH=scan.MH, delta = scan.delta, is_anti=scan.anti)


g <- fun.grid.scan(  scan.list, fun = vec.setAndComputeProb
                    ,sin2_12=sin2_12, sin2_13=sin2_13, sin2_23=sin2_23, dm2_21=dm2_21, dm2_23=dm2_23, rho=rho
                    ,L=L
                    ,nu_in=nu_in, nu_out=nu_out
                    ,E=E
                    ,.parallel=do.parallel )



g <- transform(g,delta=factor(delta,levels=c(0,pi/4,pi/2,pi),labels=c("0","pi/4","pi/2","pi")))

g <- transform(g,is_anti=factor(is_anti,levels=c(FALSE,TRUE),labels=c("neutrino","anti-neutrino")))

names(g)[names(g)=="is_anti"] <- "type"

head(g)
tail(g)


transition <- nu.trans.exp(nu_in,nu_out)

title <- substitute( paste("Oscillation Probability ",transition," (L=",L," km)"), list(L=L,transition=transition) )


p <- ggplot(g, aes(x=E,y=value, color=delta) ) + geom_line(size=0.5)
p <- p + facet_grid(MH~type,labeller=labeller(type=label_value,MH=label_both), as.table=F)
p <- p + ylim(0,0.15)

p <- p + theme_bw()
p <- p + theme(panel.grid.major = element_line(colour = "black", linetype="dotted"))
p <- p + xlab(expression(E[nu]~GeV)) + ylab(expression(Probability))

p <- p + labs(title=title)
p <- p + theme(axis.title = element_text(size = 22, colour = "black", face=NULL,angle = 0))
p <- p + theme(axis.text = element_text(size = 12, colour = "black", face=NULL,angle = 0))

p <- p + theme(title = element_text(colour="black", size=22, face="bold"))
p <- p + theme(strip.text.x = element_text(size = 20, colour = "black", face=NULL,angle = 0))
p <- p + theme(strip.text.y = element_text(size = 20, colour = "black", face=NULL,angle = 0))
p <- p + theme(strip.background = element_rect(colour="red", fill="orange"))

p <- p + theme(legend.text  = element_text(colour="black", size = 26, face = NULL),legend.text.align=0)
p <- p + theme(legend.title = element_text(colour="black", size = 26, face = NULL),legend.title.align=0.4)
p <- p + theme(legend.key=element_rect(fill="white",linetype=0), legend.position="right")
p <- p + scale_color_discrete(name=bquote(delta[CP]),labels=c(expression(0),expression(pi/4),expression(pi/2),expression(pi)))

p


