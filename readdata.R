################################

library(neutRino)


library(tikzDevice)
library(data.table)
library(ggplot2)
library(weights)
library(fields)
library(stats)
library(plyr)
library(HistogramTools)


###############################


#
## Load package-data
#

#default oscillations parameters
data("db.oscpars")
row.names(db.oscpars) <- db.oscpars[,1] #first col also as raw names

#default oscillation parameters for sterile models
data("db.oscpars.sterile")
row.names(db.oscpars.sterile) <- db.oscpars.sterile[,1]

#generic utility parameters (particle masses, earth density,)
data("db.pars")
row.names(db.pars) <- db.oscpars[,1]

# Data set of reconstructed neutrino events
data("data_nue_selected")


#
##
#


nevents <- nrow(data_nue_selected)

# add osc prob of the event (init at 1)
nue_data <- cbind(data_nue_selected, rep(1, nevents))



# Names of event variables (features)
var.names <- c("event", "iflav", "isel", "iint", "ifv", "Ereco", "Etrue", "pathLen", "weight", "IsNuE", "IsBkg", "IsGamma", "OscProb")

nvars <- length(var.names)

names(nue_data) <- var.names



#
## Add Variable
#


sample_levels <- factor(c("NuE","Bkg","Gamma"))
inter_levels  <- factor(c("CCQE","CCRES","CCDIS","CCCOH","NC1Pi0","NCOth"))
flav_levels   <- factor(c("NuE","ANuE","NuMu","ANuMu"))
fv_levels     <- factor(c("Inside","Out Gamma","Out Oth."))
selec_levels  <- factor(c("CCQE","CCnQE","Gamma")) #,"CCInc"))


nue_data$FV           <- fv_levels[nue_data$ifv+1]
nue_data$interaction  <- inter_levels[nue_data$iint+1]
nue_data$flavor       <- flav_levels[nue_data$iflav+1]
nue_data$selection    <- selec_levels[nue_data$isel+1]


c_sample <- apply(nue_data,1,function(x){
                            x<-as.list(x)
                            type <- "NA"
                            if(x$IsNuE==1)        type <- sample_levels[1]
                            else if(x$IsBkg==1)   type <- sample_levels[2]
                            else if(x$IsGamma==1) type <- sample_levels[3]
                            return(as.factor(type))
                            })

nue_data$sample <- c_sample




# convert bool to logical for convenience

nue_data$IsNuE <- as.logical(nue_data$IsNuE)
nue_data$IsBkg <- as.logical(nue_data$IsBkg)
nue_data$IsGamma <- as.logical(nue_data$IsGamma)


#
## divide in subsamples
#

data_nue   <- nue_data[nue_data$IsNuE, ]
data_bkg   <- nue_data[nue_data$IsBkg, ]
data_gamma <- nue_data[nue_data$IsGamma, ]



## Ereco energy binning units: MeV

bins.ereco <- c(200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200,
    2350, 2500, 2700, 3000, 3300, 3500, 4000, 4400, 5000, 6000, 10000)

bins.etrue <- seq(from = 0, to = 20000, length.out = 201)


bins.egamma <- c(200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1700, 1900, 2200, 2500, 2800, 4000,
    10000)


##

label.ereco <- "Ereco [MeV]"
label.etrue <- "Etrue [MeV]"

label.entries <- "Entries"

#
## create base histograms
#


# 1D

hNuE_Etrue <- histo1d(x=data_nue$Etrue,breaks=bins.etrue ,weight=data_nue$weight
                      ,main="NuE - True Energy", xlab=label.etrue, ylab=label.entries,show=F)

hNuE_Ereco <- histo1d(x=data_nue$Ereco,breaks=bins.ereco ,weight=data_nue$weight
                      ,main="NuE - Reconstructed Energy", xlab=label.ereco, ylab=label.entries,show=F)


hBkg_Etrue <- histo1d(x=data_bkg$Etrue,breaks=bins.etrue ,weight=data_bkg$weight
                      ,main="Bkg - True Energy", xlab=label.etrue, ylab=label.entries,show=F)

hBkg_Ereco <- histo1d(x=data_bkg$Ereco,breaks=bins.etrue ,weight=data_bkg$weight
                      ,main="Bkg - Reconstructed Energy", xlab=label.ereco, ylab=label.entries,show=F)


hGamma_Etrue <- histo1d(x=data_gamma$Etrue,breaks=bins.egamma ,weight=data_gamma$weight
                      ,main="Gamma - True Energy", xlab=label.etrue, ylab=label.entries,show=F)

hGamma_Ereco <- histo1d(x=data_gamma$Ereco,breaks=bins.ereco ,weight=data_gamma$weight
                      ,main="Gamma - Reconstructed Energy", xlab=label.ereco, ylab=label.entries,show=F)


# 2D

hNuE_EtrueVsEreco <- histo2d(x=data_nue$Etrue, y=data_nue$Ereco, breaks.x=bins.etrue, breaks.y=bins.ereco, weight=data_nue$weight
                      ,main="NuE - True vs Reco Energy", xlab=label.etrue, ylab=label.ereco,show=F)


hBkg_EtrueVsEreco <- histo2d(x=data_bkg$Etrue, y=data_bkg$Ereco, breaks.x=bins.etrue, breaks.y=bins.ereco, weight=data_bkg$weight
                             ,main="NuE - True vs Reco Energy", xlab=label.etrue, ylab=label.ereco,show=F)


hGamma_EtrueVsEreco <- histo2d(x=data_gamma$Etrue, y=data_gamma$Ereco, breaks.x=bins.egamma, breaks.y=bins.ereco, weight=data_gamma$weight
                             ,main="NuE - True vs Reco Energy", xlab=label.etrue, ylab=label.ereco,show=F)


#
##
#

ggplot(nue_data,aes(x=sample,fill=interaction)) + geom_bar(position="dodge") + labs(title="Composition") + theme_bw()

ggplot(nue_data,aes(x=sample,fill=FV)) + geom_bar(position="dodge") + labs(title="Composition") + theme_bw() + facet_grid(selection~FV)

ggplot(nue_data,aes(x=sample,fill=selection)) + geom_bar(position="dodge") + labs(title="Composition") + theme_bw()

ggplot(data_nue,aes(x=Etrue,y=Ereco)) + geom_bin2d(breaks=list(x=bins.etrue,y=bins.ereco)) + labs(title="Composition") + theme_bw() + facet_grid(selection~FV)

plot(hGamma_EtrueVsEreco,legend=F)


