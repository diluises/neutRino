################################

library(neutRino)

library(inline)
library(Rcpp)
library(RcppArmadillo)

library(tikzDevice)
library(data.table)
library(ggplot2)
library(weights)
library(stats)
library(plyr)



###############################



Nu.ErecoAsEle <- function(Elepton, Plepton, cos.theta) {

    Mp2 <- 1
    Mn2 <- 1
    me2 <- 1
    Eb <- 27  #MeV

    Ereco <- (Mp2 - (Mn - Eb)^2 - me2 + 2(Mn - Eb) * Elepton)
    Ereco <- Ereco/2/(Mn - Eb - Elepton + Plepton * cos.theta)

    return(Ereco)
}


Nu.ErecoAsMuon <- function(Elepton, Plepton, cos.theta) {

    Mp2 <- 1
    Mn2 <- 1
    ml2 <- 1
    Eb <- 27  # MeV

    Ereco <- (Mp2 - (Mn - Eb)^2 - ml2 + 2(Mn - Eb) * Elepton)
    Ereco <- Ereco/2/(Mn - Eb - Elepton + Plepton * cos.theta)

    return(Ereco)
}



LogLike.Histo <- function(hobs, hexp) {

    LL <- hobs - hexp - hobs * log(hobs/hexp)

    LL[!is.finite(LL)] <- 0

    LL <- sum(LL, na.rm = TRUE)

    return(-2 * LL)

}


Chi2.Histo <- function(hobs, herr, hexp) {

    chi2 <- 0.5 * sum((hobs - hexp)^2/(herr^2))

    return(chi2)
}






####

set.seed(1234)

####

Stats <- function(x) {

    Mean <- mean(x, na.rm = TRUE)
    Median <- median(x, na.rm = TRUE)
    SD <- sd(x, na.rm = TRUE)
    Min <- min(x, na.rm = TRUE)
    Max <- max(x, na.rm = TRUE)

    xx <- x[x > Median]
    Qr <- quantile(xx, 0.68)
    xx <- x[x <= Median]
    Ql <- quantile(xx, 1 - 0.68)

    return(c(Mean = Mean, SD = SD, Min = Min, Max = Max, Median = Median, Ql = Ql, Qr = Qr))

}

stat_fluct <- function(n, opt = "Poisson") {

    if (n < 0) {
        return(0)
    }
    dn <- 0
    if (opt == "Gaus")
        dn <- rnorm(1, mean = 0, sd = sqrt(n)) else dn <- rpois(1, lambda = sqrt(n))

    if (n + dn > 0)
        return(n + dn) else return(0)

}



DataRead.AsFrame <- function(file = file.name, check.firstCol = T, ...) {

    data <- read.delim(file = file.name, ...)

    if (check.firstCol) {
        # remove first colums if all NA's

        nevents <- nrow(data)

        col.1 <- data[, 1]

        nna <- sum(is.na(col.1))

        if (nna == nevents) {
            data <- data[, 2:ncol(data)]
        }
    }

    return(data)
}



####



file.name <- "data/nue_events_list.txt"

nue_data <- DataRead.AsFrame(file = file.name, check.firstCol = T, sep = " ", header = F)

nevents <- nrow(nue_data)

# add event osc prob (init at 1)
nue_data <- cbind(nue_data, rep(1, nevents))



# Names of event variables
var.names <- c("event", "iflav", "isel", "iint", "ifv", "Ereco", "Etrue", "pathLen", "weight", "IsNuE", "IsBkg", "IsGamma", "OscProb")

nvars <- length(var.names)

names(nue_data) <- var.names

## Ereco energy binning units: MeV

bins.ereco <- c(200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200,
    2350, 2500, 2700, 3000, 3300, 3500, 4000, 4400, 5000, 6000, 10000)

bins.etrue <- seq(from = 0, to = 20000, length.out = 200)


bins.egamma <- c(200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1700, 1900, 2200, 2500, 2800, 4000,
    10000)


##

labl.ereco <- "Ereco [MeV]"
labl.etrue <- "Etrue [MeV]"

labl.entries <- "Entries"


hErecoNuE <- wtd.hist(main = "NuE Selection", xlab = labl.ereco, ylab = labl.entries, x = nue_data[nue_data$IsNuE == 1, ]$Ereco,
    weight = nue_data[nue_data$IsNuE == 1, ]$weight, freq = T, breaks = bins.ereco)

hErecoBkg <- wtd.hist(main = "Background", xlab = labl.ereco, ylab = labl.entries, x = nue_data[nue_data$IsBkg == 1, ]$Ereco, weight = nue_data[nue_data$IsBkg ==
    1, ]$weight, freq = T, breaks = bins.ereco)

hErecoGamma <- wtd.hist(main = "Gamma Selection", xlab = labl.ereco, ylab = labl.entries, x = nue_data[nue_data$IsGamma == 1, ]$Ereco,
    weight = nue_data[nue_data$IsGamma == 1, ]$weight, freq = T, breaks = bins.ereco)


hEtrueNuE <- wtd.hist(main = "NuE Selection", xlab = labl.etrue, ylab = labl.entries, x = nue_data[nue_data$IsNuE == 1, ]$Etrue,
    weight = nue_data[nue_data$IsNuE == 1, ]$weight, freq = T, breaks = bins.etrue)



# apply stat fluc around the bin content value get the transpost in order ot have a column for each bin
opt <- "Gaus"
hist.fluct <- t(replicate(10000, sapply(hErecoNuE$counts, stat_fluct, opt = opt)))

hist.fluct.summary <- apply(hist.fluct, 2, Stats)

hist.fluct.Median <- apply(hist.fluct.summary, 2, function(x) as.list(x)$Median)
hist.fluct.Ql <- apply(hist.fluct.summary, 2, function(x) as.list(x)$Ql)
hist.fluct.Qr <- apply(hist.fluct.summary, 2, function(x) as.list(x)$Qr)

hist.fluct.summary

plot(hErecoNuE$mids, hist.fluct.Qr, type = "l", ylim = c(0, max(hist.fluct.Qr) * 1.1), ylab = "Entries", xlab = "Ereco [Mev]")
lines(hErecoNuE$mids, hist.fluct.Ql)
lines(hErecoNuE$mids, hist.fluct.Median, col = "red")

## Apply Oscillation

elist <- seq(0.1, 5, 0.01)

EneLen <- cbind(nue_data$Etrue/1000, nue_data$pathLen/1000)

pars <- c(sin2 = 0.7, dm2 = 50)
# pars <- c(sin2=1,dm2=2)

# baseline <- 0.280 # [km], nominal max
baseline <- 0.247  # [km], average

probs <- apply(EneLen, 1, function(x) Sterile.ProbNuEToNuE(x[1], x[2], pars = pars))

# probsLconst <- apply(EneLen,1,function(x) Sterile.ProbNuEToNuE(x[1],baseline,pars=pars))
probsLconst <- Sterile.ProbNuEToNuE(elist, baseline, pars = pars)


plot <- plot(EneLen[, 1], probs, type = "p", cex = 0.01, xlim = c(0, 5), ylab = "Probability", xlab = "True Energy [GeV]")
par(new = T)
plot <- plot(elist, probsLconst, type = "l", cex = 0.01, xlim = c(0, 5), ylab = "", xlab = "", main = "nue disappearance prob")


# Update Probs

# go here####
nue_data$OscProb <- probs

##

hErecoNuEOsc <- wtd.hist(main = "NuE Selection", xlab = labl.ereco, ylab = labl.entries, x = nue_data[nue_data$IsNuE == 1, ]$Ereco,
    weight = nue_data[nue_data$IsNuE == 1, ]$OscProb * nue_data[nue_data$IsNuE == 1, ]$weight, freq = T, breaks = bins.ereco)


ratioEreco <- hErecoNuEOsc$counts/hErecoNuE$counts

hEtrueNuEOsc <- wtd.hist(main = "NuE Selection", xlab = labl.etrue, ylab = labl.entries, x = nue_data[nue_data$IsNuE == 1, ]$Etrue,
    weight = nue_data[nue_data$IsNuE == 1, ]$OscProb * nue_data[nue_data$IsNuE == 1, ]$weight, freq = T, breaks = bins.etrue)



ratioEtrue <- hEtrueNuEOsc$counts/hEtrueNuE$counts


plot(hErecoNuE$mids, ratioEreco, type = "S")


plot(hEtrueNuE$mids, ratioEtrue, type = "S")


#

hObs.EtrueNuE <- hEtrueNuEOsc

base <- seq(0.21, 0.28, length.out = 200)

sin2 <- seq(0.5, 1, length.out = 100)

vobs <- hEtrueNuE$counts * Sterile.ProbNuEToNuE(hEtrueNuE$mids, 0.23, pars = pars)

nLL <- as.numeric()

for (i in 1:length(sin2)) {

    p <- c(sin2[i], pars[2])
    v <- hEtrueNuE$counts * Sterile.ProbNuEToNuE(hEtrueNuE$mids, 0.23, pars = p)

    # nLL[i] <- LogLike.Histo(hObs.EtrueNuE$counts,v)
    nLL[i] <- LogLike.Histo(vobs, v)
}

plot(sin2, nLL)



