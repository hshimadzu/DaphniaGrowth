Supplement for: Allometric growth model reveals temperature
unpredictability impacts growth, reproduction, and longevity of
*Daphnia*
================
Hideyasu Shimadzu and Miguel Barbosa

This supplement file provides some lines of essential R code used in the
paper for illustration purposes.

### Data preparation

``` r
data0 <- read.csv("DaphniaData.csv")

# temp treatments
cond <- c("mean", "high", "var", "low")
# IDs
cond.ids <- unique(data0[,c("ID", "f1_treatment")])
ids <- cond.ids$ID
# time
dt <- 0.5
time <- seq(0, max(data0$time, na.rm=T), by=dt)
```

## Parameter estimation

### Required packages

``` r
library(gam)
```

### Estimating ![\hat{l}\_{it}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7Bl%7D_%7Bit%7D "\hat{l}_{it}") and ![d\hat{l}\_{it}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d%5Chat%7Bl%7D_%7Bit%7D "d\hat{l}_{it}") *via* loess smoothing

``` r
#######################
### fitting a loess curve
#######################
len <- list()
for(i in 1:length(ids)){
    len[[i]] <- try(loess(length~time, data=subset(data0, ID==ids[i])))
}

#######################
### predicting a loess curve
#######################
p.l <- lapply(len, function(x)try(predict(x, newdata=data.frame(time=time))))
loc1 <- sapply(p.l, class)=="numeric" # brunches which have the predicted values.

######################
### calculating derivatives
#######################
p.dl <- lapply(p.l, function(x)try(diff(x)/dt))
```

This fitting process produces error messages because some individuals’
body-length records are too short for loess to fit/predict the model.
The code above remove such individuals from the following parameter
estimation procedure. A total of 569 individuals (High: 146; Rearing:
148; Low: 140; and Unpredictable: 135) are used.

### Estimating parameters

These are the lines of the essential code for parameter estimation. The
integration is carried out *via* the trapezium method.

``` r
## prep weights for trapezim integration
wgt <- lapply(p.dl, function(x)as.numeric(!is.na(x))*dt)
for(i in 1:length(wgt)){
  wgt[[i]][1] <- dt/2
  wgt[[i]][sum(!is.na(p.dl[[i]]))] <- dt/2
}

## specified lambda
lambda <- 0.179

#################
## Parameter estimation
#################
fit.l <- list()
for(i in 1:length(cond)){
  loc2 <- cond.ids$f1_treatment==cond[[i]]
  select.loc <- c(1:length(p.l))[loc1&loc2]
  p.dl0 <- lapply(p.dl[select.loc], function(x)data.frame(time=time[1:length(x)],length=x))
  tim <- unlist(lapply(p.dl0, function(x)x$time))
  pdl <- unlist(lapply(p.dl0, function(x)x$length))
  wg0 <- wgt[select.loc]
  wg <- unlist(wg0)[!is.na(pdl)]
  tim <- tim[!is.na(pdl)]
  pdl <- pdl[!is.na(pdl)]
  pl <- unlist(lapply(p.l[select.loc], function(x)x[-1])); pl <- pl[!is.na(pl)]

  X <- model.matrix(~-1+as.factor(tim))
  Z <- X*pl

  lm1 <- lm(pdl~-1+I(Z^lambda), weights=wg)
  r <- data.frame(time=unique(tim), r=coef(lm1))
  rownames(r) <- NULL

fit.l[[i]] <- list(r, lm1)
}

#################
## Calculating r_t
#################
for(i in 1:length(fit.l)){
fit.l[[i]][[1]]$rt <- 1-fit.l[[i]][[1]]$r/fit.l[[i]][[1]]$r[1] # q = fit.l[[i]][[1]]$r[1]
}
```

The confidence intervals for the parameters can be obtained *via*
bootstrap methods, repeating the above process thousands times.

## Predicting body length

### Required package

``` r
library(deSolve)
```

### Calculating body length ![\hat{l}\_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7Bl%7D_t "\hat{l}_t")

``` r
#######################
### Growth model
#######################
Lmodel <- function(t, state, parms){
  with(as.list(c(state, parms)),{
    r <- rt(t)
    dl <- r0*(1-r)*l^a
  list(dl=dl, r=r)
  })
}

#######################
### Calculating l_t
#######################
trj <- list()
for(i in 1:4){
  rt <- approxfun(fit.l[[i]][[1]][c("time", "rt")], rule = 2)
  t.max <- max(fit.l[[i]][[1]]$time)
  times <- seq(from=0, to=t.max, by=0.5)
  parms <- c(a=lambda)
  loc2 <- cond.ids$f1_treatment==cond[[i]]
  l0 <- mean(unlist(lapply(p.l[loc1&loc2], function(x)x[1])))
  state <- c(l=l0) ###########
  r0 <- fit.l[[i]][[1]][1, "r"]

trj[[i]] <- data.frame(ode(y=state, times=times, func = Lmodel, parms=parms))
}
```

## Figures

### Required package and settings

DaphniaBoots.RData file contains the bootstrap results to illustrate
confidence intervals.

``` r
load("DaphniaBoots.RData")

library(RColorBrewer)
library(scales)

condName <- c("Rearing", "High", "Unpredictable", "Low")
colour <- c("#FFA500", "#FF0000", "#00FF00", "#0000FF")
```

### Figure 1

``` r
par(mar = c(4, 4, 0, 1))
with(q.dframe, plot(temp, q, type="b", bty="n", pch=19, xlab="Temperature [C]", ylab="q", ylim=c(0.1, 0.22), xaxt="n", las=1))
axis(1, c(15, 20, 25,30), c(as.character(c(15, 20, 25)), "15-25"))
for(i in 1:4){
  segments(q.dframe[i,"temp"], q.dframe[i,"LCI"], q.dframe[i,"temp"], q.dframe[i,"UCI"])
}
```

![](DaphniaGrowthModelling_files/figure-gfm/Fig1-1.png)<!-- -->

### Figure 2

``` r
par(mfrow=c(2, 1))
par(mar = c(4, 4, 0, 0))
## plotting trajectories
plot(trj.CI[[1]][,1:2], xlim=c(0, 140), ylim=c(0.5, 6), bty="n", type="n", ylab="Length [mm]", xlab="", las=1)
for(i in 1:4){
  tmp <- trj.CI[[i]][!is.na(trj.CI[[i]]$Length),]
  with(tmp, lines(time, Length, lwd=2, col=colour[i]))
  with(tmp, polygon(x=c(time, rev(time)), y=c(LCI, rev(UCL)), col=alpha(colour[i], 0.2), border=NA))
}
legend("topleft", legend=c("High", "Rearing", "Low", "Unpredictable"), lty=1, lwd=2, col=c("#FF0000", "#FFA500", "#0000FF", "#00FF00"), bty = "n")

## Plotting rt
plot(fit.l[[1]][[1]][,c(1,3)], xlim=c(0, 140), ylim=c(0, 1), bty="n", type="n", xlab="Days", ylab=expression(r[t]), las=1)
for(i in 1:4){
lines(rt.CI[[i]][,1:2], col=colour[i], lwd=2)
with(rt.CI[[i]], polygon(x=c(time, rev(time)), y=c(LCI, rev(UCL)), col=alpha(colour[i], 0.2), border=NA))
}
```

![](DaphniaGrowthModelling_files/figure-gfm/Fig2-1.png)<!-- -->

### Figure 3

``` r
par(mfrow=c(2,2))
par(mar = c(4, 4, 1, 0))
for(i in c(2, 1, 4, 3)){
  condata <- subset(data0, f1_treatment==cond[i])
  with(condata, plot(time, length, col = rgb(0, 0, 0, 0.07), main=condName[i], pch=20, ylim=c(0, 6.3), xlim=c(0, 140), bty="n", xlab="Days", ylab="Length [mm]", las=2))
  with(trj.CI[[i]], lines(time, Length, col=colour[i], lwd=2))
}
```

![](DaphniaGrowthModelling_files/figure-gfm/Fig3-1.png)<!-- -->

### Figure 4

``` r
par(mfrow=c(1,1))
par(mar = c(4, 4, 0, 0))
with(data0, plot(time, cum.n, type="n", bty="n", xlab="Days", ylab="Total number of neonates", las=1))
for(i in 1:4){
  with(gt[[i]], lines(time, coef(lmlist[[i]])*gt, lwd=2, col=colour[i]))
  with(gt[[i]], polygon(x=c(time, rev(time)), y=c(CIl, rev(CIu)), col=alpha(colour[i], 0.2), border=NA))
  ttau <- tau[names(tau)==cond[i]]
  ggt <- coef(lmlist[[i]])*gt[[i]][time==ttau, "gt"]
  segments(ttau, 0, ttau, ggt, col=colour[i], lwd=2, lty=2)
  segments(ttau, ggt, 0, ggt, col=colour[i], lwd=2, lty=2)
}
legend("topleft", legend=c("High", "Rearing", "Low", "Unpredictable"), lty=1, lwd=2, col=c("#FF0000", "#FFA500", "#0000FF", "#00FF00"), bty = "n")
```

![](DaphniaGrowthModelling_files/figure-gfm/Fig4-1.png)<!-- -->

### Figure 5

``` r
par(mfrow=c(2,2))
par(mar = c(4, 4, 1, 0))
for(i in c(2, 1, 4, 3)){
  with(subset(data0, f1_treatment==cond[i]&is.na(DL)), plot(time, cum.n,
   col=1, pch=20, xlim=c(0, 140), ylim=c(0, 500), main=condName[i], bty="n", xlab="Days", ylab="Total number of neonates", las=1))
  with(gthat[[i]], lines(time, coef(lmlist[[i]])*gt, col=colour[i], lwd=2))
  ttau <- tau[names(tau)==cond[i]]
  ggt <- coef(lmlist[[i]])*gt[[i]][time==ttau, "gt"]
  segments(ttau, 0, ttau, ggt, col=colour[i], lwd=2, lty=2)
  segments(ttau, ggt, 0, ggt, col=colour[i], lwd=2, lty=2)
}
```

![](DaphniaGrowthModelling_files/figure-gfm/Fig5-1.png)<!-- -->
