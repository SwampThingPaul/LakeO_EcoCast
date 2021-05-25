## 
## Preliminary effort Lake O forecast
##
## Experimenting with EDDIE module 5 model
## https://github.com/MacrosystemsEDDIE/module5
##
## Code was compiled by Paul Julian
## contact info: pjulian@sccf.org

## BAD 
## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
## Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

## Libraries
# Data Wrangling
library(AnalystHelper);
library(plyr)
library(reshape2)

# Paths
wd="C:/Julian_LaCie/_GitHub/LakeO_EcoCast"
paths=c(paste0(wd,c("/Exports/","/Plots/","/Data/","/src/")))
# Folder.Maker(paths);#One and done. Creates folders in working directory.
export.path=paths[1]
plot.path=paths[2]
data.path=paths[3]


# Models ------------------------------------------------------------------
source("./src/eddie_mod5_models.R")

create_npz_inputs <- function(time, PAR = NULL, swr = NULL, temp = NULL) {
  
  # year <- lubridate::year(time[1])
  ndays <- round(as.numeric(difftime(time, time[1], units = "days")))
  
  # PAR calculations ----
  if(is.null(PAR) & !is.null(swr)) {
    PAR <- LakeMetabolizer::sw.to.par.base(swr)
  }
  
  if(is.null(PAR) & is.null(swr)) {
    PAR <- 0.5 * (540 + 440 * sin(2 * pi * ndays/365-1.4))
  }
  
  if(length(temp) == 2) {
    TEMP <- (temp[1] + temp[2] * sin(2*pi*ndays/365-1.4))
  } else {
    TEMP <- temp
  }
  
  out = data.frame(Day = ndays, PAR = PAR, TEMP = TEMP)

  return(out)
  
}

create_npz_inputs(time=seq(date.fun("2021-05-01"),date.fun("2021-05-12"),"1 days"),temp=rnorm(12,20))


# -------------------------------------------------------------------------
# Experimenting with EDDIE module 5 model
## https://github.com/MacrosystemsEDDIE/module5


# Parameters for NPZ model
parms <- c(
  maxUptake = 1.0, #day-1
  kspar=120, #uEinst m-2 s-1
  ksdin=0.5, #mmol m-3
  maxGrazing=1.0, # day-1
  ksphyto=1, #mmol N m-3
  pFaeces=0.3, #unitless
  mortalityRate=0.4, #(mmmolN m-3)-1 day-1
  excretionRate=0.1, #day-1
  mineralizationRate=0.1, #day-1
  Chl_Nratio = 1, #mg chl (mmolN)-1
  Q10 = 2,  #unitless
  addTEMP = 0, # added to temperature
  scaleNLOAD = 1, # multiplier for N loading
  refTEMP = 20
)

# Initial conditions for NPZ
yini <- c(
  PHYTO = 0.2, #mmolN m-3
  # ZOO = 0.4, #mmolN m-3
  # DETRITUS = 1, #mmolN m-3
  DIN = 27) #mmolN m-3

# Load parameters
# params_site_NP_model="https://raw.githubusercontent.com/MacrosystemsEDDIE/module5/main/data/params_site_NP_model.csv"
# site_parms <- read.csv(params_site_NP_model, fileEncoding = "UTF-8-BOM")
# yini_sites_NP_model="https://raw.githubusercontent.com/MacrosystemsEDDIE/module5/0e392a528a456e466729f6fb9789fe6478fad92c/data/yini_sites_NP_model.csv"
# site_yini <- read.csv(yini_sites_NP_model, fileEncoding = "UTF-8-BOM")


## 
npz_inp=data.frame(Date=seq(date.fun("2021-05-01"),date.fun("2021-05-25"),"1 days"))

npz_inputs=create_npz_inputs(time=npz_inp$Date,temp=rnorm(25,25))
times=1:nrow(npz_inputs)

res <- matrix(NA, nrow = length(times), ncol = 3)
colnames(res) <- c("time", "Phytoplankton", "Nutrients")
res[, 1] <- times
res[1, -1] <- c(yini)

for(i in 2:length(times)) {
out <- as.matrix(deSolve::ode(y = yini, times = times[(i-1):i], func = NP_model,
                              parms = parms, method = "ode45", inputs = npz_inputs))
# print((out[2, ]))
res[i, -1] <- out[2, c(2, 3)]
# print(res[i, ])
yini <- out[2, c(2:3)]
}

res <- as.data.frame(res)
res$time <- npz_inp$Date
res$Chla <- (res$Phytoplankton * 62) # Convert from mmol/m3 to ug/L # * 4.97 + 1.58
res$Nutrients <- res$Nutrients * 0.062 # Convert from mmol/m3 to mg/L
res <- res[, c("time", "Chla", "Nutrients")]

plot(Chla~time,res)
plot(Nutrients~time,res)




# Real Data ---------------------------------------------------------------
N.mw=14.0067

dates=date.fun(c("2020-05-01","2021-05-01"))

## WQ Data
parameters=data.frame(param=c("NH4","TOC","Chla","Chla","DO","TKN","NOx","pH","SRP","TP","Sal","Temp","TN","Turb","Temp"),
                      Test.Number=c(20,100,61,179,8,21,18,10,23,25,98,7,80,12,7))
parameters=subset(parameters,param%in%c("Chla","NH4","NOx","Temp"))

wq.dat=DBHYDRO_WQ(dates[1],dates[2],"LZ40",parameters$Test.Number)
wq.dat=merge(wq.dat,parameters,"Test.Number")
wq.dat$FLWY=WY(wq.dat$Date.EST)

wq.dat.xtab=dcast(wq.dat,Station.ID+Date.EST+FLWY~param,value.var="HalfMDL",mean)
wq.dat.xtab$NH4.mM_m3=(wq.dat.xtab$NH4/N.mw)*1000
wq.dat.xtab$NOx.mM_m3=(wq.dat.xtab$NOx/N.mw)*1000
wq.dat.xtab$DIN=with(wq.dat.xtab,NH4.mM_m3+NOx.mM_m3)

wq.dat.xtab$DIN*N.mw/1000
with(wq.dat.xtab,NOx+NH4)
## Sonde Data
DBKEY.val=data.frame(DBKEY=c(39942,39937),param=c("Chl","Temp"))
wq.sond=DBHYDRO_daily(dates[1],dates[2],DBKEY.val$DBKEY)
wq.sond=merge(wq.sond,DBKEY.val,"DBKEY")
wq.sond$Date.EST=date.fun(wq.sond$Date)
wq.sond.xtab=dcast(wq.sond,Station+Date.EST~param,value.var="Data.Value",mean)

plot(Chl~Date.EST,subset(wq.sond.xtab,Date.EST%in%seq(date.fun("2021-03-01"),date.fun("2021-05-01"),"1 days")))

tmp=subset(wq.sond.xtab,Date.EST%in%seq(date.fun("2021-03-23"),date.fun("2021-04-06"),"1 days"))

plot(Chl/62~Date.EST,tmp)

## Weather Data
DBKEY.val=data.frame(DBKEY=c(13080),param=c("RadT"));# Units: kW m^2
wx.dat=DBHYDRO_daily(dates[1],dates[2],DBKEY.val$DBKEY)
wx.dat$Date.EST=date.fun(wx.dat$Date)

plot(Data.Value~Date.EST,wx.dat)

### 
parms <- c(
  maxUptake = 0.25, #day-1
  kspar=120, #uEinst m-2 s-1
  ksdin=0.5, #mmol m-3
  maxGrazing=1.0, # day-1
  ksphyto=1, #mmol N m-3
  pFaeces=0.3, #unitless
  mortalityRate=0.5, #(mmmolN m-3)-1 day-1
  excretionRate=0.1, #day-1
  mineralizationRate=0.2, #day-1
  Chl_Nratio = 1, #mg chl (mmolN)-1
  Q10 = 2,  #unitless
  addTEMP = 0, # added to temperature
  scaleNLOAD = 1, # multiplier for N loading
  refTEMP = 20
)

# Initial conditions for NPZ
yini <- c(
  PHYTO = mean(tmp$Chl/62),  #mmolN m-3 Chlorophyll Data 4.4 converts ug/L to mmol N m-3
  # ZOO = 0.4, #mmolN m-3
  # DETRITUS = 1, #mmolN m-3
  DIN = mean(wq.dat.xtab$DIN[15:17],na.rm=T) #mmolN m-3
)

##
npz_inp=data.frame(Date=seq(date.fun("2021-03-23"),date.fun("2021-04-06"),"1 days"))
npz_inputs=create_npz_inputs(time=npz_inp$Date,swr=subset(wx.dat,Date.EST%in%npz_inp$Date)$Data.Value*1000,temp=tmp$Temp)
times=1:nrow(npz_inputs)

res <- matrix(NA, nrow = length(times), ncol = 3)
colnames(res) <- c("time", "Phytoplankton", "DIN")
res[, 1] <- times
res[1, c(2,3)] <- c(yini)

for(i in 2:length(times)) {
  out <- as.matrix(deSolve::ode(y = yini, times = times[(i-1):i], func = NP_model,
                                parms = parms, method = "ode45", inputs = npz_inputs))
  # print((out[2, ]))
  res[i, -1] <- out[2, c(2, 3)]
  # print(res[i, ])
  yini <- out[2, c(2:3)]
}

res <- as.data.frame(res)
res
res$time <- npz_inp$Date
res$Chla <- (res$Phytoplankton * 62) # Convert from mmol C m-3 to ug/L # 
res$DIN <- res$DIN * (N.mw/1000) # Convert from mmol N m-3 to mg/L
res <- res[, c("time", "Chla", "DIN")]
res$time=date.fun(res$time)

plot(Chla~time,res)
plot(DIN~time,res)

Model_period=seq(date.fun("2021-03-23"),date.fun("2021-05-01"),"1 days")

# png(filename=paste(plot.path,"SimpleSim_NPhytop.png",sep=""),width=5,height=3,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,3,0.5,0.5),oma=c(1,1,0.5,0.5),xpd=F)

xlim.val=date.fun(c("2021-03-23","2021-05-01"));xmaj=seq(xlim.val[1],xlim.val[2],"14 days");xmin=seq(xlim.val[1],xlim.val[2],"1 days")
ylim.val=c(0,30);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(Chla~time,res,xlim=xlim.val,ylim=ylim.val,type="n",ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col="grey")#adjustcolor("grey",0.5))
# with(res,points(time,Chla,pch=22,bg="grey",lwd=0.5))
# with(wq.sond.xtab,points(Date.EST,Chl,pch=21,bg="red",lwd=0.5))
with(res,pt_line(time,Chla,2,"grey",1,22,"grey"))
with(wq.sond.xtab,pt_line(Date.EST,Chl,2,"red",1,21,"red"))
axis_fun(1,xmaj,xmin,format(xmaj,'%m-%d'),line=-0.5)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=1,line=1.5,"Date (Month-Day)")
mtext(side=2,line=2.5,"Chlorophyll-a (\u03BCg L\u207B\u00B9)")
legend("topright",legend=c("Observed","Simulated"),
       pch=c(21,22),pt.bg=c("red","grey"),col=c("black"),
       lwd=0.05,lty=NA,
       pt.cex=1,ncol=1,cex=1,bty="n",xpd=NA,xjust=0.5,y.intersp=1,x.intersp=0.75)
dev.off()

ylim.val=c(0,0.5)
plot(DIN~time,res,xlim=xlim.val,ylim=ylim.val)
# with(wq.dat.xtab,points(Date.EST,(DIN*N.mw/1000),pch=21,bg="red"))
with(wq.dat.xtab,points(Date.EST,NOx,pch=21,bg="red"))




### 

parms <- c(
  maxUptake = 0.25, #day-1
  kspar=120, #uEinst m-2 s-1
  ksdin=0.5, #mmol m-3
  maxGrazing=1.0, # day-1
  ksphyto=1, #mmol N m-3
  pFaeces=0.3, #unitless
  mortalityRate=0.5, #(mmmolN m-3)-1 day-1
  excretionRate=0.1, #day-1
  mineralizationRate=0.2, #day-1
  Chl_Nratio = 1, #mg chl (mmolN)-1
  Q10 = 2,  #unitless
  addTEMP = 0, # added to temperature
  scaleNLOAD = 1, # multiplier for N loading
  refTEMP = 20
)

##
sdate=date.fun("2021-03-23")
npz_inp=data.frame(Date=seq(sdate,date.fun(sdate+lubridate::ddays(12)),"1 days"))
npz_inputs=create_npz_inputs(time=npz_inp$Date,swr=subset(wx.dat,Date.EST%in%npz_inp$Date)$Data.Value*1000,temp=subset(wq.sond.xtab,Date.EST%in%npz_inp$Date)$Temp)
times=1:nrow(npz_inputs)

# Initial conditions for NPZ
yini <- c(
  PHYTO = subset(wq.sond.xtab,Date.EST==sdate)$Chl/62,  #mmolN m-3 Chlorophyll Data 4.4 converts ug/L to mmol N m-3
  # ZOO = 0.4, #mmolN m-3
  # DETRITUS = 1, #mmolN m-3
  DIN = mean(wq.dat.xtab$DIN[15:17],na.rm=T) #mmolN m-3
)

sdate=date.fun(max(npz_inp$Date))

res <- matrix(NA, nrow = length(times), ncol = 3)
colnames(res) <- c("time", "Phytoplankton", "DIN")
res[, 1] <- times
res[1, c(2,3)] <- c(yini)

for(i in 2:length(times)) {
  out <- as.matrix(deSolve::ode(y = yini, times = times[(i-1):i], func = NP_model,
                                parms = parms, method = "ode45", inputs = npz_inputs))
  # print((out[2, ]))
  res[i, -1] <- out[2, c(2, 3)]
  # print(res[i, ])
  yini <- out[2, c(2:3)]
}

res <- as.data.frame(res)
res
res$time <- npz_inp$Date
res$Chla <- (res$Phytoplankton * 62) # Convert from mmol C m-3 to ug/L # 
res$DIN <- res$DIN * (N.mw/1000) # Convert from mmol N m-3 to mg/L
res <- res[, c("time", "Chla", "DIN")]
res$time=date.fun(res$time)


xlim.val=date.fun(c("2021-03-23","2021-05-01"));xmaj=seq(xlim.val[1],xlim.val[2],"14 days");xmin=seq(xlim.val[1],xlim.val[2],"1 days")
ylim.val=c(0,30);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)

# png(filename=paste(plot.path,"SimpleSim_NPhytop.png",sep=""),width=5,height=3,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,3,0.5,0.5),oma=c(1,1,0.5,0.5),xpd=F)

plot(Chla~time,res,xlim=xlim.val,ylim=ylim.val,type="n",ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col="grey")#adjustcolor("grey",0.5))
# with(res,points(time,Chla,pch=22,bg="grey",lwd=0.5))
# with(wq.sond.xtab,points(Date.EST,Chl,pch=21,bg="red",lwd=0.5))
with(res,pt_line(time,Chla,2,"grey",1,22,"grey"))
with(wq.sond.xtab,pt_line(Date.EST,Chl,2,"red",1,21,"red"))
axis_fun(1,xmaj,xmin,format(xmaj,'%m-%d'),line=-0.5)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=1,line=1.5,"Date (Month-Day)")
mtext(side=2,line=2.5,"Chlorophyll-a (\u03BCg L\u207B\u00B9)")
legend("topright",legend=c("Observed","Simulated"),
       pch=c(21,22),pt.bg=c("red","grey"),col=c("black"),
       lwd=0.05,lty=NA,
       pt.cex=1,ncol=1,cex=1,bty="n",xpd=NA,xjust=0.5,y.intersp=1,x.intersp=0.75)
dev.off()