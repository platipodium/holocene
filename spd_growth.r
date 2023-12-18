# calculates Summed Probability Density (SPD)
#   and related growth (RGR) for each region and time slice
#   using an array of methods
#
# kai wirtz (HZG) 2023
#

rm(list=ls())
library(rcarbon)  # RCARBON by Crema2017
library('R.matlab')
source("movavg.r")

cc='Europe'
##breaks=seq(3400,9800,400)
breaks=seq(3000,9800,400)
# binning size in SPD calculation
binv=c(50,50,100,100,150,150)
binv=c(100,100)
ntag=length(binv)
nbrk=length(breaks)
myargs = commandArgs(trailingOnly=TRUE)
print(paste('args:',length(myargs),' nbrk=',nbrk))
if (length(myargs)>0) {
 an = as.numeric(myargs[1])
 tagi0= 1+floor(an/nbrk)
 tagi1= tagi0
 ti0 =  1+(an%%nbrk)
 ti1 = ti0
 } else { #all sub-domains
 ti0=1
 ti1=nbrk
 tagi0=1
 tagi1=ntag
 }

# input/output directory
scdir='p3k14c/'
# list of SPD methods
tagv=c('_NoNorm_Bin100') #c(,'_NoNorm_Bin100','_NoNorm50','_NoNorm150'

# read matlab C14 data
norms=c('','No')
dat<-readMat(paste0('c14mat/C14_europe_neo.mat'))
#'lonsn','latsn','C14agesn','C14SDsn','SiteIDsn','datIDsn'

ncol=2
nrow=3 # figure outlay
nmaxclu = 45
dtl = 325
dt  = 25
colv=rainbow(5)
#if (!exists("ymv"))

# loop over SPD methods
for (tagi in tagi0:tagi1)
  {
  nplot=0 # reset of variables
#  ci=NULL, yms=NULL    tms=NULL
#  tag=tagv[tagi]
  tag=paste0('_',norms[1+(tagi%%2)],'Norm_Bin',binv[tagi])
  print(paste(tagi,tag,ti0,ti1))
  ymsa=NULL
  # loop over time slices
  #tii=1
  for (tii in seq(ti0,ti1))
    {
    ti=breaks[tii]
    #tis=breaks[ti+1]
    # read region/cluster info
    file <- paste0(scdir,'bin/clusti3_', cc[1],ti,'_120.bin')
    load(file) #=file,clusti,k,wi,cluc
    clustn=seq(k)
    nregions=k # number of regions
    nplot=0
    # loop over regions
    # load regionalized C14 dates
#ymv=array(NaN,c(nmaxclu,length(breaks),(2*dtl)/dt))

ymv =array(NaN,c(nregions,ceiling((2*dtl)/dt)))
rgrv=array(NaN,c(nregions,ceiling((2*dtl)/dt)-1))
n0=array(0,nregions);ndates=n0
dim(ymv)
if (FALSE){
  file = paste0(scdir,"mat/PrePop3",tag,'_',tii,'_',24,".mat")
  d<-readMat(file)
  i0=25
  if (d$nreg>32){
    file = paste0(scdir,"mat/PrePop3",tag,'_',tii,'_',32,".mat")
    d<-readMat(file)
    i0=33
    }
  print(paste0("read starting data from ", file))
   ymv=d$ymv
   rgrv[1:nrow(d$rgr),]=d$rgr
   n0=d$nsites
   ndates=d$ndates
  }
 else
  {
  i0=1
  }


    lat=round(dat$lat,digits=3)
    lon=round(dat$lat,digits=3)

    for (i in i0:nregions) #1:nregions)
      {
      if (i%%(ncol*nrow)==1 ) {  #& tagi==ntag new figure
        if (nplot>0) { dev.off() }
        print(paste('new figure',tii,i,nplot))
        file = paste0(scdir,"plots/spd_clu",tag,"_",ti,"_",nplot,".png")
        png(file,width=980,height=940,units='px')
        #x11(width=18,height=14)
        par(oma = c(1, 0, 1, 1), mar = c(0.1, 0.1, 0.14, 0.5), mfrow=c(nrow,ncol), cex.lab=1.5,cex.sub=0.5, cex.main=1.,cex.axis=2)
        nplot = nplot + 1
        }
      ii <- which(clusti == i)
      lati=round(dat$lat[ii],digits=3)
      loni=round(dat$lon[ii],digits=3)
      ci=c()
      sites=ci
      n0[i]=length(ii)
      for (ij in 1:n0[i]) {
        ni=which(lat==lati[ij] & lon==lati[ij])
        ci=c(ci,ni)
        sites=c(sites,rep(ij,length(ni)))
      }
      # extra region of interest (here SW Germany)
        # i3=which.min(abs(cd[,4]-8.8)+abs(cd[,5]-49))
        # ii <- which(clustv$tv == cd[i3,2]) lo=dat$lon[ii] la=dat$lat[ii]
        # i2 <- which(lo>8.2 & lo<9.56 & la>48.35 & la<49.6  )
      i2 <- which(!is.na(dat$C14agesn[ci]))
      #i2=i2[1:400]
      ii = ci[i2]
      #i1=1:min(400,length(ii))
      #ii=ii[i1]
      sites=sites[i2]
      ndates[i]=length(i2)

      # print number of valid C14 dates per region
      X <- sprintf('%2d %d sites valid dates: %d ',i,n0[i],length(ii))
      print(X)

      #  calibration using intcal20
      # running calibration over 3 cores
      eurodates <- calibrate(dat$C14agesn[ii],dat$C14SDsn[ii],calCurves='intcal20',ncores=4,normalised=(tagi%%2==0))
      print(paste('eurodates with normalised=',(tagi%%2==0), ' ready ...'))

      # creating time bins
      eurobins <- binPrep(sites=sites,ages=dat$C14agesn[ii],h=binv[tagi])

      # calculate SPD
      #for (dtle in c(dtl,dtl/2)) {
      timeRange <- ti+c(dtl,-dtl) #time range calBP, older date first   c(10200,3000)
      steps <- seq(timeRange[1],timeRange[2],-dt) #25 year slices
      tirgr=steps[2:(length(steps)-1)]
      clu_spd = spd(x = eurodates,bins=eurobins,timeRange=timeRange) #8700,3300
      clu_rgra= spd2rc(clu_spd,breaks = steps)
      rgr     = clu_rgra$roca
      #expnull <- modelTest(DK.caldates, errors=DK$C14SD, bins=DK.bins, nsim=100, timeRange=c(8000,4000), model="exponential",runm=100)
    #summary(expnull,type='roc') #plot(expnull,type='roc')
    #Test Difference between 5120 and 4920 cal BP
    #results=p2pTest(uninull,p1=5120,p2=4920) #results #$pval [1] 0.01980198
    # From the RCARBON doc: The binning process should hence be used with caution, and its implications should be explored via a sensitivity analysis. The function binsense() enables a visual assessment of how different cut-off values can modify the shape of the SPD. The example below explores six different values and show how the highest peak in the SPD changes as a function of h but the overall dynamics remains essentially the same.
    #  binsense(x=DK.caldates,y=DK$SiteID,h=seq(0,500,100),timeRange=c(8000,4000))

      y  = clu_spd$grid$PrDens ## *10
      ym = movavg(y, dt, type="n")

      # store into matrix
      mi = seq(round(0.5* dt),floor(length(clu_spd$grid$calBP)/ dt)* dt, dt)
      ym = ym[mi]
      ymsa=cbind(ymsa,ym)
      tm=clu_spd$grid$calBP #+ timscale/2
      tm2=tm[mi] #+ timscale/2

      #print(tm[1:3])
      #print(tm2[1:3])

      ymv[i,1:length(ym)]=ym # ymv[tagi,i,]=ym
      #tmv[i,tii,1:length(tm2)]=tm2
      rgrv[i,1:length(rgr)]=rgr # ymv[tagi,i,]=ym
      #trgr[i,tii,1:length(tirgr)]=tirgr
        # plot SPD with own method
      ##plot(clu_spd,type="simple",col="indianred",lwd=1,lty=1)
      plot(tm,y,col="indianred",lwd=1,lty=1)
      #plot(clu_spd,runm= timscale,add=TRUE,type="simple",lwd=1,lty=2) #using a rolling average for smoothing
      #di=di+1 } # for di
        #time grid
        #abline(v=seq(4,7)*1000, col="purple")
        # plot smoothed SPD, also of other methods
        #for (tai in 1:tagi) {  lines(tm,y,col=colv[2],lwd=1,lty=1)
      #for (tai in 1:2) {
      lines(tm2,ymv[i,],col=colv[1],lwd=2,lty=1)
      lines(tirgr,mean(ymv[i,])+5*rgrv[i,],col=colv[1],lwd=3,lty=2)

      text(ti-dtl*0.5, mean(y)*0.15, labels = paste(i,ti), cex = 3, col = NULL)
      text(ti, mean(y)*0.05, labels = length(ii), cex = 2, col = NULL)
      #text(ti+dtl*0.5, mean(y)*0.05, labels = c(as.character(cd[i,4]),':',as.character(cd[i,5])) , cex = 2, col = NULL)
      if((i)%%(ncol*nrow)==-1) {
          legend(x = "topleft",legend=tagv)
          }
      #tms=c(tms,tm)  #yms=c(yms,ym)
      #ci=c(ci,rep(i,length(ym)))
      if(i%%8==0 | i==nregions) {
        file = paste0(scdir,"mat/PrePop3",tag,'_',tii,'_',i,".mat")
        print(paste0("write data to ", file))
        writeMat(file, poptime=tm, ymv=ymv,trgr=tirgr,rgr=rgrv,nreg=nregions,nsites=n0,ndates=ndates)
        }
      } # for i

    # save population data as Matlab binary
    file = paste0(scdir,"mat/AllPop3",tag,'_',tii,".mat")
    print(paste0("write data to ", file))
    writeMat(file, poptime=tm2, ymv=ymv,trgr=tirgr,rgr=rgrv,nreg=nregions,nsites=n0,ndates=ndates) #poptime=tm

    print(paste('close figure',tii,i,nplot))
    dev.off()
    tii=tii+1
    } # for tii
  }
