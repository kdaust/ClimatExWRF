
# Plots showing the annual and monthly CONUS-II WRF climatology maps of tmin, tmax, and pr. 
# Colin Mahony colin.mahony@gov.bc.ca

library(terra)
library(data.table)
library(scales)
library(RColorBrewer)

monthcodes <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
month.abb.lowercase <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
monthdays <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 365.25)
elements <- c("tmin", "tmax", "pr")
element.names <- c("Tmin (\u00B0C)", "Tmax (\u00B0C)", "Precipitation (mm/day)")

colors <- c("simple", "extended")
color <- "simple"
for(color in colors){
  
  projections <- c("lambert", "latlon")
  
  projection <- "latlon"
  # for(projection in projections){
    
    e=3
    for(e in 1:3){
      element = elements[e]
      
      # load the ClimatEx WRF data for the variable
      dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CONUSII/"
      wrfconus2 <- rast(paste0(dir, paste("conus2_climatology_", c("tmin", "tmax", "pr")[e], "_", projection, ".tif", sep="")))
      if(e != 3) wrfconus2 <- wrfconus2 - 273.15
      
      studyarea <- ext(wrfconus2)-c(7, 5, 4, 3)
      
      png(filename=paste("results\\ClimatExEval.AnnualCycle.conus2", element, projection,color, "png",sep="."), type="cairo", units="in", width=if(projection == "latlon") 6.5 else 5.5, height=7, pointsize=12, res=600)
      
      mat <- matrix(c(2,3,4,5,13,1,1,6,12,1,1,7,11,10,9,8),4, byrow=T)   #define the plotting order
      layout(mat, widths=c(1,1,1,1), heights=c(1,1,1,1))   #set up the multipanel plot
      
      # define the upper and lower limits of daily precipitation
      if(e == 3){
        top <- 0.995
        bottom <- 0.005
        prlim <- vector()
        for(m in 1:12){
          X <- wrfconus2[[m]]/monthdays[m]
          X <- crop(X, studyarea)
            values(X) <- log2(values(X))
            values(X)[!is.finite(values(X))] <- NA
          lim <- as.numeric(unlist(quantile(values(X), c(top, bottom), na.rm=T)))
          prlim <- if(m==1) lim else c(prlim, lim)
        }
        
      }
      
      # clip the annual raster to the quantile limits
      X <- if(e == 3) wrfconus2[[13]]/monthdays[13] else wrfconus2[[13]]
      X <- crop(X, studyarea)
      if(e==3){
        values(X) <- log2(values(X))
        values(X)[!is.finite(values(X))] <- NA
      }
      top <- if(e == 3) 1 else 0.995
      lim1 <- quantile(if(e == 3) prlim else values(wrfconus2), top, na.rm=T)
      values(X)[which(values(X)>lim1)] <- lim1
      bottom <- if(e == 3) 0 else 0.005
      lim0 <- if(e == 3) -2.45 else quantile(values(wrfconus2), bottom, na.rm=T)
      values(X)[which(values(X)<lim0)] <- lim0
      values(X)[1:2] <- c(lim0, lim1)
      
      inc=0.01
      breaks=seq(lim0-inc, lim1+inc, inc)
      if(color=="extended"){
        ColScheme <- colorRampPalette(if(element=="pr") brewer.pal(9, "YlGnBu") else c(brewer.pal(5, "GnBu"), rev(brewer.pal(11, "RdYlBu")), rev(brewer.pal(5, "Greys"))))(length(breaks)-1)
      } else {
        ColScheme <- colorRampPalette(if(element=="pr") brewer.pal(9, "YlGnBu") else c("black", rev(brewer.pal(11, "RdYlBu")), "black"))(length(breaks)-1)
      }
      
      ## ANNUAL map
      par(mar=c(0.2, 0.2, 0.2, 0.2))
      plot(X, col=ColScheme, zlim = c(min(breaks), max(breaks)), mar = NA, axes = FALSE, frame = TRUE,  legend=F, xaxt="n", yaxt="n", maxcell = ncell(X))
      mtext("Annual", line=-1.5, adj=0.95, side=3, cex=1)
      # plot(bdy.bc.LCC, add=T, lwd=0.5, border=alpha("white", 0.25))
      box()
      
      ## ANNUAL legend
      if(projection=="latlon"){
        ramploc <- c(2.5, 4, 1, 8)
        boxloc <- c(1.5, 2, 0.6, 0.6)
        labelloc <- c(0.5, 0.5)
      } else {
        ramploc <- c(40, 75, 450, 670)
        boxloc <- c(35, 45, 20, 20)
        labelloc <- c(10, 10)
      }
      xl <- studyarea[1]+ramploc[1]; yb <- studyarea[3]+ramploc[3]; xr <- studyarea[1]+ramploc[2]; yt <- studyarea[3]+ramploc[4]
      rect(xl-boxloc[1],  yb-boxloc[3],  xr+boxloc[2],  yt+boxloc[4], col=alpha("white", 1), border="black")
      rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  border=NA, col=ColScheme)
      rect(xl,  yb,  xr,  yt)
      at <- (quantile(breaks, c(bottom, 0.25,0.5,0.75, top), na.rm=T))
      labels <- if(e==3) round(2^at,0) else round(at,0)
      text(rep(xr,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)-1)),labels,pos=4,cex=0.9,font=1, offset=0.3)
      text(xl-labelloc[1], mean(c(yb,yt))-labelloc[2], element.names[e], srt=90, pos=3, cex=0.85, font=2)
      
      ## MONTHLY maps
      for(m in 1:12){
        
        X <- if(e == 3) wrfconus2[[m]]/monthdays[m] else wrfconus2[[m]]
        X <- crop(X, studyarea)
        if(e==3){
          values(X) <- log2(values(X))
          values(X)[!is.finite(values(X))] <- NA
        }  
        values(X)[which(values(X)>lim1)] <- lim1
        values(X)[which(values(X)<lim0)] <- lim0
        values(X)[1:2] <- c(lim0, lim1)
        
        plot(X, col=ColScheme, zlim = c(min(breaks), max(breaks)), mar = NA, axes = FALSE, frame = TRUE,  legend=F, xaxt="n", yaxt="n", maxcell = ncell(X))
        mtext(month.abb[m], line=-1.5, adj=0.05, side=3, cex=0.8)
        # plot(bdy.bc.LCC, add=T, lwd=0.5, border=alpha("white", 0.25))
        box()
        print(m)  
      }
      
      dev.off()
      
      print(element)  
    }
    
    print(projection)
  # }
  
  print(color)
}


