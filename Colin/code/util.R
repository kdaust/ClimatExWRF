

elements <- c("Tmin", "Tmax", "Pr")
monthcodes <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
monthdays <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 365.25)
month.abb.lowercase <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")


# function for preparing data
prep <- function(x, studyarea, element){
  x <- crop(x, studyarea)
  studyarea <- project(studyarea, x)
  x <- mask(x, studyarea)
  if(element=="Pr") values(x) <- log2(values(x))
  values(x)[!is.finite(values(x))] <- NA
  return(x)
}

# function for capping the values of a raster to the valid range of color breaks. 
bound <- function(x, breaks){
  values(x)[values(x)>max(breaks)] <- max(breaks)
  values(x)[values(x)<min(breaks)] <- min(breaks)
  return(x)
}

# function for adding filled points to a taylor diagram (plotrix doesn't accomodate `bg`)
taylor.diagram.filled <- function(ref, model, pch=19, col="black", bg="red", 
                                  cex=1.5, normalize=FALSE, ...) {
  require(plotrix)
  
  # reference standard deviation
  sd.ref <- sd(ref, na.rm=TRUE)
  sd.model <- sd(model, na.rm=TRUE)
  cc <- cor(ref, model, use="complete.obs")
  
  # ---- Normalization ----
  if (normalize) {
    sd.model <- sd.model / sd.ref
    sd.ref <- 1
  }
  
  # base diagram
  taylor.diagram(ref, model, pch=NA, normalize=normalize, ...)  # setup only
  
  # compute coordinates
  x <- sd.model * cc
  y <- sd.model * sin(acos(cc))
  
  # add points with bg support
  points(x, y, pch=pch, col=col, bg=bg, cex=cex)
}

