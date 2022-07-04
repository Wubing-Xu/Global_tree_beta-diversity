library(raster)
library(tmap)
library(viridis)
library(fBasics)

## A function using tmap to make maps
# break.type: the way to divide continuous values and map them to color palette; one of :linear, log.linear, quantile
# mymidpoint: the way to decide the middle points of color palette; one of: auto, median, or a numeric value
# digits: the number of decimal places to round values and for labels
# limit: the percent of smallest and largest values to be used as lowest and highest breaks; used to avoid exteme values
mytmap <- function(myraster, myvalue=NULL, polygons, breaks=NULL, break.type="linear", mymidpoint="auto", digits=2, limit=0.00, 
                   midpoint=NULL, negative.adj=TRUE, labels=NULL, cols=NULL, style="cont", title="", main.title="", fill.col="gray90", borders.col="gray40",
                   legend.position=c(0.1,0.1), legend.text.size=1.0, inner.margins=0, outer.margins=0, frame=F, ...){
  
  # to determine breaks and cols based on settings
  if(!is.null(myvalue)) {
    x <- myvalue[!is.na(myvalue[, 2]), 2]
    if(break.type=="log.linear") x <- log2(x)
    
    if(!is.na(mymidpoint) & mymidpoint=="median"){
      if(is.null(breaks) & break.type=="quantile")
        breaks <- quantile(x, probs = c(0+limit,0.25,0.5,0.75,1-limit))
      if(is.null(cols))	cols <- viridis(option = "plasma",length(x),direction = -1) #199
    }	
    
    if(!is.na(mymidpoint) & mymidpoint=="auto"){
      if(is.null(breaks) & (break.type=="linear" | break.type=="log.linear"))
        breaks <- seq(quantile(x, probs=0+limit), quantile(x,probs=1-limit), length=5)
      if(is.null(cols))	cols <- viridis(option = "plasma",length(x),direction = -1) #199
    }	
    
    if(!is.na(mymidpoint) & is.numeric(mymidpoint)){
      if(is.null(breaks) & break.type=="quantile")
        breaks <- c(quantile(x[x<=mymidpoint],probs=c(0+limit,0.5)),mymidpoint,quantile(x[x>=mymidpoint],na.rm=T,probs=c(0.5,1-limit)))
      
      if(is.null(breaks) & (break.type=="linear"| break.type=="log.linear")){
        if(break.type=="log.linear") mymidpoint <- log2(mymidpoint)
        breaks <- c(seq(quantile(x[x<=mymidpoint], na.rm=T, probs=0+limit), mymidpoint, length=3),
                    seq(mymidpoint, quantile(x[x>=mymidpoint], na.rm=T, probs=1-limit), length=3)[-1])
      }
      if(is.null(cols))	cols <- divPalette(length(x), "RdYlBu")[length(x):1]  #199
    }
    
    if(break.type =="log.linear") breaks <- 2^breaks
    
    # round breaks
    breaks <- unique(round(breaks, digits)) 
    # define labels
    if(is.null(labels)) labels <- sprintf(paste("%.", digits, "f", sep=""), breaks)
    # define midpoint  
    # midpoint <- ifelse(!is.na(mymidpoint) & is.numeric(mymidpoint), mymidpoint, median(breaks))   
    
    # to adjust values when both negative and positive present. 
    # This is just for visualization because tm_raster assign negative and positive values into two sides of palettes, rather than follow the breaks 
    if(negative.adj){
      mins <- min(myvalue[,2],na.rm=TRUE)
      maxs <- max(myvalue[,2],na.rm=TRUE)		
      if(mins<0 & maxs>0){
        plus <- ceiling(abs(mins))
        breaks <- breaks + plus
        myvalue[,2] <- myvalue[,2] + plus	
      }
    }
    
    # pass the values to be mapping into a raster
    myvalue[,2] <- round(myvalue[,2], digits)
    myraster[myvalue[,1]] <- myvalue[,2]
  }
  
  tm_shape(shp=polygons) + 
    tm_fill(col=fill.col) +
    tm_shape(myraster) + 
    tm_raster(breaks=breaks, labels=labels, palette=cols, style=style, 
              title="", legend.reverse=TRUE, legend.is.portrait=T, midpoint=midpoint) +
    tm_shape(shp=polygons) + 
    tm_borders(lwd=0.01, col=borders.col) +
    tm_layout(title=title, main.title=main.title, inner.margins=inner.margins, outer.margins=outer.margins, frame=frame, 
              legend.position=legend.position, legend.text.size=legend.text.size,...)
}
