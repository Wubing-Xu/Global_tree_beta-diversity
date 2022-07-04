#A function to make scatter plot between two varialbes, including distance matrixs

scatter.plot <- function(x,y,coords=NULL,all.points=TRUE,npoints=10000,mantel.test=FALSE,modified.ttest=FALSE,col=adjustcolor("black",alpha=0.1),
	reg.line=TRUE,cor.coef=TRUE,pvalue=TRUE,pvalue.symbol=FALSE,r2=FALSE,text.pos="topleft",xaxis=TRUE,yaxis=TRUE,pch=16,cex.lab=1.5,cex=1,cex.axis=1,lwd.axis=1,
	col.line="black",lwd.line=2,lty.line=NA,text.cex=1.3,...){
	if(all.points) id <- which(!is.na(x) & !is.na(y))
	if(!all.points) {set.seed(100); id <- sample(which(!is.na(x) & !is.na(y)),npoints)}
	
	plot(y[id]~x[id],col=col,pch=pch,cex.lab=cex.lab,cex=cex,xaxt="n",yaxt="n",...)
	if(xaxis) axis(side=1,lwd=lwd.axis,cex.axis=cex.axis)
	if(yaxis) axis(side=2,lwd=lwd.axis,cex.axis=cex.axis)
	
	if(cor.coef){
		if(mantel.test) {
			xy.mant <- mantel(x,y,na.rm=TRUE)
			r <- xy.mant$statistic
			p <- xy.mant$signif
		}
		if(modified.ttest) {
			xy.cor <- modified.ttest(x=x,y=y,coords=coords,nclass=NULL)
			r <- xy.cor$corr
			p <- xy.cor$p.value
			if(p<0.001) p <- 0.001
		}	
		if(!mantel.test & !modified.ttest) {
			xy.cor <- cor.test(x,y)
			r <- xy.cor$estimate
			p <- xy.cor$p.value
			if(p<0.001) p <- 0.001
		}
	
		g.text <- function(p,r,r2=TRUE){
			if(r2){
				if(p<= 0.001) mytext <- bquote(atop(italic(R)^2*" = "*.(round(r^2,3))*"***"))
				if(p<= 0.01 & p> 0.001) mytext <- bquote(atop(italic(R)^2*" = "*.(round(r^2,3))*"**"))
				if(p<= 0.05 & p> 0.01) mytext <- bquote(atop(italic(R)^2*" = "*.(round(r^2,3))*"*"))
				if(p> 0.05) mytext <- bquote(atop(italic(R)^2*" = "*.(round(r^2,3))^"ns"))			
			}
			if(!r2){
				if(p<= 0.001) mytext <- bquote(atop(italic(r)*" = "*.(round(r,3))*"***"))
				if(p<= 0.01 & p> 0.001) mytext <- bquote(atop(italic(r)*" = "*.(round(r,3))*"**"))
				if(p<= 0.05 & p> 0.01) mytext <- bquote(atop(italic(r)*" = "*.(round(r,3))*"*"))
				if(p> 0.05) mytext <- bquote(atop(italic(r)*" = "*.(round(r,3))^"ns"))			
			}
			return(mytext)
		}
		if(pvalue & pvalue.symbol & r2) my.text <- g.text(p=p,r=r,r2=TRUE)
		if(pvalue & pvalue.symbol & !r2) my.text <- g.text(p=p,r=r,r2=FALSE)
		if(!pvalue & !r2) my.text <- bquote(italic(r)*" = "*.(round(r,3)))
		if(!pvalue & r2) my.text <- bquote(italic(R)^2*" = "*.(round(r^2,3)))
		if(pvalue & p > 0.001 & !pvalue.symbol & !r2) my.text <- bquote(atop(italic(r)*" = "*.(round(r,3)),italic(P)*" = "*.(round(p,3))))
		if(pvalue & p <= 0.001 & !pvalue.symbol & !r2) my.text <- bquote(atop(italic(r)*" = "*.(round(r,3)),italic(P)*" < "*.(round(p,3))))
		if(pvalue & p > 0.001 & !pvalue.symbol & r2) my.text <- bquote(atop(italic(R)^2*" = "*.(round(r^2,3)),italic(P)*" = "*.(round(p,3))))
		if(pvalue & p <= 0.001 & !pvalue.symbol & r2) my.text <- bquote(atop(italic(R)^2*" = "*.(round(r^2,3)),italic(P)*" < "*.(round(p,3))))
		
		if(any(text.pos=="bottomright"))
			text(max(x[id],na.rm=T)-0.02*(max(x[id],na.rm=T)-min(x[id],na.rm=T)),min(y[id],na.rm=T)+0.05*(max(y[id],na.rm=T)-min(y[id],na.rm=T)),my.text,pos=2,cex=text.cex)
		if(any(text.pos=="topleft"))
			text(min(x[id],na.rm=T)+0.05*(max(x[id],na.rm=T)-min(x[id],na.rm=T)),max(y[id],na.rm=T)-0.05*(max(y[id],na.rm=T)-min(y[id],na.rm=T)),my.text,pos=4,cex=text.cex)
		if(all(is.numeric(text.pos)))
			text(text.pos[1],text.pos[2],my.text,cex=text.cex)
	}
	
	if(reg.line){
		if(is.na(lty.line)){
			if(cor.coef & p>=0.05) lty.line=2
			else lty.line=1
		}
		ypre <- predict(lm(y[id] ~ x[id]))
		xobs <- x[id]
		id.order <- order(xobs)
		lines(xobs[id.order],ypre[id.order],col=col.line,lwd=lwd.line,lty=lty.line)	
	}
}
