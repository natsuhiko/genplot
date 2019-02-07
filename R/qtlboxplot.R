qtlboxplot = function(g, y, ref="A", alt="B", rs=""){
	par(family="Liberation Sans")
	n = length(y)
	x = round(g)
	s = t(matrix(unlist(lapply(split(data.frame(seq(n),y),x),function(yy){yd=density(yy[,2]); c(rbind(yy[,1],approx(yd$x,yd$y,xout=yy[,2])$y))})),2))
	s = s[order(s[,1]),]
	s = s[,2]/max(s[,2])
	par(mgp=c(2,0.5,0), mar=c(3,3,1,1))
	boxplot(y~x, axes=F, xlab="", ylab="", at=0:2, outline=F, ylim=range(y))
	points(x+rnorm(n)/7*s, y, pch=20, col=col.ba[x+1])
	mtext(c(paste(ref,ref,sep=""),paste(ref,alt,sep=""),paste(alt,alt,sep="")), 1, at=0:2, line=0)
	mtext(rs, 1, line=1.5)
	axis(2,las=2)
}
