
arplot <- function(x1, x2, z, add=F, mindist=0, th=0, lwd=1.5, xlim=NULL, ylim=c(-1,0), collim=c("#e55f5300","#e55f53ff"), ArrowHead=T, ArrowSize=1){
  if(is.null(xlim)){xlim=range(c(x1,x2))}
  
  x=cos(225:315/180*pi);
  rx=1/diff(range(x))
  x=(x-min(x))/diff(range(x));
  y=sin(225:315/180*pi);
  ry=1/diff(range(y))
  y=(y-min(y))/diff(range(y))-1;
  
  if(is.matrix(x1)){
    d=c(abs((x1[,1])-(x2[,2])), abs((x1[,2])-(x2[,1])))
    r=abs(min(y)*max(d))
    if(!add){plot(1,1,type="n",ylim=ylim,xlim=xlim, axes=F); axis(1)}
    for(i in 1:nrow(x1)){
      if(x1[i,1] < x2[i,2]){
        a=x1[i,]
        b=x2[i,]
      }else{
        a=x2[i,]
        b=x1[i,]
      }
      xcoord = c(x*(b[2]-a[1])+a[1],rev(x)*(b[1]-a[2])+a[2])
      ycoord = c(y*(b[2]-a[1]),     rev(y)*(b[1]-a[2]))
      polygon(xcoord, ycoord/r, col=colorRampPalette(c(collim),alpha=T)(100)[ceiling(z[i]*100)], border=NA)
    }
  }else{
    len=diff(xlim)*0.03*ArrowSize
    d=abs(x1-x2)
    r=abs(min(y)*max(d))
    if(!add){plot(1,1,type="n",ylim=ylim,xlim=xlim); }
    for(i in 1:length(x1)){
      coli=colorRampPalette(c(collim),alpha=T)(100)[ceiling(z[i]*100)]
      if(x1[i]<x2[i] && d[i]>mindist && z[i]>th){
        if(ArrowHead){
          lines(x[1:90]*d[i]+x1[i], y[1:90]*d[i]/r, col=coli, lwd=lwd)
          polygon(c(len*cos(200/180*pi),0,len*cos(250/180*pi),len*cos(200/180*pi))*rx+d[i]+x1[i], c(len*sin(200/180*pi),0,len*sin(250/180*pi),len*sin(200/180*pi))*ry/r ,col=coli, border=NA)
        }else{
          lines(x*d[i]+x1[i], y*d[i]/r, col=coli, lwd=lwd)
        }
      }else if(x1[i]>x2[i] && d[i]>mindist && z[i]>th){
        if(ArrowHead){
          lines(x[2:91]*d[i]+x2[i], y[2:91]*d[i]/r, col=coli, lwd=lwd)
          polygon(c(-len*cos(200/180*pi),0,-len*cos(250/180*pi),-len*cos(200/180*pi))*rx+x2[i], c(len*sin(200/180*pi),0,len*sin(250/180*pi),len*sin(200/180*pi))*ry/r ,col=coli, border=NA)
        }else{
          lines(x*d[i]+x2[i], y*d[i]/r, col=coli, lwd=lwd)
        }
      }
    }
  }
}
