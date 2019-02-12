covplot <-
function(chr, xlim, snpreg="1:12345-12345", vcf="/path/to/your/tabixed.vcf.gz", bedg="/path/to/your/bedgraph.txt", GTF="/path/to/your/tabixed.gtf.gz", MaxNT=1, GM12878=F, LINE=F, BED=list()){
	
	col=c(rgb(25,50,75,max=255),rgb(50,201,233,max=255),rgb(245,26,87,max=255))
	a=xlim[1]
	b=xlim[2]
	y=NULL
	Ksmooth<-function(x,y,bw=50){
        	x.dens=density(x,weight=y/sum(y,na.rm=T),bw=bw,from=min(x),to=max(x))
        	return(cbind(x.dens$x,x.dens$y*sum(y)))
	}
	N = length(scan(bedg,""))
	tmp = tempfile()
	com = paste("tabix ", vcf, " ", snpreg, " | cut -f 10- > ", tmp, sep="")
	system(com)
	gen = unlist(lapply(strsplit(scan(tmp,""),":"),function(x){eval(parse(text=gsub("\\|","+",x[1])))}))
	if(length(gen)==0){gen=rep(0,N)}
	gen = gen[1:N]
	unlink(tmp)
	for(i in scan(bedg, "")){
		print(i)
		tmp = tempfile()
		com = paste("tabix ", i, " ", chr, ":", a, "-", b, " > ", tmp, sep="")
		system(com)
		x = read.table(tmp)
		x[1,2]=a
		x[nrow(x),3]=b+1
		y=cbind(y, rep(x[[4]], x[,3]-x[,2]))
	}

	z=NULL
	for(i in 0:2){
		if(sum(gen==i)>0){
			z = cbind(z, apply(y[,gen==i,drop=F],1,mean))
		}else{
			z = cbind(z, rep(0,nrow(y)))
		}
	}
	genplot(a:b, z[,1], type="l", chr=chr, xlim=xlim, col=NA, ylim=c(0,max(z)), biotype="all", GTF=GTF, MaxNT=MaxNT, BED=BED)
	rz = c(1,3)[rank(-apply(z[,-2],2,mean,na.rm=T))]
	rz = c(rz[1],2,rz[2])
	if(LINE){
		for(i in rz)lines(a:b, z[,i], col=col[i], lwd=1, type="s")
		abline(0,0)
	}else{
		for(i in rz)polygon(c(a,a:b,b), c(0,z[,i],0), col=col[i], border=col[i], lwd=2)
	}
}
