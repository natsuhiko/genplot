#' Plot with genes Function
#'
#' Given a set of minimum P-values under the different numbers of tests,
#' returns p-values adjusted using the Beta correction that
#' allows you to control the family wise error rates.
#' @param p a vector of the minimum P-values.
#' @param ntests the numbers of tests for the set of P-values.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param fixParam an indicator vector (length = 3) of 0 (estimate) or 1 (fix) for alpha, beta and gamma parameters.
#' @param trace printing convergence status
#' @examples
#' p=NULL
#' for(i in 1:100){p = c(p, min(runif(i)))}
#' b.adjust(p, 1:100)

Tick.label <-
function (x)
{
    nx = nchar(format(floor(x),scientific=F))
    x = format(x,scientific=F)
    x[nx > 3] = paste(substring(x[nx > 3], 1, nx[nx > 3] - 3), substring(x[nx > 3], nx[nx > 3] - 2, nx[nx > 3]), sep = ",")
    x[nx > 6] = paste(substring(x[nx > 6], 1, nx[nx > 6] - 6), substring(x[nx > 6], nx[nx > 6] - 5, nx[nx > 6]+1), sep = ",")
    x[nx > 9] = paste(substring(x[nx > 9], 1, nx[nx > 9] - 9), substring(x[nx > 9], nx[nx > 9] - 8, nx[nx > 9]+2), sep = ",")
    x
}

getLine1<-
function(gm, key){
    if(is.null(gm)||length(gm)==0||nrow(gm)==0){return(NULL)}
    gab=lapply(split(gm[,c(3:4)],key),function(x){c(min(x[[1]]),max(x[[2]]))})
    gname=names(gab)
    gab=data.frame(t(matrix(unlist(gab),2)))
    grank=rank( rank(gab[[1]])*max(rank(gab[[2]])+1) + rank(gab[[2]]), ties="first")
    names(grank)=gname
    return(grank)
}
getLine<-
function(gm){
    grank=getLine1(gm,gm$gid)
    if(is.null(grank)){return(NULL)}
    trank=NULL
    for(i in names(grank)){
        trank=c(trank, getLine1(gm[gm$gid==i,],gm[gm$gid==i,]$tid))
    }
    arank = grank[gm$gid]*max(trank+1) + trank[gm$tid]
    match(arank,sort(unique(arank)))
}

getNline <-
function(gm,l,a,b){
    if(is.null(gm)||length(gm)==0||nrow(gm)==0){return(0)}
    x=t(matrix(unlist(lapply(split(gm[,3:4],l),function(xx){c(min(xx[,1]),max(xx[,2]))})),2))
    if(nrow(x)==1){return(1)}
    x[,2]=x[,2]+(b-a)/7
    nline=0
    for(i in 1:(nrow(x)-1)){
        ll=(seq(nrow(x))-1)%%i+1
        iline=cumsum(c(1,diff(ll)<=0))
        flag=0
        for(j in 2:max(iline)){
            flag=flag+sum(x[iline==(j-1),2][1:sum(iline==j)]>x[iline==j,1])
            if(flag>0){break}
        }
        if(flag==0){return(i)}
    }
    return(nrow(x))
}

narm<-
function(x){
	x[!is.na(x)]
}

gmin<-
function(x,y){
    unlist(lapply(split(x,y), min, na.rm=T))
}

gmax<-
function(x,y){
    unlist(lapply(split(x,y), max, na.rm=T))
}

expandRange<-
function(x,r){
	x=range(x[!is.na(x)&x>(-Inf)&x<Inf])
	diff(x)*(r-1)/2
	x+c(-r,r)
}


# GRCh38 "/nfs/users/nfs_n/nk5/s109/Gencode/gencode.v24.annotation.gtf.sorted.gz"

genplot <-
function (x, y=NULL, chr=NULL, GTF="/nfs/users/nfs_n/nk5/s117/Gencode/gencode.v24.annotation.gtf.sorted.gz", MaxNTranscripts=0, GeneID=NULL, gname="all", biotype="all", LabGene=T, BED=NULL, xlim=NULL, ylim=NULL, ...)
{
    #BED=list("/nfs/users/nfs_n/nk5/s117/ATAC100/Peaks/Peak10M/peaksForPlot.gz", "/nfs/users/nfs_n/nk5/s117/Encode/Segmentations/segway.forPlot.gz")
    colAnnot=c(rgb(92,147,189,max=255),rgb(165,56,108,max=255), rgb(103,118,35,max=255), rgb(92,147,189,max=255), "#1B9E77","#7570B3","#D95F02","#E7298A","#66A61E","#E6AB02","#A6761D","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#666666", "#F17ABF", "#32A077", "#776CB2", "#166F2D", "#741D82", "#7B240B", "#BF007E", "#4F72B3")
    colAnnot=c("#f9a851","#b8b8b9","#925c9e","#e55f53","#0aaa87","#00aee0","#fcd804","#002d6d","#ac193c","#cdc3ba","#93bfe9","#555559")
    
    if(is.null(chr)){
	#invisible(get("plot",match("package:graphics",search()))(x, y, ...))
    }else{chr=paste("chr",gsub("chr","",as.character(chr)),sep="")}
    if(is.null(xlim)){
        if(is.vector(x)){
            xlim=expandRange(x, 1.05)
        }else{
            xlim=expandRange(x[,1], 1.05)
        }
    }#else{
    #    xlim=list(...)$xlim
    #}
    if(is.null(ylim)){
        if(!is.null(y)){
            ylim=expandRange(y,1.05)
        }else{
            if(is.vector(x)){
		y=rep(1,length(x))
                ylim=c(0.5,1.5);
            }else{
                ylim=expandRange(x[,2],1.05)
            }
        }
    }#else{
    #    ylim=list(...)$ylim
    #}
    
    
    
    # GTF filtering
    reg = paste(chr, ":", floor(xlim[1]), "-", ceiling(xlim[2]), sep = "")
    gtf = as.data.frame(as.list(.Call("loadGTF", as.character(GTF), reg)), stringsAsFactors=F)
    gtf = as.data.frame(as.list(.Call("loadGTF", as.character(GTF), as.character(gsub("chr","",reg)))), stringsAsFactors=F)
    if(!is.null(gtf)&length(gtf)>0){
        names(gtf)=c("gid", "tid", "start", "end", "strand", "source", "type","gname","biotype")
        if(biotype!="all"){gtf=gtf[gtf$biotype%in%biotype,]}
        if(gname!="all"){gtf=gtf[gtf$gname%in%gname,]}
        if(nrow(gtf)==0){
            nline=0
        }else{
            gtf.t=data.frame(names(gmin(gtf$start,gtf$tid)), gmin(gtf$start,gtf$tid), gmax(gtf$end,gtf$tid) )
            noltid=!(xlim[2]<gtf.t[[2]] | gtf.t[[3]]<xlim[1])
            gtf=gtf[gtf$tid%in%gtf.t[[1]][noltid>0],]
            if(nrow(gtf)==0){
                nline=0
            }else{
                if(MaxNTranscripts>0){
                    utid=narm(unlist(lapply(split(gtf$tid,gtf$gid),function(tid1){unique(sort(as.character(tid1)))[seq(MaxNTranscripts)]})))
                    gtf=gtf[gtf$tid%in%utid,]
                }
                # prep
                lall = getLine(gtf)
                if(is.null(lall)){
                    nline=0
                }else{
                    #print(lall)
                    nline=getNline(gtf,lall,xlim[1],xlim[2])
                    #print(nline)
                    if(LabGene){
                        nline=nline*2
                        lall=lall*2
                    }
                }
            }
        }
    }else{
        nline = 0;
    }
    
    
    # Annotations
    nannot = length(BED); #2;
    cat("N annot=");print(nannot)
    #BED=list("/nfs/users/nfs_n/nk5/s117/ATAC100/Peaks/Peak10M/peaksForPlot.gz", "/nfs/users/nfs_n/nk5/s117/Encode/Segmentations/segway.forPlot.gz")
    annot = lapply(BED, function(bed1){ as.data.frame(as.list(.Call("loadBed", as.character(bed1), as.character(reg))), stringsAsFactors=F) })
    #print(annot)
    
    
    # ploting main fig
    
    NT = 1
    nlheader=2
    nlaxis=3
    bottomline = nline + 0.2 + nannot
    yNmax=ylim[2]-ylim[1]
    
    par(family="Liberation Sans")
    par(mar=c(bottomline+nlaxis,2,nlheader,2), xpd=F, cex=1, cex.axis=1, cex.lab=1, mgp=c(2,.5,0))
    get("plot",match("package:graphics",search()))(x, y, axes=F, xaxs="i", yaxs="i", ylab="", xlab="", xlim=xlim, ylim=ylim, ...)
    #bottom=ylim[1]-(par()$din[2]-par()$fin[2]*NT)/par()$fin[2]/(bottomline+nlaxis+nlheader)*(bottomline)*yNmax
    #gwidth= (par()$din[2]-par()$fin[2]*NT)/par()$fin[2]/(bottomline+nlaxis+nlheader)*yNmax
    bottom=ylim[1]-(par()$fin[2]-par()$pin[2]*NT)/par()$pin[2]/(bottomline+nlaxis+nlheader)*(bottomline)*yNmax
    gwidth= (par()$fin[2]-par()$pin[2]*NT)/par()$pin[2]/(bottomline+nlaxis+nlheader)*yNmax
    hwidth = (1)*gwidth
    #yl=axis(2, label=F);
    ### axis 4
    ###yl=axis(4,line=-(sum(par()$mar[c(2,4)])+sum(par()$oma[c(2,4)]))/(par()$din[1]-par()$pin[1])*par()$pin[1], label=F)
    ###mtext(yl, 2, at=yl, las=2, line=-1, cex=0.8)
    xl=axis(1, label=F, line=bottomline);
    tpr=1
    if (diff(xlim)/1e+06 > 1) {
        xl=unique(sort(floor(xl/1e6)))*1e6
        Unit = "Mb"
        tpr=1e6
    } else if (diff(xlim)/1000 > 1) {
        Unit = "kb"
        xl=unique(sort(floor(xl/1000)))*1000
        tpr=1000
    } else {
        Unit = "bp"
    }
    axis(1, at = xl, lab = Tick.label(xl/tpr), line=bottomline, tick=F)
    mtext(paste("Chromosome ", gsub("chr", "", as.character(chr)), " position (", Unit, ")", sep = ""), 1, bottomline+1.5, cex=1)
    par(xpd=NA);clip(xlim[1], xlim[2], bottom,yNmax*NT*1.05) #1.05)
    axis(1, at = xl, lab = F, line=bottomline)
    segments(xlim[1],bottom,xlim[2],bottom)
    
    # ploting annotations
    if(nannot>0){
        for(j in 1:nannot){
            yannot = ylim[1] + (1 - j ) * gwidth - hwidth
            if(nrow(annot[[j]])>0){if(ncol(annot[[j]])>2){
                rect(annot[[j]][,1], yannot, annot[[j]][,2], yannot+gwidth*0.8, col=annot[[j]][,3], border=annot[[j]][,3])
            }else{
                rect(annot[[j]][,1], yannot, annot[[j]][,2], yannot+gwidth*0.8, col=colAnnot[j], border=colAnnot[j])
            }}
        }
    }
    

    # ploting genes
    if(nline>0){
        #par(xpd=NA)
        
        #clip(xlim[1],xlim[2],bottom,yNmax*NT)
        arrowwid=gwidth*0.4
        arrowlen=(par()$fin[2]-par()$pin[2]*NT)/(bottomline+nlaxis+nlheader)/par()$pin[1]*diff(xlim)
        
        # coordinates for exons
        utr = as.matrix(gtf[gtf$type%in%c("exon"),3:4])
        cds = as.matrix(gtf[gtf$type%in%c("CDS"),3:4])
        lutr= lall[gtf$type%in%c("exon")]-1
        lcds= lall[gtf$type%in%c("CDS")]-1
        ycds = ylim[1]-((c(lcds)%%nline)+(bottomline-0.2)-nline) * gwidth - hwidth
        yutr = ylim[1]-((c(lutr)%%nline)+(bottomline-0.2)-nline) * gwidth - hwidth
        
        # coordinates for transcripts
        orient =   unlist(lapply(split(gtf$strand,lall),function(xx){xx[1]}))
        tlab =     unlist(lapply(split(gtf$tid,lall),function(xx){xx[1]}))
        glab =     unlist(lapply(split(gtf$gid,lall),function(xx){xx[1]}))
        nlab =     unlist(lapply(split(gtf$gname,lall),function(xx){xx[1]}))
        #ygene = -((unlist(lapply(split(lall,         lall),function(xx){xx[1]}))-1)%%nline) * gwidth - hwidth - ifelse(nested, gwidth/5, 0)
        ygene = ylim[1]-(((unlist(lapply(split(lall,         lall),function(xx){xx[1]}))-1)%%nline)+(bottomline-0.2)-nline) * gwidth - hwidth
        x0 =       unlist(lapply(split(gtf[,3],lall),function(xx){min(xx)})) - as.numeric(orient==1)*arrowlen
        x1 =       unlist(lapply(split(gtf[,4],lall),function(xx){max(xx)})) + as.numeric(orient==0)*arrowlen
        y0 = ygene+gwidth*0.4
        
        # plotting
        gcol=rep("#074987", nrow(gtf))
	    if(!is.null(GeneID)){  gcol[gtf$gid%in%GeneID | gtf$gname%in%GeneID | gtf$tid%in%GeneID]="#FF0000"  }
        gcol.text=rep("#074987", length(nlab))
	    if(!is.null(GeneID)){  gcol.text[nlab%in%GeneID | glab%in%GeneID | tlab%in%GeneID]="#FF0000"  }
        segments(x0, y0, x1, y0, col=gcol.text)
        segments(c(x0[orient==1]+arrowlen/2,x0[orient==1]+arrowlen/2,x1[orient==0]-arrowlen/2,x1[orient==0]-arrowlen/2),
                 c(y0[orient==1]+arrowwid,y0[orient==1]-arrowwid,y0[orient==0]+arrowwid,y0[orient==0]-arrowwid),
                 c(x0[orient==1],x0[orient==1],x1[orient==0],x1[orient==0]),
                 c(y0[orient==1],y0[orient==1],y0[orient==0],y0[orient==0]), col=c(gcol.text[orient==1],gcol.text[orient==1],gcol.text[orient==0],gcol.text[orient==0]))
        rect(utr[,1], yutr, utr[,2], yutr+gwidth*0.8, border=gcol[gtf$type%in%c("exon")], col="white")
        rect(cds[,1], ycds, cds[,2], ycds+gwidth*0.8, border=NA, col=gcol[gtf$type%in%c("CDS")])
        

        x0[x0<xlim[1]]=xlim[1]
        x0[x0>xlim[2]]=xlim[2]
	if(LabGene){if(MaxNTranscripts>1){
            text(x0-arrowlen/2,y0+gwidth*0.8,paste(nlab,tlab,sep=" : "),pos=4,cex=1.0,col=gcol.text)
        }else{
            text(x0-arrowlen/2,y0+gwidth*0.8,nlab,pos=4,cex=1.0,col=gcol.text)
        }}
    }
    clip(xlim[1]-diff(xlim), xlim[2]+diff(xlim), bottom,yNmax*NT*1.05+diff(ylim))
    
    invisible(gtf)
}
