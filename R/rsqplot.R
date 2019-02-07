rsqplot <-
function(x, y, r, var.lead=NULL, ...){
	col.ld = c("#000080", "#18B6FF", "#76FF02", "#FF8714", "#FF3226", "#925c9e")
	if(is.null(var.lead)){
		ind = imax(y)
	}else{
		ind = seq(length(y))[var.lead]
	}
	col = col.ld[c(1:5,5)][floor(r^2*5)+1]
	col[ind]=col.ld[6]
	pch = rep(21, length(x))
	pch[ind] = 23
	genplot(x, y, ..., bg=col, col="#666666", pch=pch)
	axis(2,las=2)
}
