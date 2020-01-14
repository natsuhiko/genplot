tabix2charmat = function(f, reg){
	res = as.list(.Call("tabix2charmat", as.character(f), as.character(reg)))
	t(matrix(res[[2]], res[[1]][2]))
	#gtf = as.data.frame(, stringsAsFactors=F)
}
