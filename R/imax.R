imax <-
function(x){seq(length(x))[x==max(x,na.rm=T)][1]}
