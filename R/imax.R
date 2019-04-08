imax <-
function(x){seq(length(x))[x%in%max(x,na.rm=T)][1]}
