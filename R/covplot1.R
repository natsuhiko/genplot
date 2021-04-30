covplot1=
function (chr, xlim,
    bedg = "/path/to/your/bedgraph.txt",
    MaxNT = 1, GM12878 = F, LINE = F, BED = list(), biotype="protein_coding", gname)
{
print(chr)
print(xlim)
print(bedg)
    col = c(rgb(25, 50, 75, max = 255), rgb(50, 201, 233, max = 255),
        rgb(245, 26, 87, max = 255))
    a = xlim[1]
    b = xlim[2]
    y = NULL
    Ksmooth <- function(x, y, bw = 50) {
        x.dens = density(x, weight = y/sum(y, na.rm = T), bw = bw,
            from = min(x), to = max(x))
        return(cbind(x.dens$x, x.dens$y * sum(y)))
    }
    N = length(bedg)
    for (i in bedg) {
        print(i)
        tmp = tempfile()
        com = paste("tabix ", i, " ", chr, ":", a, "-", b, " > ",
            tmp, sep = "")
        system(com)
        x = read.table(tmp)
	unlink(tmp)
        x[1, 2] = a
        x[nrow(x), 3] = b + 1
        z =  rep(x[[4]], x[, 3] - x[, 2])
    }
    genplot(a:b, z, type = "l", chr = chr, xlim = xlim,
        gname=gname, col = NA, ylim = c(0, max(z)), biotype = biotype,
        MaxNT = MaxNT, BED = BED)
    if (LINE) {
        lines(a:b, z, col = col[i], lwd = 1, type = "s")
        abline(0, 0)
    }
    else {
        polygon(c(a, a:b, b), c(0, z, 0),  col = col[1], border = col[1], lwd = 2)
    }
}
