export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/htslib/lib

rm src/genplot.so
rm src/main.o

CFLAGS="-I/path/to/htslib\ -I/path/to/gsl\ -I/usr/include" #/Users/nk5/Applications/gsl-1.16/gsl\ -I/Users/nk5/Applications/gsl-1.16"
LDFLAGS="-L/usr/lib\ -L/usr/local/lib\ -L/path/to/htslib/lib"

LIBR="-lhts\ -lgsl\ -llapack\ -lm"

echo $LD_LIBRARY_PATH
echo $CFLAGS
echo $LDFLAGS

MAKEFLAGS="CFLAGS=$CFLAGS LDFLAGS=$LDFLAGS LIBR=$LIBR" R CMD INSTALL ./ --no-multiarch
#MAKEFLAGS="CFLAGS=$CFLAGS LDFLAGS=$LDFLAGS LIBR=$LIBR" R CMD INSTALL master --no-multiarch
