export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/htslib/lib

rm src/genplot.so
rm src/main.o

CFLAGS="-I/path/to/htslib\ -I/usr/include" 
LDFLAGS="-L/usr/lib\ -L/usr/local/lib\ -L/path/to/htslib/lib"

LIBR="-lhts\ -lm"

echo $LD_LIBRARY_PATH
echo $CFLAGS
echo $LDFLAGS

MAKEFLAGS="CFLAGS=$CFLAGS LDFLAGS=$LDFLAGS LIBR=$LIBR" R CMD INSTALL ./ --no-multiarch
