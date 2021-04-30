
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include "htslib/tbx.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/regidx.h"

#include <zlib.h>
#include <math.h>

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

int verbose;

typedef struct
{
    char* gene_id;          //	ENSGXXXXXXXXXXX *
    char* transcript_id;    //	ENSTXXXXXXXXXXX *
    char* gene_type;        //	list of biotypes
    char* gene_status;      //	{KNOWN, NOVEL, PUTATIVE}
    char* gene_name;        //	string;
    char* transcript_type;  //	list of biotypes
    char* transcript_status;//	{KNOWN, NOVEL, PUTATIVE}
    char* transcript_name;  //	string
    int exon_number;        //	indicates the biological position of the exon in the transcript
    char* exon_id;          //	ENSEXXXXXXXXXXX *
    int level;
}GTFATTRIB;


//SEXP loadGTF(char** Rfname, char** reg, char** sources, char** ftypes, SEXP Rfstarts, int* fends, int* strands, char** gids, char** gnames, char** tids, char** btypes);
SEXP loadGTF(SEXP Rfname, SEXP Rreg);//, SEXP Rsources, SEXP Rftypes, SEXP Rfstarts, SEXP Rfends, SEXP Rstrands, SEXP Rgids, SEXP Rgnames, SEXP Rtids, SEXP Rbtypes);
//int loadGTF(char* vcf, char* reg);
