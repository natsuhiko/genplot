#include "loadGTF.h"




/*
int main0(int argc, char** argv){
    char** sources;
    char** ftypes;
    int* fstarts;
    int* fends;
    int* strands;
    char** gids;
    char** gnames;
    char** tids;
    char** btypes;
    loadGTF(&(argv[1]), &(argv[2]), sources, ftypes, fstarts, fends, strands, gids, gnames, tids, btypes);
    return 0;
}
*/


int parseAttrib(char* sattr, GTFATTRIB* attr){
    int i;
    int flag=0;
    int flag_double_quote=0;
    char* key;
    char* val;
    int n = strlen(sattr);
    attr->gene_id = attr->transcript_id = attr->gene_type = attr->gene_status = attr->gene_name = attr->transcript_type = attr->transcript_status = attr->transcript_name = attr->exon_id = NULL;
    attr->exon_number = attr->level = -1;
    for(i=0; i<n; i++){
        //fprintf(stderr, "%d", flag);
        if(flag==0 && sattr[i]!=' ' && sattr[i]!='\"' && sattr[i]!=';'){
            key  = sattr+i;
            flag = 1;
        }
        if(flag==2 && sattr[i]!=' ' && sattr[i]!='\"'){
            if(sattr[i]!='\"'){flag_double_quote=1;}
            val  = sattr+i;
            flag = 3;
        }
        if(flag==1 && sattr[i]==' '){
            sattr[i] = '\0';
            //fprintf(stderr, "%s\n", key);
            flag=2;
        }
        if(flag==3 && (sattr[i]=='\"' || sattr[i]==';')){
            sattr[i] = '\0';
            if(strcmp(key, "gene_id")==0){
                attr->gene_id = val;
            }else if(strcmp(key, "transcript_id")==0 ){
                attr->transcript_id = val;
            }else if(strcmp(key, "gene_type")==0 ){
                attr->gene_type = val;
            }else if(strcmp(key, "gene_status")==0 ){
                attr->gene_status = val;
            }else if(strcmp(key, "gene_name")==0 ){
                attr->gene_name = val;
            }else if(strcmp(key, "transcript_type")==0 ){
                attr->transcript_type = val;
            }else if(strcmp(key, "transcript_status")==0 ){
                attr->transcript_status = val;
            }else if(strcmp(key, "transcript_name")==0 ){
                attr->transcript_name = val;
            }else if(strcmp(key, "exon_number")==0 ){
                attr->exon_number = atoi(val);
            }else if(strcmp(key, "exon_id")==0 ){
                attr->exon_id = val;
            }else if(strcmp(key, "level")==0 ){
                attr->level = atoi(val);
            }
            //fprintf(stderr, "%s=%s\n", key, val);
            flag=0;
        }
    }
    return 0;
}

// vcf : tabixed fragment file
// hid is shifted for sample i
SEXP loadGTF(SEXP Rfname, SEXP Rreg){//, SEXP Rsources, SEXP Rftypes, SEXP Rfstarts, SEXP Rfends, SEXP Rstrands, SEXP Rgids, SEXP Rgnames, SEXP Rtids, SEXP Rbtypes){
    
    Rfname = coerceVector(Rfname, STRSXP);
    Rreg   = coerceVector(Rreg,   STRSXP);
    
    const char* fname; fname = CHAR(STRING_ELT(Rfname, 0));
    const char* reg;   reg   = CHAR(STRING_ELT(Rreg, 0));
    
    verbose=0;
    
    char* regchr; regchr=(char*)calloc(1000, sizeof(char));
    int regstart, regend;
    sscanf(reg, "%[^:]:%d-%d", regchr, &regstart, &regend);
    
    //fprintf(stderr, "%s %d %d", regchr, regstart, regend);
    
    //return 0;
    
    int i;
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) fprintf(stderr, "Could not read tabixed file %s\n", fname);
    //enum htsExactFormat format = hts_get_format(fp)->format;
    
    char *fnidx = calloc(strlen(fname) + 5, 1);
    strcat(strcpy(fnidx, fname), ".tbi");
    
    //regidx_t *reg_idx = NULL;
    
    tbx_t *tbx = tbx_index_load(fnidx);
    if ( !tbx ) fprintf(stderr, "Could not load .tbi index of %s\n", fnidx);
    
    kstring_t str = {0,0,0};

    //int nseq;
    //const char **seq = NULL;
    //if ( reg_idx ) seq = tbx_seqnames(tbx, &nseq);
    
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    
    
    int gstart=270000000, gend=0;
    char* chr;    chr   =(char*)calloc(1000, sizeof(char));
    char* source; source=(char*)calloc(1000, sizeof(char));
    char* ftype;  ftype =(char*)calloc(1000, sizeof(char));
    int fstart;
    int fend;
    char* score;  score =(char*)calloc(1000, sizeof(char));
    char* strand; strand=(char*)calloc(1000, sizeof(char));
    char* phase;  phase =(char*)calloc(1000, sizeof(char));
    char* attrib; attrib=(char*)calloc(1000, sizeof(char));
    int nchar;
    int ngene=0;
    int nfeature=0;
    
    
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        //if ( reg_idx && !regidx_overlap(reg_idx,seq[itr->curr_tid],itr->curr_beg,itr->curr_end, NULL) ) continue;
        
        sscanf(str.s, "%[^\t]\t%[^\t]\t%[^\t]\t%d\t%d\t%[^\t]\t%[^\t]\t%[^\t]\t%n", chr, source, ftype, &fstart, &fend, score, strand, phase, &nchar);
        attrib = str.s+nchar;
        
        nfeature++;
        if(strcmp(ftype, "gene")==0){
            ngene++;
            if(fstart<gstart){gstart=fstart;}
            if(fend  >gend  ){gend  =fend;  }
            if(verbose>1)puts(str.s);
        }
    }
    tbx_itr_destroy(itr);
    
    char* reg2; reg2 = calloc(1000, sizeof(char));
    sprintf(reg2, "%s:%d-%d", regchr, gstart, gend);
    
    if(verbose>0){fprintf(stderr, "\nExpanded region: %s\n\n", reg2);}
    
    
    // #
    // #  counting exon/UTR/CDS in region specified by genes overlapping with initial region
    // #
    itr = tbx_itr_querys(tbx, reg2);
    nfeature=0;
    GTFATTRIB attr;
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        //if ( reg_idx && !regidx_overlap(reg_idx,seq[itr->curr_tid],itr->curr_beg,itr->curr_end, NULL) ) continue;
        
        sscanf(str.s, "%[^\t]\t%[^\t]\t%[^\t]\t%d\t%d\t%[^\t]\t%[^\t]\t%[^\t]\t%n", chr, source, ftype, &fstart, &fend, score, strand, phase, &nchar);
        attrib = str.s+nchar;
        if(strcmp(ftype, "exon")==0 || strcmp(ftype, "CDS")==0 || strcmp(ftype, "UTR")==0){
            nfeature++;
        }
    }
    tbx_itr_destroy(itr);
    
    
    if(verbose>0){fprintf(stderr, "\nN of features = %d\n\n", nfeature);}
    
    if(nfeature==0){return R_NilValue;}
    
    // #
    // #  load gtf in region specified by genes overlapping with initial region
    // #
    SEXP Rsources = PROTECT(allocVector(STRSXP, nfeature));
    SEXP Rftypes  = PROTECT(allocVector(STRSXP, nfeature));
    SEXP Rfstarts = PROTECT(allocVector(INTSXP, nfeature));
    SEXP Rfends   = PROTECT(allocVector(INTSXP, nfeature));
    SEXP Rstrands = PROTECT(allocVector(INTSXP, nfeature));
    SEXP Rgids    = PROTECT(allocVector(STRSXP, nfeature));
    SEXP Rgnames  = PROTECT(allocVector(STRSXP, nfeature));
    SEXP Rtids    = PROTECT(allocVector(STRSXP, nfeature));
    SEXP Rbtypes  = PROTECT(allocVector(STRSXP, nfeature));
    
    //return R_NilValue;
    /*char** sources = (char**)calloc(nfeature, sizeof(char*));
    char** ftypes  = (char**)calloc(nfeature, sizeof(char*));
    int*   fstarts = (int*)calloc(nfeature, sizeof(int));
    int*   fends   = (int*)calloc(nfeature, sizeof(int));
    int*   strands = (int*)calloc(nfeature, sizeof(int));
    char** gids    = (char**)calloc(nfeature, sizeof(char*));
    char** gnames  = (char**)calloc(nfeature, sizeof(char*));
    char** tids    = (char**)calloc(nfeature, sizeof(char*));
    char** btypes  = (char**)calloc(nfeature, sizeof(char*));*/
    itr = tbx_itr_querys(tbx, reg2);
    
    int l=0;
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        //if ( reg_idx && !regidx_overlap(reg_idx,seq[itr->curr_tid],itr->curr_beg,itr->curr_end, NULL) ) continue;
        
        sscanf(str.s, "%[^\t]\t%[^\t]\t%[^\t]\t%d\t%d\t%[^\t]\t%[^\t]\t%[^\t]\t%n", chr, source, ftype, &fstart, &fend, score, strand, phase, &nchar);
        attrib = str.s+nchar;
        if(strcmp(ftype, "exon")==0 || strcmp(ftype, "CDS")==0 || strcmp(ftype, "UTR")==0){
            //puts(str.s);
            parseAttrib(attrib, &attr);
            SET_STRING_ELT(Rsources, l, mkChar(source));
            SET_STRING_ELT(Rftypes,  l, mkChar(ftype));
            INTEGER(Rfstarts)[l] = fstart;
            INTEGER(Rfends)[l]   = fend;
            INTEGER(Rstrands)[l] = strcmp(strand, "+")==0 ? 0 : 1;
            SET_STRING_ELT(Rgids,    l, mkChar(attr.gene_id));
            SET_STRING_ELT(Rgnames,  l, mkChar(attr.gene_name));
            SET_STRING_ELT(Rtids,    l, mkChar(attr.transcript_id));
            SET_STRING_ELT(Rbtypes,  l, mkChar(attr.transcript_type));
            
            /*sources[l] = (char*)calloc(strlen(source), sizeof(char)); strcpy(sources[l], source);
            ftypes[l]  = (char*)calloc(strlen(ftype), sizeof(char)); strcpy(ftypes[l], ftype);
            
            fstarts[l] = fstart;
            fends[l]   = fend;
            strands[l] = strcmp(strand, "+")==0 ? 0 : 1;
            
            gids[l]    = (char*)calloc(strlen(attr.gene_id), sizeof(char)); strcpy(gids[l], attr.gene_id);
            gnames[l]  = (char*)calloc(strlen(attr.gene_name), sizeof(char)); strcpy(gnames[l], attr.gene_name);
            tids[l]    = (char*)calloc(strlen(attr.transcript_id), sizeof(char)); strcpy(tids[l], attr.transcript_id);
            btypes[l]  = (char*)calloc(strlen(attr.transcript_type), sizeof(char)); strcpy(btypes[l], attr.transcript_type);*/
            l++;
        }
    }
    
    UNPROTECT(9);
    
    return CONS(Rgids, CONS(Rtids, CONS(Rfstarts, CONS(Rfends, CONS(Rstrands, CONS(Rsources, CONS(Rftypes, CONS(Rgnames, CONS(Rbtypes, R_NilValue)))))))));
    
    
    tbx_itr_destroy(itr);
    
    //free(seq);
    free(regchr);
    free(fnidx);
    free(reg2);
    
    free(chr);
    free(source);
    free(ftype);
    free(score);
    free(strand);
    free(phase);
    free(attrib);
    
    free(str.s);
    tbx_destroy(tbx);
    
    //if ( reg_idx ) regidx_destroy(reg_idx);
    if ( hts_close(fp) ) fprintf(stderr, "hts_close returned non-zero status: %s\n", fname);
    
    
    UNPROTECT(1);
    
    return Rsources;
}










SEXP loadBed(SEXP Rfname, SEXP Rreg){//, SEXP Rsources, SEXP Rftypes, SEXP Rfstarts, SEXP Rfends, SEXP Rstrands, SEXP Rgids, SEXP Rgnames, SEXP Rtids, SEXP Rbtypes){
    
    Rfname = coerceVector(Rfname, STRSXP);
    Rreg   = coerceVector(Rreg,   STRSXP);
    
    const char* fname; fname = CHAR(STRING_ELT(Rfname, 0));
    const char* reg;   reg   = CHAR(STRING_ELT(Rreg, 0));
    
    verbose=0;
    
    char* regchr; regchr=(char*)calloc(1000, sizeof(char));
    int regstart, regend;
    sscanf(reg, "%[^:]:%d-%d", regchr, &regstart, &regend);
    if(verbose>0){ fprintf(stderr, "Reg=%s\n", reg); }
    int i, k;
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) fprintf(stderr, "Could not read tabixed file %s\n", fname);
    
    char *fnidx = calloc(strlen(fname) + 5, 1);
    strcat(strcpy(fnidx, fname), ".tbi");
    
    tbx_t *tbx = tbx_index_load(fnidx);
    if ( !tbx ) fprintf(stderr, "Could not load .tbi index of %s\n", fnidx);
    
    kstring_t str = {0,0,0};
    
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    
    int gstart=270000000, gend=0;
    char* chr;    chr   =(char*)calloc(1000, sizeof(char));
    int fstart;
    int fend;
    char* attrib; //attrib=(char*)calloc(1000, sizeof(char));
    int nchar;
    int nfeature=0; // n of rows
    
    int ncol=0; // n of extra columns in the bed file
    //int* vartype;// 0: char;  1: color;  2: double;  3: int;
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%d%n", chr, &fstart, &fend, &nchar);
        attrib = str.s+nchar;
        // count ncol
        if(nfeature==0){
            if(attrib[0]=='\t'){for(k=0; k<strlen(attrib); k++){if(attrib[k]=='\t'){ncol++;}}}
            //if(ncol>0){vartype = (int*)calloc(ncol, sizeof(int));}
        }
        nfeature++;
    }
    tbx_itr_destroy(itr);
    
    if(verbose>0){fprintf(stderr, "\nN of features = %d\n\n", nfeature);}
    
    if(nfeature==0){return R_NilValue;}
    
    // #
    // #  load bed
    // #
    SEXP Rfstarts = PROTECT(allocVector(INTSXP, nfeature));
    SEXP Rfends   = PROTECT(allocVector(INTSXP, nfeature));
    SEXP Raddcol  = PROTECT(allocVector(STRSXP, nfeature));
    
    itr = tbx_itr_querys(tbx, reg);
    
    int l=0;
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%d%n", chr, &fstart, &fend, &nchar);
        
        INTEGER(Rfstarts)[l] = fstart;
        INTEGER(Rfends)[l]   = fend;
        if(ncol>0){
            attrib = str.s+nchar+1;
            SET_STRING_ELT(Raddcol, l, mkChar(attrib));
        }
        l++;
    }
    UNPROTECT(3);
    
    if(ncol>0){
        return CONS(Rfstarts, CONS(Rfends, CONS(Raddcol, R_NilValue)));
    }
    return CONS(Rfstarts, CONS(Rfends, R_NilValue));
}






























