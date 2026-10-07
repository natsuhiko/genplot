#include "loadGTF.h"

/*
 * Safer R/.Call interface for tabix-backed GTF/BED readers.
 *
 * Main fixes relative to the original code:
 *   - Protect coerced R arguments.
 *   - Check hts_open(), tbx_index_load(), and tbx_itr_querys().
 *   - Handle regions containing no gene records.
 *   - Never call mkChar(NULL); missing GTF attributes become NA.
 *   - Keep output objects protected until the result pairlist is complete.
 *   - Do not free pointers that point inside kstring_t::s.
 *   - Bound sscanf() string fields to avoid buffer overflow.
 *   - Release htslib/kstring resources before returning.
 */

static SEXP mkCharOrNA(const char *s)
{
    return s == NULL ? NA_STRING : mkChar(s);
}

int parseAttrib(char *sattr, GTFATTRIB *attr)
{
    int i;
    int flag = 0;
    char *key = NULL;
    char *val = NULL;
    int n;

    if (attr == NULL) return -1;

    attr->gene_id = NULL;
    attr->transcript_id = NULL;
    attr->gene_type = NULL;
    attr->gene_status = NULL;
    attr->gene_name = NULL;
    attr->transcript_type = NULL;
    attr->transcript_status = NULL;
    attr->transcript_name = NULL;
    attr->exon_id = NULL;
    attr->exon_number = -1;
    attr->level = -1;

    if (sattr == NULL) return 0;

    n = (int)strlen(sattr);

    for (i = 0; i < n; i++) {
        if (flag == 0 && sattr[i] != ' ' && sattr[i] != '\t' &&
            sattr[i] != '"' && sattr[i] != ';') {
            key = sattr + i;
            flag = 1;
        }

        if (flag == 2 && sattr[i] != ' ' && sattr[i] != '\t' && sattr[i] != '"') {
            val = sattr + i;
            flag = 3;
        }

        if (flag == 1 && (sattr[i] == ' ' || sattr[i] == '\t')) {
            sattr[i] = '\0';
            flag = 2;
        }

        if (flag == 3 && (sattr[i] == '"' || sattr[i] == ';')) {
            sattr[i] = '\0';

            if (key != NULL && val != NULL) {
                if (strcmp(key, "gene_id") == 0) {
                    attr->gene_id = val;
                } else if (strcmp(key, "transcript_id") == 0) {
                    attr->transcript_id = val;
                } else if (strcmp(key, "gene_type") == 0) {
                    attr->gene_type = val;
                } else if (strcmp(key, "gene_status") == 0) {
                    attr->gene_status = val;
                } else if (strcmp(key, "gene_name") == 0) {
                    attr->gene_name = val;
                } else if (strcmp(key, "gene_biotype") == 0) {
                    attr->transcript_type = val;
                } else if (strcmp(key, "transcript_type") == 0) {
                    attr->transcript_type = val;
                } else if (strcmp(key, "transcript_status") == 0) {
                    attr->transcript_status = val;
                } else if (strcmp(key, "transcript_name") == 0) {
                    attr->transcript_name = val;
                } else if (strcmp(key, "exon_number") == 0) {
                    attr->exon_number = atoi(val);
                } else if (strcmp(key, "exon_id") == 0) {
                    attr->exon_id = val;
                } else if (strcmp(key, "level") == 0) {
                    attr->level = atoi(val);
                }
            }

            key = NULL;
            val = NULL;
            flag = 0;
        }
    }

    return 0;
}


SEXP loadGTF(SEXP Rfname, SEXP Rreg)
{
    SEXP fnameS = PROTECT(coerceVector(Rfname, STRSXP));
    SEXP regS   = PROTECT(coerceVector(Rreg, STRSXP));

    if (XLENGTH(fnameS) < 1 || XLENGTH(regS) < 1 ||
        STRING_ELT(fnameS, 0) == NA_STRING || STRING_ELT(regS, 0) == NA_STRING) {
        UNPROTECT(2);
        return R_NilValue;
    }

    const char *fname = CHAR(STRING_ELT(fnameS, 0));
    const char *reg   = CHAR(STRING_ELT(regS, 0));

    verbose = 0;

    char regchr[1000];
    int regstart = 0, regend = 0;
    if (sscanf(reg, "%999[^:]:%d-%d", regchr, &regstart, &regend) != 3) {
        fprintf(stderr, "Invalid region: %s\n", reg);
        UNPROTECT(2);
        return R_NilValue;
    }

    htsFile *fp = hts_open(fname, "r");
    if (fp == NULL) {
        fprintf(stderr, "Could not read tabixed file %s\n", fname);
        UNPROTECT(2);
        return R_NilValue;
    }

    size_t fnidx_len = strlen(fname) + 5;
    char *fnidx = (char *)calloc(fnidx_len, 1);
    if (fnidx == NULL) {
        hts_close(fp);
        UNPROTECT(2);
        return R_NilValue;
    }
    snprintf(fnidx, fnidx_len, "%s.tbi", fname);

    tbx_t *tbx = tbx_index_load(fnidx);
    if (tbx == NULL) {
        fprintf(stderr, "Could not load .tbi index of %s\n", fnidx);
        free(fnidx);
        hts_close(fp);
        UNPROTECT(2);
        return R_NilValue;
    }

    kstring_t str = {0, 0, NULL};
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    if (itr == NULL) {
        fprintf(stderr, "Could not create tabix iterator for region %s\n", reg);
        tbx_destroy(tbx);
        free(fnidx);
        hts_close(fp);
        UNPROTECT(2);
        return R_NilValue;
    }

    int gstart = 0;
    int gend = 0;
    int ngene = 0;

    char chr[1000];
    char source[1000];
    char ftype[1000];
    char score[1000];
    char strand[1000];
    char phase[1000];
    int fstart = 0, fend = 0, nchar = 0;

    while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
        int nr = sscanf(str.s,
                        "%999[^\t]\t%999[^\t]\t%999[^\t]\t%d\t%d\t%999[^\t]\t%999[^\t]\t%999[^\t]\t%n",
                        chr, source, ftype, &fstart, &fend,
                        score, strand, phase, &nchar);
        if (nr != 8) continue;

        if (strcmp(ftype, "gene") == 0) {
            if (ngene == 0) {
                gstart = fstart;
                gend = fend;
            } else {
                if (fstart < gstart) gstart = fstart;
                if (fend > gend) gend = fend;
            }
            ngene++;
        }
    }
    tbx_itr_destroy(itr);
    itr = NULL;

    if (ngene == 0) {
        if (verbose > 0) fprintf(stderr, "No gene overlaps region %s\n", reg);
        free(str.s);
        tbx_destroy(tbx);
        free(fnidx);
        hts_close(fp);
        UNPROTECT(2);
        return R_NilValue;
    }

    char reg2[1200];
    int nw = snprintf(reg2, sizeof(reg2), "%s:%d-%d", regchr, gstart, gend);
    if (nw < 0 || nw >= (int)sizeof(reg2)) {
        fprintf(stderr, "Expanded region is too long\n");
        free(str.s);
        tbx_destroy(tbx);
        free(fnidx);
        hts_close(fp);
        UNPROTECT(2);
        return R_NilValue;
    }

    if (verbose > 0) fprintf(stderr, "Expanded region: %s\n", reg2);

    /* Count exon/CDS/UTR records. */
    itr = tbx_itr_querys(tbx, reg2);
    if (itr == NULL) {
        fprintf(stderr, "Could not create tabix iterator for expanded region %s\n", reg2);
        free(str.s);
        tbx_destroy(tbx);
        free(fnidx);
        hts_close(fp);
        UNPROTECT(2);
        return R_NilValue;
    }

    int nfeature = 0;
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
        int nr = sscanf(str.s,
                        "%999[^\t]\t%999[^\t]\t%999[^\t]\t%d\t%d\t%999[^\t]\t%999[^\t]\t%999[^\t]\t%n",
                        chr, source, ftype, &fstart, &fend,
                        score, strand, phase, &nchar);
        if (nr != 8) continue;

        if (strcmp(ftype, "exon") == 0 || strcmp(ftype, "CDS") == 0 || strcmp(ftype, "UTR") == 0)
            nfeature++;
    }
    tbx_itr_destroy(itr);
    itr = NULL;

    if (verbose > 0) fprintf(stderr, "N of features = %d\n", nfeature);

    if (nfeature == 0) {
        free(str.s);
        tbx_destroy(tbx);
        free(fnidx);
        hts_close(fp);
        UNPROTECT(2);
        return R_NilValue;
    }

    SEXP Rsources = PROTECT(allocVector(STRSXP, nfeature));
    SEXP Rftypes  = PROTECT(allocVector(STRSXP, nfeature));
    SEXP Rfstarts = PROTECT(allocVector(INTSXP, nfeature));
    SEXP Rfends   = PROTECT(allocVector(INTSXP, nfeature));
    SEXP Rstrands = PROTECT(allocVector(INTSXP, nfeature));
    SEXP Rgids    = PROTECT(allocVector(STRSXP, nfeature));
    SEXP Rgnames  = PROTECT(allocVector(STRSXP, nfeature));
    SEXP Rtids    = PROTECT(allocVector(STRSXP, nfeature));
    SEXP Rbtypes  = PROTECT(allocVector(STRSXP, nfeature));

    itr = tbx_itr_querys(tbx, reg2);
    if (itr == NULL) {
        fprintf(stderr, "Could not create tabix iterator for expanded region %s\n", reg2);
        free(str.s);
        tbx_destroy(tbx);
        free(fnidx);
        hts_close(fp);
        UNPROTECT(11); /* 9 outputs + 2 inputs */
        return R_NilValue;
    }

    int l = 0;
    GTFATTRIB attr;

    while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
        int nr = sscanf(str.s,
                        "%999[^\t]\t%999[^\t]\t%999[^\t]\t%d\t%d\t%999[^\t]\t%999[^\t]\t%999[^\t]\t%n",
                        chr, source, ftype, &fstart, &fend,
                        score, strand, phase, &nchar);
        if (nr != 8) continue;

        if (strcmp(ftype, "exon") != 0 && strcmp(ftype, "CDS") != 0 && strcmp(ftype, "UTR") != 0)
            continue;

        if (l >= nfeature) {
            fprintf(stderr, "Internal error: more features encountered than counted\n");
            break;
        }

        char *attrib = str.s + nchar;  /* points inside str.s: never free separately */
        parseAttrib(attrib, &attr);

        if (verbose > 0) {
            fprintf(stderr, "%s %s %d %d %s | gene_id=%s gene_name=%s transcript_id=%s transcript_type=%s\n",
                    source, ftype, fstart, fend, strand,
                    attr.gene_id ? attr.gene_id : "NULL",
                    attr.gene_name ? attr.gene_name : "NULL",
                    attr.transcript_id ? attr.transcript_id : "NULL",
                    attr.transcript_type ? attr.transcript_type : "NULL");
        }

        SET_STRING_ELT(Rsources, l, mkChar(source));
        SET_STRING_ELT(Rftypes,  l, mkChar(ftype));
        INTEGER(Rfstarts)[l] = fstart;
        INTEGER(Rfends)[l]   = fend;
        INTEGER(Rstrands)[l] = strcmp(strand, "+") == 0 ? 0 : 1;
        SET_STRING_ELT(Rgids,   l, mkCharOrNA(attr.gene_id));
        SET_STRING_ELT(Rgnames, l, mkCharOrNA(attr.gene_name));
        SET_STRING_ELT(Rtids,   l, mkCharOrNA(attr.transcript_id));
        SET_STRING_ELT(Rbtypes, l, mkCharOrNA(attr.transcript_type));
        l++;
    }

    tbx_itr_destroy(itr);
    itr = NULL;
    free(str.s);
    tbx_destroy(tbx);
    free(fnidx);
    if (hts_close(fp) != 0)
        fprintf(stderr, "hts_close returned non-zero status: %s\n", fname);

    /* Preserve the original return type/order: a pairlist of 9 objects. */
    SEXP ans = PROTECT(allocList(9));
    SEXP p = ans;
    SETCAR(p, Rgids);     p = CDR(p);
    SETCAR(p, Rtids);     p = CDR(p);
    SETCAR(p, Rfstarts);  p = CDR(p);
    SETCAR(p, Rfends);    p = CDR(p);
    SETCAR(p, Rstrands);  p = CDR(p);
    SETCAR(p, Rsources);  p = CDR(p);
    SETCAR(p, Rftypes);   p = CDR(p);
    SETCAR(p, Rgnames);   p = CDR(p);
    SETCAR(p, Rbtypes);

    UNPROTECT(12); /* 2 inputs + 9 output vectors + ans */
    return ans;
}


SEXP loadBed(SEXP Rfname, SEXP Rreg)
{
    SEXP fnameS = PROTECT(coerceVector(Rfname, STRSXP));
    SEXP regS   = PROTECT(coerceVector(Rreg, STRSXP));

    if (XLENGTH(fnameS) < 1 || XLENGTH(regS) < 1 ||
        STRING_ELT(fnameS, 0) == NA_STRING || STRING_ELT(regS, 0) == NA_STRING) {
        UNPROTECT(2);
        return R_NilValue;
    }

    const char *fname = CHAR(STRING_ELT(fnameS, 0));
    const char *reg   = CHAR(STRING_ELT(regS, 0));
    verbose = 0;

    htsFile *fp = hts_open(fname, "r");
    if (fp == NULL) {
        fprintf(stderr, "Could not read tabixed file %s\n", fname);
        UNPROTECT(2);
        return R_NilValue;
    }

    size_t fnidx_len = strlen(fname) + 5;
    char *fnidx = (char *)calloc(fnidx_len, 1);
    if (fnidx == NULL) {
        hts_close(fp);
        UNPROTECT(2);
        return R_NilValue;
    }
    snprintf(fnidx, fnidx_len, "%s.tbi", fname);

    tbx_t *tbx = tbx_index_load(fnidx);
    if (tbx == NULL) {
        fprintf(stderr, "Could not load .tbi index of %s\n", fnidx);
        free(fnidx);
        hts_close(fp);
        UNPROTECT(2);
        return R_NilValue;
    }

    kstring_t str = {0, 0, NULL};
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    if (itr == NULL) {
        fprintf(stderr, "Could not create tabix iterator for region %s\n", reg);
        tbx_destroy(tbx);
        free(fnidx);
        hts_close(fp);
        UNPROTECT(2);
        return R_NilValue;
    }

    char chr[1000];
    int fstart = 0, fend = 0, nchar = 0;
    int nfeature = 0;
    int ncol = 0;

    while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
        int nr = sscanf(str.s, "%999[^\t]\t%d\t%d%n", chr, &fstart, &fend, &nchar);
        if (nr != 3) continue;

        if (nfeature == 0) {
            char *attrib = str.s + nchar;
            if (*attrib == '\t') {
                char *q;
                for (q = attrib; *q != '\0'; q++)
                    if (*q == '\t') ncol++;
            }
        }
        nfeature++;
    }
    tbx_itr_destroy(itr);
    itr = NULL;

    if (verbose > 0) fprintf(stderr, "N of features = %d\n", nfeature);

    if (nfeature == 0) {
        free(str.s);
        tbx_destroy(tbx);
        free(fnidx);
        hts_close(fp);
        UNPROTECT(2);
        return R_NilValue;
    }

    SEXP Rfstarts = PROTECT(allocVector(INTSXP, nfeature));
    SEXP Rfends   = PROTECT(allocVector(INTSXP, nfeature));
    SEXP Raddcol  = R_NilValue;
    if (ncol > 0)
        Raddcol = PROTECT(allocVector(STRSXP, nfeature));

    itr = tbx_itr_querys(tbx, reg);
    if (itr == NULL) {
        fprintf(stderr, "Could not create tabix iterator for region %s\n", reg);
        free(str.s);
        tbx_destroy(tbx);
        free(fnidx);
        hts_close(fp);
        UNPROTECT(ncol > 0 ? 5 : 4);
        return R_NilValue;
    }

    int l = 0;
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
        int nr = sscanf(str.s, "%999[^\t]\t%d\t%d%n", chr, &fstart, &fend, &nchar);
        if (nr != 3) continue;
        if (l >= nfeature) break;

        INTEGER(Rfstarts)[l] = fstart;
        INTEGER(Rfends)[l]   = fend;

        if (ncol > 0) {
            char *attrib = str.s + nchar;
            if (*attrib == '\t') attrib++;
            SET_STRING_ELT(Raddcol, l, mkChar(attrib));
        }
        l++;
    }

    tbx_itr_destroy(itr);
    free(str.s);
    tbx_destroy(tbx);
    free(fnidx);
    if (hts_close(fp) != 0)
        fprintf(stderr, "hts_close returned non-zero status: %s\n", fname);

    int nout = ncol > 0 ? 3 : 2;
    SEXP ans = PROTECT(allocList(nout));
    SEXP p = ans;
    SETCAR(p, Rfstarts); p = CDR(p);
    SETCAR(p, Rfends);
    if (ncol > 0) {
        p = CDR(p);
        SETCAR(p, Raddcol);
    }

    UNPROTECT(ncol > 0 ? 6 : 5); /* inputs + output vectors + ans */
    return ans;
}


SEXP tabix2charmat(SEXP Rfname, SEXP Rreg)
{
    SEXP fnameS = PROTECT(coerceVector(Rfname, STRSXP));
    SEXP regS   = PROTECT(coerceVector(Rreg, STRSXP));

    if (XLENGTH(fnameS) < 1 || XLENGTH(regS) < 1 ||
        STRING_ELT(fnameS, 0) == NA_STRING || STRING_ELT(regS, 0) == NA_STRING) {
        UNPROTECT(2);
        return R_NilValue;
    }

    const char *fname = CHAR(STRING_ELT(fnameS, 0));
    const char *reg   = CHAR(STRING_ELT(regS, 0));
    verbose = 0;

    htsFile *fp = hts_open(fname, "r");
    if (fp == NULL) {
        fprintf(stderr, "Could not read tabixed file %s\n", fname);
        UNPROTECT(2);
        return R_NilValue;
    }

    size_t fnidx_len = strlen(fname) + 5;
    char *fnidx = (char *)calloc(fnidx_len, 1);
    if (fnidx == NULL) {
        hts_close(fp);
        UNPROTECT(2);
        return R_NilValue;
    }
    snprintf(fnidx, fnidx_len, "%s.tbi", fname);

    tbx_t *tbx = tbx_index_load(fnidx);
    if (tbx == NULL) {
        fprintf(stderr, "Could not load .tbi index of %s\n", fnidx);
        free(fnidx);
        hts_close(fp);
        UNPROTECT(2);
        return R_NilValue;
    }

    kstring_t str = {0, 0, NULL};
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    if (itr == NULL) {
        fprintf(stderr, "Could not create tabix iterator for region %s\n", reg);
        tbx_destroy(tbx);
        free(fnidx);
        hts_close(fp);
        UNPROTECT(2);
        return R_NilValue;
    }

    int nfeature = 0;
    int ncol = 0;

    while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
        if (nfeature == 0) {
            char *q;
            for (q = str.s; *q != '\0'; q++)
                if (*q == '\t') ncol++;
        }
        nfeature++;
    }
    tbx_itr_destroy(itr);
    itr = NULL;

    if (verbose > 0) fprintf(stderr, "Dim = %d x %d\n", nfeature, ncol + 1);

    if (nfeature == 0) {
        free(str.s);
        tbx_destroy(tbx);
        free(fnidx);
        hts_close(fp);
        UNPROTECT(2);
        return R_NilValue;
    }

    int total_cols = ncol + 1;
    SEXP Raddcol = PROTECT(allocVector(STRSXP, (R_xlen_t)nfeature * total_cols));
    SEXP Rdim    = PROTECT(allocVector(INTSXP, 2));
    INTEGER(Rdim)[0] = nfeature;
    INTEGER(Rdim)[1] = total_cols;

    itr = tbx_itr_querys(tbx, reg);
    if (itr == NULL) {
        fprintf(stderr, "Could not create tabix iterator for region %s\n", reg);
        free(str.s);
        tbx_destroy(tbx);
        free(fnidx);
        hts_close(fp);
        UNPROTECT(4); /* 2 inputs + Raddcol + Rdim */
        return R_NilValue;
    }

    R_xlen_t l = 0;
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
        char *field = str.s;
        int c;

        for (c = 0; c < total_cols; c++) {
            char *tab = strchr(field, '\t');

            if (tab != NULL && c < total_cols - 1) {
                char saved = *tab;
                *tab = '\0';
                SET_STRING_ELT(Raddcol, l++, mkChar(field));
                *tab = saved;
                field = tab + 1;
            } else {
                SET_STRING_ELT(Raddcol, l++, mkChar(field));
                field += strlen(field);
            }
        }
    }

    tbx_itr_destroy(itr);
    free(str.s);
    tbx_destroy(tbx);
    free(fnidx);
    if (hts_close(fp) != 0)
        fprintf(stderr, "hts_close returned non-zero status: %s\n", fname);

    /* Preserve the original return type/order: pairlist(Rdim, Raddcol). */
    SEXP ans = PROTECT(allocList(2));
    SEXP p = ans;
    SETCAR(p, Rdim); p = CDR(p);
    SETCAR(p, Raddcol);

    UNPROTECT(5); /* 2 inputs + Raddcol + Rdim + ans */
    return ans;
}

