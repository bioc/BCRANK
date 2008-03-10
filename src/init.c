#include <R.h>
#include <Rinternals.h>
#include "match_consensus.h"
#include <R_ext/Rdynload.h>

static R_NativePrimitiveArgType match_consensus_t[8] = {
	STRSXP,  /* 		char** sequences */
	STRSXP,  /* 		char** motif */
	STRSXP,  /* 		char** RCmotif */
	INTSXP,  /* 		int* motifLength */
	INTSXP,  /* 		int* nrSeq */
	INTSXP,  /* 		int* seqLengths */
	INTSXP,  /* 		int* maxOccurences */
	INTSXP   /* 		int* totalMatch */
};

static const R_CMethodDef CEntries[]  = {
		{"match_consensus", (DL_FUNC) &match_consensus, 8, match_consensus_t},
		{NULL, NULL, 0}
};


void R_init_BCRANK(DllInfo *dll)
{
		R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
		R_useDynamicSymbols(dll, FALSE);
}
