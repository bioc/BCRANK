#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>

#ifndef MATCH_CONSENSUS
#define MATCH_CONSENSUS

void match_consensus(char** sequences, char** motif, char** RCmotif, int* motifLength, int* nrSeqs, int* seqLengths, int* maxOccurences, int* total_match);
	
#endif
