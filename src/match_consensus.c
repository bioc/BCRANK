// Functions for matching IUPAC consensus to DNA sequences 

#include "match_consensus.h"

int IUPACmatch[256][256];

// Returns 1 if a nucleotide is matching an IUPAC symbol, 0 otherwise
int match_DNA_and_IUPAC(char DNA, char IUPAC){
  int match = 0;
  
  if(IUPAC == 'N'){
    return 1;
  }
  
  switch(DNA){
  case 'A':
    match = (IUPAC == 'A' || IUPAC == 'M' || IUPAC == 'W' || IUPAC == 'R' || IUPAC == 'D' || IUPAC =='H' || IUPAC =='V');
    break;
  case 'C':
    match = (IUPAC == 'C' || IUPAC == 'M' || IUPAC == 'S' || IUPAC == 'Y' || IUPAC == 'B' || IUPAC =='H' || IUPAC =='V');
    break;
  case 'G':
    match = (IUPAC == 'G' || IUPAC == 'R' || IUPAC == 'S' || IUPAC == 'K' || IUPAC == 'B' || IUPAC =='D' || IUPAC =='V');
    break;
  case 'T':
    match = (IUPAC == 'T' || IUPAC == 'W' || IUPAC == 'K' || IUPAC == 'Y' || IUPAC == 'B' || IUPAC =='D' || IUPAC =='H');
    break;
  }

  return match;

}

void initIUPACmatch(){
	int i,j;
	char DNA[4] = {'A','C','G','T'};
	char IUPAC[15] = {'A','C','G','T','R','Y','K','M','S','W','B','D','H','V','N'};

	for(i=0; i<4; i++){
		for(j=0; j<15; j++){
			IUPACmatch[(int) DNA[i]][(int) IUPAC[j]] = match_DNA_and_IUPAC(DNA[i],IUPAC[j]);
		}
	}
	
}


// Takes a consensus and its reverse complement and a sequence. Returns the total number of
// motif occcurences in the sequence.
int match_occurences(char* motif, char* RCmotif, int motifLength, int seqLength, int maxOccurences, char* line){
  int pos, matching, matchingRC,someMatching;
  int nr_match=0, i=0;
  char inputChar;

  for(i=0; i<=seqLength-motifLength; i++){

    matching = 1;
    matchingRC = 1;
	someMatching = 1;
	pos = 0;
	
    while((pos<motifLength) && someMatching){
    		
		inputChar = line[i+pos];
		
		if(matching){
			matching = IUPACmatch[(int) inputChar][(int) motif[pos]];
		}
		if(matchingRC){
			matchingRC = IUPACmatch[(int) inputChar][(int) RCmotif[pos]];
		}
		
		someMatching = (matching||matchingRC);
		
		pos++;
    }
    
    if(someMatching){
		nr_match++;
		if(nr_match>=maxOccurences){
			return(maxOccurences);
		}
    }
  }
  
  return nr_match;
}

// Matches all occurences of a motif and its reverse complement in a fasta file.
void match_consensus(char** sequences, char** motif, char** RCmotif, int* motifLength, int* nrSeqs, int* seqLengths, int *maxOccurences, int* total_match){
	int i, nr_match; 

	initIUPACmatch();
	
	for(i=0; i<*nrSeqs; i++){
		nr_match = match_occurences(*motif,*RCmotif,*motifLength,seqLengths[i],*maxOccurences, sequences[i]);
		total_match[i] = nr_match;
	}

}
