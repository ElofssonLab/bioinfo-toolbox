
#if !defined(STRUCTURE_H)
#define STRUCTURE_H



#include "defines.h"

struct datatable //this structure contains all the info read from thermodynamic
							//data files
{
 int poppen [5],maxpen,eparam[11],dangle[6][6][6][3],inter[31],bulge[31],
		hairpin[31],stack[6][6][6][6],tstkh[6][6][6][6],tstki[6][6][6][6],
		tloop[maxtloop+1][2],numoftloops,iloop22[6][6][6][6][6][6][6][6],
		iloop21[6][6][6][6][6][6][6],iloop11[6][6][6][6][6][6],
      coax[6][6][6][6],tstackcoax[6][6][6][6],coaxstack[6][6][6][6],
      tstack[6][6][6][6],tstkm[6][6][6][6],auend,gubonus,cint,cslope,c3,
      efn2a,efn2b,efn2c,triloop[maxtloop+1][2],numoftriloops,init,gail;

float prelog;
datatable();

};



//////////////////////////////////////////////////////////////////////
struct structure //this structure contains all the info for a structure
{
int numofbases,numofstructures;
int pair[maxforce][2],npair;
int *numseq,*hnumber;
int **basepr;
int ndbl, dbl[maxforce];
int energy[maxstructures+1],inter[3];
int nnopair,nopair[maxforce];
int ngu,gu[maxgu];
char ctlabel[maxstructures+1][ctheaderlength],*nucs;
bool intermolecular,allocated,templated;
bool **tem;
//int **fce;//[maxbases+1][2*maxbases]
structure();
~structure();
void allocate(int size = maxbases);
void allocatetem();
/*structure is set up to hold many possible structures of the same sequence
	numofbases = number of bases in sequence
	numofstructures = number of alternative structures of the sequence
				that is held by structure
	numseq[i] = a numeric that stands for the base in the ith position
				of the sequence,
			A = 1
			C = 2
			G = 3
			U = 4
	basepr[i][j] = base to which the jth base is paired in the ith structure
	force[i] = any information about the way the ith pair is treated when
				folding; eg: forced single, etc
	energy[i] = 10 * the Gibb's free energy of the ith structure, this
				is done so that integer math can be performed
				and the nearest tenth of a kcal/mol can be
				followed
	ctlabel = a string of information for each of the structures
	fce = an array that can keep track of how each i,j pair is being
				forced
   hnumber array stores the historical numbering of a sequence
   nucs is a character array to store the sequence information --
   	this will allow the program to keep T and U from getting confused
*/

};

#endif
