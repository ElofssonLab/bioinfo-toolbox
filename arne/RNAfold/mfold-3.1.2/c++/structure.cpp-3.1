#if !defined(STRUCTURE_CPP)
#define STRUCTURE_CPP

#include "structure.h"



//***********************************code for Structures:


datatable::datatable()
{
int a,b,c,d,e,f,g,h;






for (a=0;a<=5;a++) {
	for (b=0;b<=5;b++) {
		for (c=0;c<=5;c++) {
			for (d=0;d<=5;d++) {
				for (e=0;e<=5;e++) {
					for (f=0;f<=5;f++) {
               	iloop11[a][b][c][d][e][f] = 0;
						for (g=0;g<=5;g++) {
			iloop21[a][b][c][d][e][f][g] = infinity;
         for (h=0;h<=5;h++) {
          	iloop22[a][b][c][d][e][f][g][h] = infinity;
         }
						}
					}
				}
			}
		}
	}
}








};



structure::structure()
{
	int i;
	for (i=1;i<=maxstructures;i++) {
		energy[i]=0;
   allocated = false;
	}
	nnopair=0;
	npair=0;
	ndbl=0;
   intermolecular = false;
   ngu = 0;
   templated = false;
}

structure::~structure()
{
	int i;

   if (allocated) {
		delete[] numseq;
   	for (i=0;i<=maxstructures;i++) {
    		delete[] basepr[i];
   	}

   	delete[] basepr;
      delete[] hnumber;
      delete[] nucs;
   }
   if (templated) {
		for (i=0;i<=numofbases;i++) {
    		delete[] tem[i];
   	}

   	delete[] tem;
   }

}



void structure::allocate(int size)

{
	int i;
	//Size = size;//save the size of the array so that the destructor can
   				//deallocate the space
   numseq = new int [2*size+1];
   hnumber = new int [size+1];
   nucs = new char [size+2];
   basepr = new int *[maxstructures+1];
   for (i=0;i<=maxstructures;i++) {
    	basepr[i] = new int [size+1];
   }
   allocated = true;

}


//this allocates space in an array that is used for folding with phylogenetic data
void structure::allocatetem()

{
	int i;
	//Size = size;//save the size of the array so that the destructor can
   				//deallocate the space

   tem = new bool *[numofbases+1];
   for (i=0;i<=numofbases;i++) {
    	tem[i] = new bool [i+1];
   }
   templated = true;

}

#endif
