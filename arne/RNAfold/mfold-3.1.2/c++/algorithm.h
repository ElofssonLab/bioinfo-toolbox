
#if !defined(ALGORITHM_H)
#define ALGORITHM_H

#include "structure.h"
#include "defines.h"
#include "interface.h"


//***********************************Structures:

/////////////////////////////////////////////////////////////////////

struct stackstruct //this structure contains a stack of data, used by
						//	functions that analyze a structure piecewise
{
int stk[51][4],sp;
};

////////////////////////////////////////////////////////////////////

struct arraystruct //array used during the creation of a line out file (linout)
{
char array[7][amax][2];
};
//////////////////////////////////////////////////////////////////

void de_allocate (int **v,int i);//deallocates memory for a 2d array
void de_allocate (bool **v,int i);//alternative form of de_allocate
void de_allocate (int **v,int i);//alternative form of de_allocate

/////////////////////////////////////////////////////////////////////////
////dotarray encapsulates the array needed to store dot plot information

class dotarray{

		int **array;
      int store;

   public:

   	dotarray(int size);
      int &dot(int i, int j) {

      	return array[j][i];

      }
   	~dotarray();


};


////////////////////////////////////////////////////////////////////////
//arrayclass encapsulates the large 2-d arrays of w and v, used by the dynamic
//	algorithm
class arrayclass {
   int Size;

   public:
   	//int **dg;
      int k;
      int **dg;
      int infinite;
      //ofstream out;

      //the constructor allocates the space needed by the arrays
   	arrayclass(int size) {

      	//out.open("temp.out");

      	///k = infinity;
         infinite = infinity;

      	Size = size;
      	register int i,j;
      	dg = new int *[size+1];

   		for (i=0;i<=(size);i++)  {
   			dg[i] = new int [size+1];
   		}
         for (i=0;i<=size;i++) {
         	for (j=0;j<size+1;j++) {
            	dg[i][j] = infinity;
            }
         }

      }

      //the destructor deallocates the space used
      ~arrayclass() {

      	//out.close();
      	int i;
       	//de_allocate (dg,Size+2);
         for (i=0;i<=Size;i++) {
         	delete[] dg[i];
         }
         delete[] dg;
      }

      //f is an integer function that references the correct element of the array
   	int &f(int i, int j) {
      	//out << i<<"\t"<<j<<"\n";
      	if (i>j) {
         	return infinite;
         }
         else if (i>Size) return f(i-Size,j-Size);//dg[i-Size][j-Size];
         else return dg[i][j-i];
         //else return dg[i][j];
      }
};


//output is a union that is used by the dynamic algorithm to produce a save file
//	and used by the function opensav() to recall the data in a save file
union output {
  int i;
  char ch[4];
  float f;

};






//**********************************prototypes:

int opendat (char *loop2,char *stackf,char *tstackh,char *tstacki,
		char *tloop,char *miscloop, char *danglef, char *int22, char *int21,
      char *coax,char *tstackcoax,char *coaxstack,
      char *tstack,char *tstackm, char *triloop, char *int11, datatable* data);
		//gets thermodynamic data from data files
void tonum(char *base,structure *ct,int count); //converts base to a numeric
void getout (char *energyfile);//get name of file to output
										//	energy info
void openct(structure *ct,char *ctfile);//reads ct file
void efn2(datatable *data,structure *ct, int structnum = 0);//energy calculator
void push(stackstruct *stack,int a,int b,int c,int d);//push info onto the stack
void pull(stackstruct *stack,int *i,int *j,int *open,int *null,int *stz);
														//pull info from the stack
int erg1(int i,int j,int ip,int jp,structure *ct,datatable *data);
		//calculates energy of stacked base pairs
int erg2(int i,int j,int ip,int jp,structure *ct,datatable *data,int a,
	int b);
		//calculates energy of a bulge/internal loop
int erg3(int i,int j,structure *ct,datatable *data,int dbl);
		//calculates energy of a hairpin loop
int erg4(int i,int j,int ip,int jp,structure *ct,datatable *data,
	bool lfce);
		//calculates energy of a dangling base
int penalty(int i,int j,structure *ct,datatable *data);
	//calculates end of helix penalty
int penalty2(int i, int j, datatable *data);
void energyout(structure *ct,char *enrgyfile);
int openseq (structure *ct,char *seqfile);//inputs a sequence from file
															// seqfile
void ctout (structure *ct,char *ctoutfile);//outputs a ct file
void dynamic (structure *ct,datatable *data,int cntrl6,int cntrl8,int cntrl9,
	TProgressDialog* update=0,char* savfile = 0);
			//this is the dynamic folding algorithm of Zuker
         //cntrl6 = #tracebacks
         //cntrl8 = percent sort
         //cntrl9 = window
inline void swap(int *a,int *b);//Swap two variables
inline void swap(int *a,int *b);//swap two variables
int tonumi(char *base); //converts base to a numeric
void linout(structure *ct,char *file);//line printer output of structure
void digit (int row,int column,int pos,arraystruct* table);//used by linout
													//to write a digit into the line
                                       //printer output file
char *tobase (int i);//convert a numeric value for a base to the familiar
								//character
void sortstructures (structure *ct);//this routine resorts the structures
					//predicted by dynamic according to energies calculated by efn2

void errmsg(int err,int err1);//function for outputting info in case of an error
void update (int i);//function informs user of progress of fill algorithm

//the force... functions are used to initialize the arrays used to apply
//constraints during the folding process
void forcepair(int x,int y,structure *ct,arrayclass *v);
void forcesingle(int x,structure* ct,arrayclass *v);
void forcedbl(int dbl,structure* ct,int **w,bool *v);
void forceinter(int dbl,structure* ct,int **w);
void forceinterefn(int dbl,structure* ct,int **w);

//filter is used to choose structures to output after efn2
//	this can make the number of structures more reasonable for inspection
//	it takes a structure, ct, which also contains the final output,
//	percent sort, maximum number of structures, and window size
void filter(structure* ct, int percent, int max, int window);

//force is used to prepare arrays for the function dynamic, used during the
//	fill routines - it coordinates the force...() functions above
void force(structure *ct,arrayclass *v, int **fce, bool *lfce);
void opensav(char* filename, structure* ct, arrayclass* w, arrayclass* v,
	int *w3, int *w5,int *vmin,datatable *data);//opens a save file with information filled by
   									//fill algorithm
void traceback(structure *ct, datatable *data, arrayclass *v, arrayclass *w, int *w3, int *w5, int **fce,
	bool *lfce,int vmin, int cntrl6, int cntrl8, int cntrl9);//uses fill information to
   																			// make a suboptimal structures

//this function is used to calculate the values of all the dots in a dot plot
void dotefn2(structure *ct, datatable *data, arrayclass *v, arrayclass *w,
	int *w3, int *w5, int **fce, bool *lfce,int vmin,dotarray *dots,
   TProgressDialog* PD = 0);
void calcpnum(dotarray *dots, int *pnum, int increment, int numofbases,
	TProgressDialog *PD = 0);
void savefile(int i, ofstream* sav);//this function is used to make a save file
											//after the fill algorithm
int readfile(ifstream *read);//this function is used to read save files
void savedot(dotarray *dots,structure *ct, char *filename); //save dot plot info
void readdot(dotarray *dots, structure *ct, char *filename);//read a dot plot file
void dpalign(dotarray *dots1,dotarray *dots2,structure* ct1,structure *ct2,int *align);
int getbestdot(dotarray *dots1,dotarray *dots2, structure* ct1,
	structure *ct2, int i, int j);//return the best dot for base i
   //in dots1 and j in dots2
//dpalign will align two dot plots and store the info in the array align
void energydump (structure *ct, arrayclass *v, int n,char *filename);
//energydump will spit out the composite free energies for a traceback
void energydump (structure *ct, datatable *data,arrayclass *v, int n,char *filename);
//energydump2 will spit out the composite free energies for a traceback -- with
//the au penalty associated with the correct entity
#endif
