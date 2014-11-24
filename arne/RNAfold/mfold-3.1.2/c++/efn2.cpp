

#include "algorithm.cpp"


/*	Function getdat

	Function gets the names of data files to open

*/

void getdat(char *loop, char *stackf, char *tstackh, char *tstacki,
		char *tloop, char *miscloop, char *danglef, char *int22,
      char *int21,char *coax, char *tstackcoax,
      char *coaxstack, char *tstack, char *tstackm, char *triloop,
      char *int11,char *Path)

{

strcpy (loop,Path);
strcat (loop,"loop.dat");
strcpy (stackf,Path);
strcat (stackf,"stack.dat");
strcpy (tstackh,Path);
strcat (tstackh,"tstackh.dat");
strcpy (tstacki,Path);
strcat (tstacki,"tstacki.dat");
strcpy (tloop,Path);
strcat (tloop,"tloop.dat");
strcpy (miscloop,Path);
strcat (miscloop,"miscloop.dat");
strcpy (danglef,Path);
strcat (danglef,"dangle.dat");
strcpy (int22,Path);
strcat (int22,"int22.dat");
strcpy (int21,Path);
strcat (int21,"int21.dat");
strcpy (triloop,Path);
strcat (triloop,"triloop.dat");
strcpy (coax,Path);
strcat (coax,"coaxial.dat");
strcpy (tstackcoax,Path);
strcat (tstackcoax,"tstackcoax.dat");
strcpy (coaxstack,Path);
strcat (coaxstack,"coaxstack.dat");
strcpy (tstack,Path);
strcat (tstack,"tstack.dat");
strcpy (tstackm,Path);
strcat (tstackm,"tstackm.dat");
strcpy (int11,Path);
strcat (int11,"int11.dat");

}




int main(int argc, char* argv[]) {
	structure ct;
   datatable data;
   char loop2[maxfil],stackf[maxfil],tstackh[maxfil],tstacki[maxfil],
   	tloop[maxfil],miscloop[maxfil],danglef[maxfil],
      triloop[maxfil],int11[maxfil],
      int22[maxfil],int21[maxfil],coax[maxfil],DataPath[maxfil],
      tstackcoax[maxfil],coaxstack[maxfil],tstack[maxfil],tstackm[maxfil];
   char infile[maxfil],outfile[maxfil];
   char temp1[ctheaderlength],temp2[ctheaderlength];
   int i,k;

	if (argc!=3) {
   	cout << "EFN2 recalculates the free energy of an RNA structure in a CT file.\n";
    	cout << "Command Line Usage: efn2 [in ctfile] [out ctfile]\n";
      cout << "Enter the name of an existing CT file:";
      cin >> infile;
      cout << "Enter the name of an output ct file:";
      cin >> outfile;
   }
   else {
   	strcpy(infile,argv[1]);
      strcpy(outfile,argv[2]);

   }

   strcpy(DataPath,"");//require the data files to be located with the program

   getdat (loop2,stackf,tstackh,tstacki,tloop,miscloop,danglef,
   		int22,int21,coax,tstackcoax,coaxstack,tstack,tstackm,triloop,
         int11,DataPath); //get the names of the files

   //open the data files and check for error:
   if (opendat (loop2,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
   	coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,&data)==0) {


      cout << "One of the data files is missing\n";
   }

   //open the ct file
   openct(&ct,infile);

   //calculate the free energies
   efn2(&data,&ct);

   //re-sort the structures by free energy
   sortstructures(&ct);


   //remove the energy from the ct label (ct.ctlabel[1]) and
   //add a carriage return if needed

   for (k=1;k<=ct.numofstructures;k++) {
   	strcpy(temp2,ct.ctlabel[k]);
   	strcpy(ct.ctlabel[k],temp2+19);
   }


   //add a /n if needed to correct for PC-UNIX differences:
   for (k=1;k<=ct.numofstructures;k++) {
   	i = strlen(ct.ctlabel[k]);
   	strcpy(temp1,ct.ctlabel[k]+(i-1));
   	if (strcmp(temp1,"\n")) {
    		//need to add a carriage return:
      	strcat(ct.ctlabel[k],"\n");

   	}
   }




   //output the structures, re-ordered, in the out file
   ctout(&ct,outfile);


   return 0;
}

void errmsg(int err,int erri) {

if (err==30) {
	cout << "End Reached at traceback #"<<erri<<"\n";
   
}
if (err==100) {
	cout << "error # "<<erri;
   
}
switch (err) {
	case '1':
   	cout << "Could not allocate enough memory";
      break;
   case '2':
   	cout << "Too many possible base pairs";
      break;
   case '3':
   	cout << "Too many helixes in multibranch loop";
   case '4':
   	cout << "Too many structures in CT file";
   default:
   	cout << "Unknown error";
}
cin >> err;

return;

}
