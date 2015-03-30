

#include "algorithm.cpp"





void main(int argc, char* argv[]) {
	structure ct;
   datatable data;
   char infile[maxfil],outfile[maxfil];
   char temp1[ctheaderlength],temp2[ctheaderlength];
   int i,window,k;

	if (argc!=4) {
   	cout << "Filter will sort through a CT file to apply a Zuker Mfold Window.\n";
    	cout << "Command Line Usage: filter [in ctfile] [out ctfile] [window] \n";
      cout << "Enter the name of an existing CT file:";
      cin >> infile;
      cout << "Enter the name of an output ct file:";
      cin >> outfile;
      cout  << "Enter the window size:";
      cin >> window;
   }
   else {
   	strcpy(infile,argv[1]);
      strcpy(outfile,argv[2]);
      window = atoi(argv[3]);

   }


   //open the ct file
   openct(&ct,infile);


   filter(&ct, 99, 100000, window);

   //add a /n if needed to correct for PC-UNIX differences:
   for (k=1;k<=ct.numofstructures;k++) {
   	i = strlen(ct.ctlabel[k]);
   	strcpy(temp1,ct.ctlabel[k]+(i-1));
   	if (strcmp(temp1,"\n")) {
    		//need to add a carriage return:
      	strcat(ct.ctlabel[k],"\n");

   	}
   }

   //output the structures, filtered, in the out file
   ctout(&ct,outfile);



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
