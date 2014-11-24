

//SGI platform code for RNAstructure
# define bool int
# define true 1
# define false 0
# define TRUE 1
# define FALSE 0
//# define try /* not supported on SGI */
//# define catch /* not supported on SGI */
//# define (xalloc) /* not defined on SGI */
//#define floor ffloor
//#define floor /*floor*/
#define sgifix strcat(line," ");
#define binary in


#include <math.h>


int min(int i, int j) {

if (i>j) return j;

return i;

}

int max (int i,int j) {

if (i<j) return j;

return i;

}

int pow10(int i) {
int j,n;

if (i==0) return 1;
j = 1;

for (n=1;n<=i;n++) {
	j = 10*j;
}

return j;
}


/*void itoa (int column,char *number,int i )
{

int check,n,m,digit[10];

check = column;
for (n=1;check;n++) {
	digit[n] = column - (int (pow10(n)))*int (floor(column/int (pow10(n))));
	for (m=1;m<=n-1;m++) digit[n] = digit[n] - int (pow10(m))*digit[m];
	check = check - digit[n];
	digit[n] = digit[n]/int (pow10(n-1));
   number[n-1]=digit[n];
}

}*/



// removed 2002-11-03 NRM
/*void itoa (int x,char *ch,int i) 
{
float y;

y = float (x);
gcvt(y,6,ch);

}*/





/*double log(double i) {


return (log (double(i)));


}*/

/*int floor(float i) {

return (ffloor (i));

}*/


