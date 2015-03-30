/* copyright_ps.inc */
/* Darrin Stewart and Michael Zuker */
/* Feb 5, 1999 */

#include <time.h>

#define TRUE    1
#define FALSE   0

char copyright_message1_ps[]="by D. Stewart and M. Zuker";
char copyright_message2_ps[]= " Washington University";

void display_time_year_ps(char *year)
{char time_data[28];
 int i,length;
 time_t now;
 now=time(NULL);
 strcpy(time_data,ctime(&now));
 length=strlen(time_data);
 for(i=0;i<=3;i++)
   {year[i]=time_data[length-(5-i)];
   }
 year[4]='\0';
 
}
void extra_ps_char(FILE* postfp)
{/* for copyright, delta, degree */
  fprintf(postfp,"\n/ISOLatin1Encoding where {pop save true}{false} ifelse");
  fprintf(postfp,"\n/reencodeISO {"); 
  fprintf(postfp,"\n   dup length dict begin");
  fprintf(postfp,"\n        {1 index /FID ne {def}{pop pop} ifelse} forall"); 
  fprintf(postfp,"\n         /Encoding ISOLatin1Encoding def");
  fprintf(postfp,"\n        currentdict"); 
  fprintf(postfp,"\n    end");
  fprintf(postfp,"\n} def"); 
  fprintf(postfp,"\n/findISO {");
  fprintf(postfp,"\n    dup /FontType known {"); 
  fprintf(postfp,"\n        dup /FontType get 3 ne {");
  fprintf(postfp,"\n              dup /CharStrings known {"); 
  fprintf(postfp,"\n                 dup /CharStrings get /Thorn known {");
  fprintf(postfp,"\n                     true"); 
  fprintf(postfp,"\n                }{ false } ifelse");
  fprintf(postfp,"\n            }{ false } ifelse");
  fprintf(postfp,"\n       }{ false } ifelse");
  fprintf(postfp,"\n     }{ false } ifelse"); 
  fprintf(postfp,"\n} def");
  fprintf(postfp,"\n");
}


void display_copyright_message_ps(FILE *postfp,int portrait_flag,
    char *program_name)
{ char year[8];
  float ver_offset;
  if(portrait_flag==TRUE)
    ver_offset=700.;
  else
    ver_offset=560;
  extra_ps_char(postfp);
  display_time_year_ps(year);
  /* Try not to make text overwrite image below */
  fprintf(postfp,"/sf 7.5 def\n");
  fprintf(postfp,"/Helvetica findfont findISO { reencodeISO /Symbol-ISO exch definefont \n");
  fprintf(postfp,"               sf scalefont setfont }\n");
  fprintf(postfp,"               { sf scalefont setfont } ifelse\n");
  fprintf(postfp,"15 %f moveto",ver_offset+66.);
  fprintf(postfp," (%s %s) show\n",program_name,copyright_message1_ps);
  fprintf(postfp,"26 %f moveto",ver_offset+55);
  fprintf(postfp," (\251 %s  %s) show\n",year,copyright_message2_ps);

}
