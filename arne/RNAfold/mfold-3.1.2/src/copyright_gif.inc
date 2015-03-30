/* copyright_gif.inc */
/* Darrin Stewart and Michael Zuker */
/* Feb 5, 1999 */

#include <time.h>

char copyright_message1_gif[]="by D. Stewart and M. Zuker";
char copyright_message2_gif[]= " Washington University";

void display_time_year_gif(char *year)
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
void display_copyright_message_gif(FILE *giffp,char *program_name)

{ char year[8];
  display_time_year_gif(year);
  /* Try not to make text overwrite image below */
  fprintf(giffp,"colorallocate im copyr  0 0 50\n"); /* copyright */
  fprintf(giffp,"string im gdFont6x9 %d %d \"%s %s\" copyr\n",
                           6,1,program_name,copyright_message1_gif);
  fprintf(giffp,"string im gdFont6x9 %d %d \"  %s%s\" copyr\n",
                           6,18,year,copyright_message2_gif);
  /* put in copyright symbol */
  fprintf(giffp,"string im gdFont12x24 5 9 \"O\" copyr\n");
  fprintf(giffp,"string im gdFont6x13 8 15 \"C\" copyr\n"); 
}

