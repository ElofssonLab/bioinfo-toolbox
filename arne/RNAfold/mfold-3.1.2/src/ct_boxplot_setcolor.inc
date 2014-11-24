/* ct_boxplot_setcolor.h */
/* July 6, 1998 */
/* Darrin Stewart and Michael Zuker */
/* Read ct_boxplot.col to set colors and probabilities for each color */

#define COLOR_PATH1 "MFOLDLIB" /* environmental variable containing path */
/* to location of COLOR_FILE */
/* the working directory is searched first */
#define COLOR_PATH2 "MFOLD" /* secondary option */
#define COLOR_FILE "ct_boxplot.col"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>



 /* 0 for red, 1 for green, 2 for blue */

/* color 0 is overlap */
/* color 1 is partial overlap */
/* colors 2 through 16 are the 15 available colors for dots */
/* 202 is for dot labels */
/* 203 is background */


FILE *cfp;

char color_filename[80];

void read_line(char *rec)
{ if(fgets(rec,90,cfp)==NULL)
    {printf("\n ct_boxplot.col file is too short.\n");
     exit(1);
    }
}
void open_error(void)
{printf("\n !!!!!!!!!!!!!!!!Error !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
 printf("\n \n The file %s could not be found",COLOR_FILE);
 printf("\n Either ");
 printf("\n1. Place %s in the current working directory",COLOR_FILE);
 printf("\n\n or \n\n2. Set environmental variable %s to point to the",
        COLOR_PATH1);
 printf("\n location of %s",color_filename);
 printf("\n\n or \n\n3. Set environmental variable %s ",COLOR_PATH2);
 printf("\n and %s will be searched for within /dat of the path of %s\n",
     COLOR_FILE,COLOR_PATH2);
 exit(1);
}

void open_color(void)
{ char *color_path;
  strcpy(color_filename,COLOR_FILE);
  if((cfp=fopen(color_filename,"r"))!=NULL)
     return; /* COLOR_FILE was found in pwd */
  printf("\n color file was not in working directory");
  color_path=getenv(COLOR_PATH1);
  if(color_path!=NULL)
    { strcpy(color_filename,color_path);
      strcat(color_filename,"/");
      strcat(color_filename,COLOR_FILE);
      if((cfp=fopen(color_filename,"r"))!=NULL)
        return; /* COLOR_FILE was found in COLOR_PATH1 */
    }
 color_path=getenv(COLOR_PATH2);
 if(color_path!=NULL)
    { strcpy(color_filename,color_path);
      strcat(color_filename,"/dat/");
      strcat(color_filename,COLOR_FILE);
      if((cfp=fopen(color_filename,"r"))!=NULL)
        return; /* COLOR_FILE was found in COLOR_PATH2 */
    }
  printf("\n near end of open_color");
   open_error(); /* was not found any one of the places */
}

void initialize_colors(float color_table[][3],int global_ct_files)
{ char rec[90];
  int i;
  open_color();
  read_line(rec); /* read first line , ignore it*/
  for(i=0;i<=16;i++) /* read and set colors */
    {read_line(rec);
     sscanf(rec,"%f %f %f",
        &color_table[i][0],&color_table[i][1],&color_table[i][2]);
    }
  for(i=17;i<=201;i++)
    {color_table[i][0]=color_table[16][0];
     color_table[i][1]=color_table[16][1];
     color_table[i][2]=color_table[16][2];
    }
  for(i=202;i<=205;i++) /* read and set colors */
    {read_line(rec);
     sscanf(rec,"%f %f %f",
        &color_table[i][0],&color_table[i][1],&color_table[i][2]);
    }
  /* if there are only 2 ct files, copy 5 to 2 */
  if(global_ct_files<3)
   {
    color_table[2][0]=color_table[5][0];
    color_table[2][1]=color_table[5][1]; 
    color_table[2][2]=color_table[5][2];
   }
  fclose(cfp);
}










