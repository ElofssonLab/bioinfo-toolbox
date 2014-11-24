/* boxplot_setcolor.inc  Version 2.0 */
/* May 12, 2000 */
/* Darrin Stewart and Michael Zuker */
/* Read boxplot-3.0.col to set colors and probabilities for each color */
#define TOTAL_COLORS 13
#define LAST_PROB_COLOR 8
#define COLOR_PATH1 "MFOLDLIB" /* environmental variable containing path */
/* to location of COLOR_FILE */
/* the working directory is searched first */
#define COLOR_PATH2 "MFOLD" /* secondary option */
#define COLOR_FILE "boxplot-3.0.col"


#define COLOR_BACKGROUND 10
#define COLOR_LABEL 9
#define COLOR_TEXT  0
#define COLOR_GRID  11
#define COLOR_OPTIMAL 1
#define COLOR_COPY  12

float global_color[TOTAL_COLORS][3]; /* 0 for red, 1 for green, 2 for blue */
/* color 0 is for text and outlines */
/* colors 1 through 7 are the 7 available colors for dots */
/* 8 is for dot labels */
/* 9 is background */
float global_prob_color[9][3];/* for 8 color probability */

float global_prob[5][8];/* define color for each range of probabilities */
/* 4 color scheme for 0 */
/* 5 color scheme for 1 */
/* 6 color scheme for 2 */
/* 7 color scheme for 3 */
/* 8 color scheme for 4 */
/* values 0 through 5 are the range of prob. for each color 1 to 7 */
/* 5 is used only with 7 colors */
FILE *cfp;

char color_filename[80];


const char *gif_color[]={"b0","bk","re","gr","ye","pu","br","bl","fl",
               "la","bak","gl","co"};

float global_ps_dot_min[8][2]; /* 1-7 are colors,*/
/*for 2nd  0 is dot magnifier, 1 is miniminum size in points */

int global_number_of_colors; /* assumes values from 4 to 7 */
/* It starts at 4, but can be changed */

void read_line(char *rec)
{ if(fgets(rec,90,cfp)==NULL)
    {printf("\n boxplot-3.0.col file is too short.\n");
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
   open_error(); /* was not found any one of the places */
}

void initialize_colors(void)
{ char rec[90];
  int i;
  float prob;
  float size_mag,size_min;
  open_color();
  read_line(rec); /* read first line , ignore it*/
  for(i=0;i<TOTAL_COLORS;i++) /* read and set colors */
    {read_line(rec);
     sscanf(rec,"%f %f %f",
        &global_color[i][0],&global_color[i][1],&global_color[i][2]);
     /*  printf("\n colors %f %f %f",
            global_color[i][0],global_color[i][1],global_color[i][2]); */
    }
  read_line(rec); /* read and ignore line */
  /* set colors 4 color scheme probability */
  for(i=0;i<=2;i++)
    {read_line(rec);
     sscanf(rec,"%f",&prob);
     global_prob[0][i]=prob;
     /* printf("\n global_prob[0][%d] = %f",i,prob); */
    }
/* set colors 5 color scheme probability */
  read_line(rec);
  for(i=0;i<=3;i++)
    {read_line(rec);
     sscanf(rec,"%f",&prob);
     global_prob[1][i]=prob;
    } 
/* set colors 6 color scheme probability */
  read_line(rec);
  for(i=0;i<=4;i++)
    {read_line(rec);
     sscanf(rec,"%f",&prob);
     global_prob[2][i]=prob;
    }
 /* set colors 7 scheme probability */
  read_line(rec);
  for(i=0;i<=5;i++)
    {read_line(rec);
     sscanf(rec,"%f",&prob);
     global_prob[3][i]=prob;
    }
/* set colors 8 scheme probability */
  read_line(rec);
  for(i=0;i<=6;i++)
    {read_line(rec);
     sscanf(rec,"%f",&prob);
     global_prob[4][i]=prob;
    }
/* set dot magfier and size for postscript */
  read_line(rec);
  read_line(rec);
  for(i=1;i<=7;i++)
    {read_line(rec);
     sscanf(rec,"%f%f",&size_mag,&size_min);
     global_ps_dot_min[i][0]=size_mag;
     global_ps_dot_min[i][1]=size_min;
     }
/* read 8 probability colors */
 read_line(rec); /* read first line , ignore it*/
  for(i=1;i<=LAST_PROB_COLOR;i++) /* read and set colors */
    {read_line(rec);
     sscanf(rec,"%f %f %f",
        &global_prob_color[i][0],&global_prob_color[i][1],
        &global_prob_color[i][2]);
    }
 fclose(cfp);
}

/* set color for probability */
/* probabilities have been altered by: */
/* energy(int)(-1*(10000.*rec_energy+.5));*/
int find_color_prob(int energy)
{int scheme; /* 0 for 4 colors, 3 for 7 colors */
 float prob; /* convert back from energy */
 int current_range;
 scheme=global_number_of_colors-4;
 prob=((float)energy)/-10000.; /* convert back to float */
 /* The following test is for speed */
 /* It eliminates the need for the bottom line */
 if(prob<=global_prob[scheme][global_number_of_colors-2])
   return global_number_of_colors;
 for(current_range=0;current_range<=(global_number_of_colors-2);
          current_range++)
    {if(prob>global_prob[scheme][current_range])
        {return current_range+1;
        }
    }
    /* the bottom line should never be reached */
 return global_number_of_colors; /* within range of lowest probability */
}





int find_color(int energy,int contains_prob,float color_increment,
           int optimal_energy)
{  int result;
   if(contains_prob==TRUE)
     return (find_color_prob(energy));
   if(energy==optimal_energy)
     return (1);
   result=(int)
     (((float)(energy-optimal_energy))/color_increment+
                                                      1.9999);
    return (result);
  
}




