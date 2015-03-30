/* Darrin Stewart */
/* July 20, 1998 */
/* plt22ps.c */
/* execute with plt22ps name */
/* This version does automatic centering. */
/*  */
/*  */
/*  */
/*  */
/* Change Jan 30 to add -s option */
/* Aug 18, added -z option */
/* name is a plt2 file that already exists */
/* name must be the last paramater!*/
/* The suffix .plt2 will be added if necessary */
/* -------------------------------------------------------*/
/* Valid lower case parameters  */
/* See pltt2ps.doc  also   */
/*                         */
/* -c Use appropriate colors for bases based on .ann or .ss-count */
/*      */
/* -d Substitute appropriately colored dots for the bases */
/*    Position is the position of the previous base, or end of */
/*    connecting line segments */
/*                              */
/* -n Read ss-count file rather than ann file for -c or -d option above */
/*                                         */
/* -o name_out                             */
/*     Creates a ps file named name_out.ps */
/*     Default is name.ps                  */
/*     The default name is produced from name.plt2 by removing the .plt2 */
/*       */
/*------------------------*/
/* -p forces portrait mode */
/*    The default is to guess between portrait and landscape  by reading */
/*     the plt2 file first to decide which works better.*/
/* -q  with -c or -d use a probability annoation file */
/*     file is name_k.ann or name.ann sucessively     */
/*    */
/* -s real_number */
/*    scales the postscript image in the x and y directions by          */
/*    real number. */  
/*     example:  plt22ps -s .25 example */
/*         will produce a file .25 * 8.5 inches wide and */
/*                             .25 * 11 inches tall */
/*                        */                                
/*                                      */
/* -z mag_fac ix iy, zoom with scale mag fac about point (ix,iy) */
/*                  0,0 is bottom left, 612,792 is top right */ 
/* Note:  The file name must be the last parameter */
/* Other parameters will be searched for and processed whenever found */
/* Invalid parameters will probably be ignored and default values used */
/* ----------------------------------------------------------- */


/* To compile cc -o plt22ps plt22ps.c -lm */

/* This program is similar to plt22gif.c */

/* all centimeters are scaled according to SA within the plt2 file */

typedef int boolean;

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>   
#include <time.h> /* for time within postscript file */
#include "color_ann.h"
 
#define global_line_length 1.1 /* defines the distance from label to base */
#define global_character_size 1.60 /* adjust character size, 1.3 is good */

/* It is in centimeters */
/* $$$$ */

#define TRUE    1
#define FALSE   0
#define global_init_view_x_p 0.83 /* p for portrait, l for landscape */
#define global_init_view_y_p 0.83 /* make them bigger to shrink margins */
#define global_init_view_x_l 1.03 /* make them smaller to enlarge margins */
#define global_init_view_y_l 1.03 /* x and y can be different, but it skews*/

/* Messages for top left corner of output */

char global_title_message[120];

char copyright_message1[]="plt22ps by D. Stewart and M. Zuker";
char copyright_message2[]= " Washington University";
char copyright_error_message[120]=" ";
                                 


int global_portrait_flag; /* TRUE or FALSE */

int global_extra_colors; /* true when -c or -d option is used */
int global_char_colors;  /* true when -c is used */
int global_char_dot;  /* true when -d is used */
int global_ann;  /* true with -c or -d unless -n flag is used */
int global_prob_ann; /* true when -q flag is used with -c or -d */     
int global_radius_has_been_set;
int global_ann_radius_has_been_set;
char   global_current_color[3];
float  global_font_width;
float  global_line_width;

float global_origin_x,global_origin_y,global_user_tran_x,global_user_tran_y;
float global_view_x,global_view_y;

float global_current_user_x; /* This set only by MOV  */
float global_current_user_y; /* Perhaps it should be set by LI too */

float global_ppcm; /* using points per cm for postscript */

/* The program reads the plt2 file twice */
/* when global_first_run is TRUE, largest x,y, smallest,x,y are found */
/* They are then used to compute offsets for horizontal and vertical */
/* centering */
int global_first_run;
float global_largest_x;
float global_largest_y;
float global_smallest_x;
float global_smallest_y;
float global_hor_cen; /* offsets for centering */
float global_ver_cen; /* offsets for centering */

/* screen cooridinates=view matrix*(user coordinates-user_tran_vector)+ */
/* screen_tran_vec */
/*   */
/* screen_tran_vector is set by ORI */
/* user_translation_vector is set by OD*/
/* View matrix is scaled by dividing by SA */

/* Only the first three letters matter */
/* DEV  specifies a device and is ignored */
/* SA divides the view matrix for scaling */
/* ORI Sets the origin of screen coordinates */
/* CSZ sets the character size in centimeters */
/* OD set user translation vector */
/* BRI set brightness of line from 1 to 5, controls width */
/*            Ignored for now, line width determined by the view matrix  */
/* CTA display text */
/* CIA circle */
/* LI line */
/* CM comment */
/* CTX display text */
/* MOV move to user cooridinates */
/* COLOR sets color to RED, MAGENTA, or WHITE, which is shown as black */
/*  */
/* Variables that could easily be changed for personal preference are */
/* usually accompanied by $$$.  Colors in particular                */

/* ______________________________________________________________*/   

/* File Variables _______________________________________________*/
char rec[90];
FILE *plt2fp;   /* for .plt2 file */
FILE *psfp;  /* for postscript file */
char ps_filename[80];
char ann_filename[80];
char plt2_filename[80];

#include "plt22gif_or_ps.h"

/* _________________________________________________________________*/
/* The program needs max and min of all drawing */
/* It runs twice, this checks max and min as each point is drawn */
/* the first time */
/* These are used to find an offset for centering */




void set_max_min(float x,float y)
     /* update largest and smallest x and y for centering */
{ if(x<global_smallest_x)
    global_smallest_x=x;
  if(x>global_largest_x)
    global_largest_x=x;
  if(y<global_smallest_y)
    global_smallest_y=y;
  if(y>global_largest_y)
    global_largest_y=y;
}

void set_color_ps(int col)
{if(col>=0)
   fprintf(psfp," c_%d\n",col);
 else
   fprintf(psfp," c_%c\n",main_color_table_let[-1*col]);
}



/* strip off .plt2 from end of title */
void fix_title(char *title)
{int i;
 i=strlen(title);
 if(i>5)
   {i=i-5;
    if((title[i]=='.')&&(title[i+1]=='p')&&(title[i+2]=='l')&&
          (title[i+3]=='t')&&(title[i+4]=='2'))
          {title[i]='\0';
           return;
          }
    } 
    return;
}

void clean_title_ps(char *title)
{/* place backslash in front of ( or ) or " or ' for postscript file */
 char clean_title[100];
 int i;
 int j;
 int length;
 length=strlen(title);
 j=0;
 for(i=0;i<length;i++)
   {if((title[i]==40)||(title[i]==41)||(title[i]==34)||(title[i]==39))
     {clean_title[j]=92; /* should be backslash */
      j++;
     }
      clean_title[j]=title[i];
      j++;
   }
 clean_title[j]='\0';
 strcpy(title,clean_title); 
}
void open_plt2()
{

 printf("\n Trying to open  %s\n",plt2_filename);   
 if((plt2fp=fopen(plt2_filename,"r"))==NULL)
   {printf("\n Could not open plt2 file: %s \n",plt2_filename);
    exit(1);
   }   
}

/*_______________________________________________________________________*/
 
void try_exit(int exit_choice)
{ if(exit_choice==1)
   { printf("\n Normal exit from plt22ps \n\n");
     exit(0);
   }
}


int open_ps(void)
{    /* open specified output file */
   if ((psfp = fopen(ps_filename,"w")) == NULL)
                                                 /* open a file */
         {printf ("\n * Could not open file:  %s", ps_filename);
          exit(1);
          }
    return (0);
}

/*___________________________________________________________________________*/
/* create a ps */
/*___________________________________________________________________________*/
void remove_quotes(char *string) /* also removes line feeds  */
{int i,len;
 len=strlen(string);
 for(i=0;i<len;i++)
   {if(string[i]=='"')
     {string[i]=' ';
      if(i==(len-1))
         string[i]='\0';
     }
    if(string[i]==10)
     string[i]='\0';
   }
}

void fix_color(char *color)
{int test;      /* WHITE is converted to black */
 char string[4];
 strncpy(string,color,3);
 string[3]='\0';
 test=strcmp(string,"RED");
 if(test==0)
     {strcpy(global_current_color,"re");
      set_color_ps(COLOR_RED_DOT);
      return;
     }
 test=strcmp(string,"MAG");
 if(test==0)
     {strcpy(global_current_color,"ma");
      set_color_ps(COLOR_BLUE_DOT);
      return;
     }
 test=strcmp(string,"GRE");
 if(test==0)
    {strcpy(global_current_color,"gr");
     fprintf(psfp,"0.0 0.8 0.0 setrgbcolor\n");
     return;
    }
 test=strcmp(string,"BLU");
 if(test==0)
     {strcpy(global_current_color,"bl");
      fprintf(psfp,"0.0 0.0 .99 setrgbcolor\n");
      return;
     }
 test=strcmp(string,"YEL");
 if(test==0)
     {strcpy(global_current_color,"ye");
      fprintf(psfp,"0.8 0.8 0.0 setrgbcolor\n");
      return;
     }
 test=strcmp(string,"CYA");
 if(test==0)
    {strcpy(global_current_color,"cy");
     fprintf(psfp,"0.0 .78 .88 setrgbcolor\n");
     return;
    }
 strcpy(global_current_color,"bk");
 set_color_ps(COLOR_TEXT);

}
void open_files(void)
{ 
  if(strlen(ps_filename)==0)
      {printf("\n \n zero length filename from plt22ps");
       exit(1);
      }
  open_ps();
}

void display_time(void)
{time_t now;
 now=time(NULL);
 fprintf(psfp,"%s %s\n","%%",ctime(&now));
}
void display_time_year(char *year)
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



void display_title_message(int portrait_flag)
{ char year[8];
  float ver_offset;
  if(portrait_flag==TRUE)
       ver_offset=700.;
  else
       ver_offset=560.;
  display_time_year(year);
  /* Try not to make text overwrite image below */
  fprintf(psfp,"%s Copyright message ","%%");
  extra_ps_char(psfp);
  fprintf(psfp,"/Helvetica findfont\n");
  fprintf(psfp,"7.5 scalefont\n");
  fprintf(psfp,"setfont\n");
  fprintf(psfp,"15 %f moveto",ver_offset+66.);
  fprintf(psfp," (%s) show\n",copyright_message1);
  fprintf(psfp,"26 %f moveto",ver_offset+55);
  fprintf(psfp,"/sf 7.5 def\n");
  fprintf(psfp,"/Helvetica findfont findISO { reencodeISO /Symbol-ISO exch definefont \n");
  fprintf(psfp,"               sf scalefont setfont }\n");
  fprintf(psfp,"               { sf scalefont setfont } ifelse\n");
  fprintf(psfp," (\251 %s  %s) show\n",year,copyright_message2);
  if(global_prob_ann==TRUE)
    strcpy(global_title_message,"probability annotation");
  else
   {if(global_extra_colors==TRUE)
     {if(global_ann==TRUE)
        strcpy(global_title_message,"p-num annotation");
       else
        strcpy(global_title_message,"ss-count annotation");
     }
   } 
  /* Try not to make text overwrite image below */


  fprintf(psfp,"/Helvetica findfont\n");
  fprintf(psfp,"14 scalefont\n");
  fprintf(psfp,"setfont\n");
  fprintf(psfp,"15 %f moveto",ver_offset+40.);
  fprintf(psfp," (%s) show\n",global_title_message);
  if(strlen(copyright_error_message)!=0)
          { fprintf(psfp,"15 %f moveto",ver_offset+22.);
            fprintf(psfp," (%s) show\n",copyright_error_message);
          }
 fprintf(psfp,"%s End of Copyright message \n","%%");
  
}

void define_colors_ps(void)
{int stop_point;
 int i;
 if(global_extra_colors)
   {
    if(global_prob_ann)
     stop_point=LOG_NUM_COLORS;
    else
      stop_point=NUM_COLORS;
   }
 fprintf(psfp,"%s Define colors for postscript \n","%%"); 
 if(global_extra_colors)
   fprintf(psfp,"%s Edit c_0 through c_%d for annotation \n","%%",stop_point);
 fprintf(psfp,"%s c_k is background, c_t is for general text\n","%%");
 fprintf(psfp,"%s c_r is red dots, c_b is for blue dots\n","%%");
 fprintf(psfp,"%s c_w , c_l annotation color for letters on bases \n","%%");
 fprintf(psfp,"%s c_c for copyright \n","%%");
 fprintf(psfp,"%s Colors originally set in color_ann.h\n","%%");
 fprintf(psfp,"%s For dot basepairs adjust radius with c_rad far below\n","%%");
 if(global_extra_colors)
   fprintf(psfp,"%s For annotated base on circle, adjust radius with c_drad far below\n","%%");
 if(global_extra_colors)
   {
    for(i=0;i<=stop_point;i++)
      {fprintf(psfp,"/c_%d { %.3f %.3f %.3f setrgbcolor} def\n",
             i,  color_table_ps[i].red,
                 color_table_ps[i].green,
                 color_table_ps[i].blue);
      }
   }
 for(i=1;i<MAIN_COLORS;i++)
    {fprintf(psfp,"/c_%c { %.3f %.3f %.3f setrgbcolor} def\n",
             main_color_table_let[i],  
                 main_color_table_ps[i].red,
                 main_color_table_ps[i].green,
                 main_color_table_ps[i].blue);
    }

}

void initialize_ps(int portrait_flag,float scale,int no_name_change,
                   char *specific_ann_file,int zoom_flag,int zoom_x,
                   int zoom_y, float zoom_scale,int web_zoom_flag,
                   int table_flag,char color_table_type)
{int bb_right,bb_left,bb_top,bb_bottom;
 int zoom_hor_cent;
 int zoom_ver_cent;
 int temp;
 float bb_scale;
 if(web_zoom_flag)
   {if(portrait_flag)
     {zoom_y=792-zoom_y;
     }
    else
      {temp=zoom_x;
       zoom_x=zoom_y;
       zoom_y=temp;
      }
   }
 if(global_extra_colors==TRUE)
     {define_map_ann(ann_filename,global_ann,global_prob_ann,
             no_name_change,specific_ann_file,table_flag,
             color_table_type,FALSE,FALSE);
      set_extra_ps_colors(global_prob_ann);
     }
  else
     set_main_ps_colors();
  open_files(); /* this opens ps file*/
  if(zoom_flag!=TRUE)
    {bb_left=0;
     bb_bottom=0;
    }
  else
    {zoom_hor_cent=(int)(zoom_x*zoom_scale+.5);
     zoom_ver_cent=(int)(zoom_y*zoom_scale+.5);
     bb_left=zoom_hor_cent-306;
     bb_bottom=zoom_ver_cent-396;
    }
  bb_right=bb_left+612;
  bb_top=bb_bottom+792;
  bb_right=bb_right*scale+.5;
  bb_top=bb_top*scale+.5;
  bb_bottom=bb_bottom*scale+.5;
  bb_left=bb_left*scale+.5;
  fprintf(psfp,"%c!PS-Adobe-3.0 EPSF-2.0 \n",'%');
  if(scale!=1)
     {fprintf(psfp,"%sBoundingBox: %d %d %d %d\n","%%",
      bb_left,bb_bottom,bb_right,bb_top);
     }
 else
   {fprintf(psfp,"%sBoundingBox: 0 0 612 792\n","%%");
    bb_scale=612./(bb_right-bb_left);
    fprintf(psfp,"%f %f scale\n",bb_scale,bb_scale);
    fprintf(psfp,"%d %d translate\n",-1*bb_left,-1*bb_bottom);
   }
 fprintf(psfp,"%s Created from plt2 file by plt22ps.c\n","%%");
  fprintf(psfp,"%s Darrin Stewart and Michael Zuker\n","%%");
  display_time();
  if(zoom_flag==TRUE)
    { fprintf(psfp,"%s                This scale is for zooming\n","%%");
      fprintf(psfp," %f %f scale\n",zoom_scale,zoom_scale);
      fprintf(psfp,"%s This file contains the complete image\n","%%");
      fprintf(psfp,"%s Omit obove scales and translate and change bound box to 0 0 612 792\n",
         "%%");
      fprintf(psfp,"%s to return to unzoomed image\n","%%");
    }
  fprintf(psfp,"%s\n%s\n","%%","%%");
  fprintf(psfp,"%s           This scale controls size of output\n","%%");
  fprintf(psfp," %f %f scale\n",scale,scale);
  fprintf(psfp,"%s\n","%%");
  fprintf(psfp,"%s To change image size ....\n","%%");
  fprintf(psfp,"%s For full size image Use Bound ngBox:0 0 612 792\n","%%");
  fprintf(psfp,"%s                     and 1.0 1.0 scale\n","%%");
  fprintf(psfp,"%s For image .25 of full size, replace 612 with .25*612\n","%%");
  fprintf(psfp,"%s                replace 792 with .25*792 and replace\n","%%");
  fprintf(psfp,"%s                1.0 with .25: Other scales are similar.\n","%%");
  fprintf(psfp,"%s                Valid scales are .05 to 1.00\n","%%");
  fprintf(psfp,"%s ONLY scale and Bound ngBox lines above should be changed.\n",
       "%%");
  fprintf(psfp,"%s\n","%%");
  define_colors_ps();
  set_color_ps(COLOR_BACKGROUND);
  fprintf(psfp,"%s Set background color\n","%%");
  fprintf(psfp,"0 0 moveto\n");
  fprintf(psfp,"611 0 rlineto\n");
  fprintf(psfp,"0 791 rlineto\n");
  fprintf(psfp,"-611 0 rlineto\n");
  fprintf(psfp,"0 -791 rlineto\n");
  fprintf(psfp,"closepath\n");
  fprintf(psfp,"fill\n");
  if(portrait_flag==FALSE)
    {fprintf(psfp,"645 0 translate\n");
     fprintf(psfp,"90 rotate\n");
    }
  set_color_ps(COLOR_COPYRIGHT);  
  display_title_message(portrait_flag);
  set_color_ps(COLOR_TEXT);
}

void finish_ps_file(void)
{/* execute tgd to convert file and then remove it */
  fprintf(psfp,"showpage\n");
  fprintf(psfp,"%sEOF\n","%%");
  fclose(psfp);
  printf("\n postscript  file: %s\n",ps_filename);
}
/* conversion functions -------------------------------------------*/
float user_to_screen_x(float x)
{ float x_screen;
  x_screen=global_view_x*(x-global_user_tran_x)+global_origin_x;
  return x_screen;
}

float user_to_screen_y(float y)
{ float y_screen;
  y_screen=global_view_y*(y-global_user_tran_y)+global_origin_y; 
/*  y_screen=global_view_y*(y)-global_user_tran_y+global_origin_y; */
  return y_screen;
}

float screen_to_ps_x(float x)
{ float x_screen;
  x_screen=global_ppcm*x; /* convert centimeters from left edge */
  /* to points from left edge */
  x_screen=x_screen+global_hor_cen; /* 72  moves right 1 inch */;
  return x_screen;
}
float screen_to_ps_y(float y)
{  float y_screen;
  y_screen=global_ppcm*y; /* convert centimeters to points from bottom*/
  y_screen=y_screen+global_ver_cen;/* 105 moves up 105 points */ 
  /* origin of plt2 is at bottom */
  /* move up 1 inch from bottom by adding 72 here*/
  return y_screen;
}

float user_to_ps_x(float x)
{/* convert to screen first */
 float x_screen;
 x_screen=user_to_screen_x(x);
 /* convert from screen to postscript */
 return screen_to_ps_x(x_screen);
}

float user_to_ps_y(float y)
{/* convert to screen first */
 float y_screen;
 y_screen=user_to_screen_y(y);
 /* convert from screen to postscript */
 return screen_to_ps_y(y_screen);

}

/* process each line ______________________________________________*/
/*_________________________________________________________________*/
void process_OD(float float1,float float2)
{/* this point in user coordinates will be mapped to the origin */
  /* of screen coordinates */
  /* rotations will be performed about this point */
 global_user_tran_x=float1;
 global_user_tran_y=float2; /* sets user translation vector */
}
void process_LI(float x1,float y1,float x2,float y2,
                int *current_base,int outline_flag)
{/* draw a line in user coodinates from x1,y1 to x2,y2 */
 float ps_x1,ps_x2,ps_y1,ps_y2;
 float line_width;
 int test;
 ps_x1=user_to_ps_x(x1);
 ps_x2=user_to_ps_x(x2);
 ps_y1=user_to_ps_y(y1);
 ps_y2=user_to_ps_y(y2);
 /* printf("\n line %f,%f to %f,%f",ps_x1,ps_y1,ps_x2,ps_y2); */ 
 if(global_first_run==TRUE)    /* stop here when in */
   {set_max_min(ps_x1,ps_y1); /* check for center and portrait mode */
    set_max_min(ps_x2,ps_y2);
    return;
   }
 test=strcmp(global_current_color,"bk");
 if(test==0)
   line_width=.69*global_view_x+.1; /* set width for black lines */
 else /* $$$$$$$ */
   line_width=1.2*global_view_x+.15; /* set width for non black lines */
 if(global_line_width!=line_width)
   {global_line_width=line_width;
    fprintf(psfp,"%.2f setlinewidth\n",line_width);
   }
 if((outline_flag==TRUE)&&(*current_base<=global_total_bases_in_file))
    {global_location_base_x[*current_base]=ps_x2;
     global_location_base_y[*current_base]=ps_y2;
     *current_base=*current_base+1;
    }
 fprintf(psfp,"%.2f %.2f moveto ",ps_x1,ps_y1);
 fprintf(psfp,"%.2f %.2f lineto stroke \n",ps_x2,ps_y2);
} 
void process_CSZ(float char_size)
{/* sets the character size in centimeters */ 
  /* probabably sets the width */
 float width;/* in points */
 /* printf("\n character size is %f",char_size); */
 /* $$$$ make middle sizes a little larger */
 if(char_size<.15)
   width=2.1*char_size*global_ppcm*global_view_x*global_character_size+1.4;
 else
    if (char_size<.24)
      width=1.5*char_size*global_ppcm*global_view_x*global_character_size+1.4;
  else 
   width=char_size*global_ppcm*global_view_x*global_character_size+1.4;
 global_font_width=width;
 if(global_font_width>29)
   global_font_width=29;
 /* Above puts an upper bound on the font size */
 fprintf(psfp,"/Helvetica findfont\n");
 fprintf(psfp,"%.2f scalefont\n",global_font_width);
 fprintf(psfp,"setfont\n");
}

void process_SA(float scale)
{/* scales the view matrix */
  if(scale!=0)
   {global_view_x=global_view_x/scale; 
    global_view_y=global_view_y/scale;
   }
 else
   printf("\n attempt to divide by zero when setting scale with SA");
}
void process_ORI(float x,float y)
{ /* sets the origin of the screen cooridinates */
    global_origin_x=x;
    global_origin_y=y;
}
void process_CTA(float x,float y,char *string)
{ /* print the text centered at x y in screen coordinates */ 
  float len;
  float ps_x;
  float y_offset;
  float size_of_letters;
  size_of_letters=15.;/* $$$$ size of Title characters at bottom */
  remove_quotes(string);
  len=((float)strlen(string));
  /* center vertically */
  if(global_portrait_flag==TRUE)
     {y_offset=30.; /* display text based on center */
      ps_x=612./2.; /* $$$ page size */
     }
  else
     {y_offset=60;
      ps_x=792/2.; /* $$$$ location of title */
     }
  ps_x=ps_x-(len-2)*size_of_letters*.5*.51;/* .51 for character width? */
  /* .5 for using half of the length of the string */
  fprintf(psfp,"/Helvetica findfont\n");
  fprintf(psfp,"%.2f scalefont\n",size_of_letters);
  fprintf(psfp,"setfont \n");
  fprintf(psfp,"%.2f %.2f moveto",ps_x,y_offset);
  clean_title_ps(string);
  fprintf(psfp," (%s) show\n",string);
}
void process_COL(char *color)
{
 fix_color(color);
}
void process_MOV(float x,float y)
{/* set the current user coordinates */
  /* seems to be used primarily before circles */
 global_current_user_x=x;
 global_current_user_y=y;
 if(global_first_run==TRUE)
   set_max_min(user_to_ps_x(x),user_to_ps_y(y));
}
/* define radius of circle */
void define_radius(void)
{ fprintf(psfp,"%s Adjust the radius of dot base pairs here\n","%%");
  fprintf(psfp,"/c_rad { %.4f } def \n",global_ppcm*.179*global_view_x);
  global_radius_has_been_set=TRUE;
}

void process_CIA(float radius)
{ float ps_x,ps_y,radius_in_p;
  /* draws a circle at current user cooridinates */
  /* radius of zero is converted to look nice */
  /* radius is not implemented */
  ps_x=user_to_ps_x(global_current_user_x);
  ps_y=user_to_ps_y(global_current_user_y);
  /* radius of zero is drawn at .5 cm */
  /* basepair circles result from zero radius */ 
  if(radius!=0.0)
     {/* printf("\n weird radius is %.2f",radius); */
      radius_in_p=global_ppcm*radius*global_view_x;
      fprintf(psfp,"%.2f %.2f moveto\n",ps_x,ps_y);
      fprintf(psfp,"%.2f %.2f %.3f 0 360 arc closepath fill \n",
         ps_x,ps_y,radius_in_p);
      return;
     }
  /* set radius of circle below .21 is good */
  if(!global_radius_has_been_set) 
    {define_radius();
    }
  fprintf(psfp,"%.2f %.2f c_rad 0 360 arc closepath fill \n",
         ps_x,ps_y);
}
void process_CTX(float x,float y,char *string,int *current_base)
{/* draw text centered at user coordinates x,y */
  /* This code is used for most text */
  /* convert to screen coordinates*/
  /* convert to to ps coordinates  */
 float screen_x,screen_y;
 float len,ps_x,ps_y;
 /*  printf("\n pair is %d",pair); */
 screen_x=user_to_screen_x(x);
 screen_y=user_to_screen_y(y);
 remove_quotes(string);
 len=((float)strlen(string));
 ps_x=screen_to_ps_x(screen_x);
 ps_y=screen_to_ps_y(screen_y);
 if(global_char_dot!=TRUE)
   {
   ps_x=ps_x-len*global_font_width*3/8.;
   ps_y=ps_y-3*global_font_width/8;
   }
 /* convert and center, multiply by 1/2 seems better*/
 /* but 3/8 looks better */
 /* center vertically */
 /* The ps_x line affects centering $$$$$$$$$$$$$$$$$*/
 /* The above seems to center the text over its position */

 if(global_first_run==TRUE)
    {set_max_min(ps_x,ps_y);
     if(global_char_colors==TRUE)
       *current_base=*current_base+1;
     return;
    }
 if(*current_base<=global_total_bases_in_file)
  {if(global_char_colors==TRUE)
      {global_char_base[*current_base]=string[0];
      }
    if(global_extra_colors==TRUE)
     {global_location_base_x[*current_base]=ps_x;
      global_location_base_y[*current_base]=ps_y;
      *current_base=*current_base+1;
     }
    else  
     {fprintf(psfp,"%.2f %.2f moveto",ps_x,ps_y);
      fprintf(psfp," (%s) show\n",string);
     }
  }
 else  
  {fprintf(psfp,"%.2f %.2f moveto",ps_x,ps_y);
   fprintf(psfp," (%s) show\n",string);
  }
}

void process_TEX(float angle,float x,float y,char *string)
{/* plot string with left character at x,y in user coordinates*/
  /* given user coordinates */
 float x1,y1; /* location to place label */
 float ps_sx,ps_x1,ps_x2;
 float ps_sy,ps_y1,ps_y2;
 float line_width;
 /* gif_s is x any y location for string */
 int i,len;
 char number[15];
 /* place label global_line_length away based on angle*/
 x1=x+global_line_length*(float)cos((double)(angle*-3.14159/180.));
            /* convert to radians */
 y1=y+global_line_length*(float)sin((double)(angle*-3.14159/180.));
 ps_sx=user_to_ps_x(x1);
 ps_sy=user_to_ps_y(y1);
 remove_quotes(string);
 len=strlen(string);
 if((string[0]==' ')&&(string[1]==' '))
    {for(i=2;i<len;i++)
      {number[i-2]=string[i]; /* copy nonspace characters */
      }
     number[len-2]='\0'; /* TEX lines have 2 spaces at start? */
    }
 else
   { strcpy(number,string);
   }
 len=strlen(number);
 if(global_first_run==TRUE)
   {set_max_min(ps_sx-global_font_width*len/5.,
                 ps_sy-3*global_font_width/8.);
    return;
   }
 x1=x+.65*global_line_length*(float)cos((double)(angle*-3.14159/180.));
           /* convert to radians */
 y1=y+.65*global_line_length*(float)sin((double)(angle*-3.14159/180.));
 ps_x1=user_to_ps_x(x1);
 ps_y1=user_to_ps_y(y1);
 x1=x+.25*global_line_length*(float)cos((double)(angle*-3.14159/180.));
          /* convert to radians */
 y1=y+.25*global_line_length*(float)sin((double)(angle*-3.14159/180.));
 ps_x2=user_to_ps_x(x1);
 ps_y2=user_to_ps_y(y1);
 line_width=.4*global_view_x*+.10;/* sets line  width for labels $$$$$*/
 if(global_line_width!=line_width)
   {global_line_width=line_width;
    fprintf(psfp,"%.2f setlinewidth \n",line_width);
   }
 fprintf(psfp,"%.2f %.2f moveto ",
                 ps_sx-global_font_width*len/5.,
                 ps_sy-3*global_font_width/8.);
 /* The above adjust the centering of the text at a position */
 /* len/4 would make more sense as len/2, i.e center the number */
 /* The above affects labels mostly $$$$$$$$$$$$$$$$$$$$ */
 fprintf(psfp,"(%s) show \n",number);
 fprintf(psfp,"%.2f %.2f moveto ",ps_x1,ps_y1);
 fprintf(psfp,"%.2f %.2f lineto stroke\n",ps_x2,ps_y2);
} 

/* __________________________________________________________________*/
/*___________________________________________________________________*/
void uncode_line(char *rec,int *current_text_base,int outline_flag)
{ int len;
  int test;
  float float1,float2,float3,float4,float5,float6;
  char junk[80];
  char *string_pointer;
  char char_string[225];
  char char_string2[225];
  float angle; 
  len=strlen(rec);
  if(len<5)
    {/*printf("\n short line of length %d  was not processed. It was:  %s\n",
           len,rec);*/
     return;
    }
  /* only TEX,CTX,ORI,SA,OD,LI are needed on the first run */
  test=strncmp("SA ",rec,3);
  if(test==0)
    {sscanf(rec,"%s%f",junk,&float1);
     process_SA(float1);
     return;
    }
  test=strncmp("ORI",rec,3);
  if(test==0)
    {sscanf(rec,"%s%f%f",junk,&float1,&float2);
     process_ORI(float1,float2);
     return;
    }
 test=strncmp("OD ",rec,3);
  if(test==0)
    {sscanf(rec,"%s%f%f",junk,&float1,&float2);
     process_OD(float1,float2);
     return;
    }
  test=strncmp("LI ",rec,3);
  if(test==0)
    {sscanf(rec,"%s%f%f%f%f%f%f",junk,&float1,&float2,&float3,
                                       &float4,&float5,&float6);
     process_LI(float1,float2,float4,float5,current_text_base,outline_flag);
     return;
    }
  test=strncmp("TEX ",rec,3);
  if(test==0)
    {sscanf(rec,"%s%f%f%f",junk,&angle,&float1,
                         &float2);
     string_pointer=strchr(rec,'"');
     process_TEX(angle,float1,float2,string_pointer);
     return;
    }
  test=strncmp("CTX",rec,3);
  if(test==0)
    {sscanf(rec,"%s%f%f",junk,&float1,&float2);
     string_pointer=strchr(rec,'"');
     strcpy(char_string,string_pointer);
     remove_quotes(char_string);

     sscanf(char_string,"%s",char_string2);
     process_CTX(float1,float2,char_string2,current_text_base);
     return;
    }
   if(global_first_run==TRUE)
     return;
  /* commands below are not needed for first run */
  /* centering is based on lines and circles */
  test=strncmp("MOV",rec,3);
  if(test==0)
    {sscanf(rec,"%s%f%f",junk,&float1,&float2);
     process_MOV(float1,float2);
     return;
    }
  test=strncmp("DEV",rec,3);
  if(test==0)
    {/* ignore device */
     return;
    }
  test=strncmp("CSZ",rec,3);
  if(test==0)
    {sscanf(rec,"%s%f",junk,&float1);
     process_CSZ(float1);
     return;
    }
  test=strncmp("CTA",rec,3);
  if(test==0)
    {sscanf(rec,"%s%f%f",junk,&float1,&float2);
     string_pointer=strchr(rec,'"');
     strcpy(char_string,string_pointer);
     process_CTA(float1,float2,char_string);
     return;
    }
  test=strncmp("COL",rec,3);
  if(test==0)
    {sscanf(rec,"%s%s",junk,char_string);
     process_COL(char_string);
     return;
    }
  test=strncmp("BRI",rec,3);
  if(test==0)
    {/* not used yet */
     return;
    }
  test=strncmp("CIA",rec,3);
  if(test==0)
    {sscanf(rec,"%s%f",junk,&float1);
     process_CIA(float1);
     return;
    }
  test=strncmp("CM ",rec,3);
  if(test==0)
    {/* comment */

     return;
    }
  printf("\n unknown command in plt2 file ---->  %s\n",rec);
  return;
}
void process_plt2lines(int outline_flag,int *bases_used)
{char rec[91];
  while(fgets(rec,90,plt2fp)!=NULL)
       {/* printf("%s\n",rec); */ 
        uncode_line(rec,bases_used,outline_flag);
       }
}

/* automatic check for portrait or landscape mode */
/*_____________________________________________________________________*/

int set_portrait_flag(int *current_text_base)
{/* based on the location of lines from LI */
  /* check smallest and largest for values of x and y */
  /* if the x difference is largest, use landscape mode */
  /* if the y difference is largest, use portrait mode */
 float dif;
 global_smallest_x=10000000;
 global_smallest_y=10000000;
 global_largest_x=-10000000;
 global_largest_y=-10000000; /* set these so that first read resets them; */
 open_plt2();
 while(fgets(rec,90,plt2fp)!=NULL)
      {uncode_line(rec,current_text_base,FALSE);
      }
 rewind(plt2fp);
 /* printf("\n x dif was %f",largest_x-smallest_x); */
 /* printf("\n y dif was %f",largest_y-smallest_y); */
 dif=(global_largest_x-global_smallest_x)/
              (global_largest_y-global_smallest_y);
 printf("\n Portrait test of width/height= %.3f   indicates ",dif);
 if(dif<1.2)/* 1.1 creates a tendency to be portrait, 1 is equal $$$$$ */
    {printf(" Portrait mode\n");
     printf(" Values < 1.2 are  portrait\n"); 
     return TRUE;
    }
 else
    {printf(" Landscape mode.\n");
     printf(" Values < 1.2 are  portrait\n"); 
     return FALSE;
    }
}
/*_____________________________________________________________________________*/


void define_ann_dot_radius(void)
{ fprintf(psfp,"%s Adjust the radius annotated base here\n","%%");
  fprintf(psfp,"/c_drad { %.4f } def \n",global_ppcm*.33*global_view_x);
  fprintf(psfp,"stroke\n");
  global_ann_radius_has_been_set=TRUE;
}

void create_colored_dots_characters(int bases_used,int outline_flag)
     /* This is done last in order to place them on top of the */
     /* previously drawn lines */
{int i;
 float ps_x,ps_y;
 float offset;
 offset=11./32.*global_font_width;
 if(global_char_dot==TRUE)
   { /* set radius */
	/* draw dots */
     if(!global_ann_radius_has_been_set)
         define_ann_dot_radius();
      for(i=1;i<=bases_used;i++)
       {    ps_x=global_location_base_x[i];
            ps_y=global_location_base_y[i];
            /* set the color */
            set_color_ps(ann_to_color[i]);
            fprintf(psfp,"%.2f %.2f c_drad 0 360 arc closepath fill \n",
                ps_x,ps_y);
            if(global_char_colors==TRUE)
               {if(ann_to_color[i]<WHITE_BLACK_SWITCH)
                  set_color_ps(COLOR_BLACK_LETTER);
                else
                  set_color_ps(COLOR_WHITE_LETTER);
                fprintf(psfp,"%.2f %.2f moveto",ps_x-offset,
                         ps_y-offset);
                fprintf(psfp," (%c) show\n",global_char_base[i]);
                if(i<bases_used)
                   fprintf(psfp," %.2f %.2f moveto\n",
                     global_location_base_x[i+1],global_location_base_y[i+1]);
               } 
       } 
   }
 else if(global_char_colors==TRUE)
   {  for(i=1;i<=bases_used;i++)
	  {ps_x=global_location_base_x[i];
            ps_y=global_location_base_y[i];
            set_color_ps(ann_to_color[i]);
            fprintf(psfp,"%.2f %.2f moveto",ps_x,ps_y);
            fprintf(psfp," (%c) show\n",global_char_base[i]);
          }
   }
}

void make_the_ps(int portrait_flag,float scale,
              int no_name_change,char *specific_ann_file,int  zoom_flag,
               int zoom_x,int zoom_y,float zoom_scale,int web_zoom_flag,
               int table_flag,char color_table_type)
{ /* The initialization of global var. is only a precaution */
  /* against a poor ct file in most cases */
  /* Perform initialization based on portrait flag */
  /* Next, check to determine whether portrait mode is needed */
  /* and check to find offsets for centering */
  /* the image is plotted twice. Once to find center, once for real */
  /* on first run, portrait mode is assumed */
  float image_width,image_height;/* in points */
  int bases_used; /* needed for making dots for bases in outline mode */
  int outline_flag; /* is the image only an outline */
  global_portrait_flag=TRUE;
  global_view_x=global_init_view_x_p; /* fix_view */
  /* These values are for portrait mode */
  global_view_y=global_init_view_y_p; /* depend on naview input,*/
  /* page size, and orientation */
  global_origin_x=0.; /* smaller views fit on page easier */    
  global_origin_y=0.;       
  global_user_tran_x=0.;
  global_user_tran_y=0.;
  global_ppcm=28.3;
  /* produces points per cm, for postscript 72 points per inch /2.54 */
  /* produced 28.346 or 28.3 points per cm */
  /* This is true cm not scaled cm as most are. */
  global_hor_cen=0;/* for centering */
  global_ver_cen=0;
  global_line_width=0.0;
  global_first_run=TRUE;
  bases_used=1;
  global_portrait_flag=set_portrait_flag(&bases_used);/* also sets center */
  if((bases_used==1)&&(global_char_dot==TRUE))
    {/* printf("\n outline_flag set to true" ); */
        outline_flag=TRUE;
    }
  if((bases_used==1)&&(global_char_colors==TRUE))
    {printf("\n -----Error   No bases can be colored in -----");
     printf("\n ---- an outline structure.              ----- \n");
     global_char_colors=FALSE;
     outline_flag=FALSE;
     global_char_colors=FALSE;
     strcat(copyright_error_message,"Cannot color bases in an outline structure, try the Dot option ");
    }
  global_line_width=0.0; /* This one is important */
  if(portrait_flag==TRUE)
     global_portrait_flag=TRUE;
  image_height=global_largest_y-global_smallest_y;
  image_width=global_largest_x-global_smallest_x;
  if(global_portrait_flag==TRUE)
     {global_hor_cen=(612-image_width)/2;  /* 72*8.5=612 */
      global_ver_cen=(792-image_height)/2; /* 72*11=792 */
      /* above is the location the text should start */
      /* smallest x smallest y is where it presently starts */
    /* if smallest_x is smaller than calculated cen position, */
    /* Use an offset for the difference */
      /* The calculation below is used as an offset to correct for the */
      /* difference */
      global_hor_cen=global_hor_cen-global_smallest_x;
      global_ver_cen=global_ver_cen-global_smallest_y;
      global_view_x=global_init_view_x_p;
      global_view_y=global_init_view_y_p;
     }
  else
    {global_hor_cen=(792-
        image_width/global_init_view_x_p*global_init_view_x_l)/2;
     global_ver_cen=(612-
        image_height/global_init_view_y_p*global_init_view_y_l)/2;
     global_hor_cen=global_hor_cen-global_smallest_x/
         global_init_view_x_p*global_init_view_x_l;
     global_ver_cen=global_ver_cen-global_smallest_y/
         global_init_view_y_p*global_init_view_y_l+30;
     /* The above 30 helps for some reason $$$$$$*/
     global_view_x=global_init_view_x_l;
      /* The page is wider for landscape mode */
     global_view_y=global_init_view_y_l;/* use 1.05 later */
     /* for horizontal, */
    }
  initialize_ps(global_portrait_flag,scale,no_name_change,specific_ann_file,
         zoom_flag,zoom_x,zoom_y,zoom_scale,web_zoom_flag,
         table_flag,color_table_type);
  /* Opens  ps file */
  /* Sets defines colors, and writes them */
  /* to ps file */
  /* The above are highly dependent on the values sent to Naview */
  /* Use .67 for older ct files*/
  /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
  /* Smaller values make the viewed image smaller and fit on the page easier */
  /* This seems to work, depends on values sent to naview */
  /* width of page specified in naview affects this. */
  /* When naview is correct, above should be 1, but .9 looks better */
  global_font_width=10.; /* in points */
  strcpy(global_current_color,"bk");
  global_first_run=FALSE;
  bases_used=1;
  process_plt2lines(outline_flag,&bases_used);
  bases_used--;
  if(bases_used>1)
    create_colored_dots_characters(bases_used,outline_flag);
  fclose(plt2fp);
  finish_ps_file(); /* Write end of ps file and close it */
}  

/*________________________________________________________________*/
 
/* ______________________________________________________________*/


void call_file_makers(int argc,char **argv)
{  int i,test;
  int portrait_flag;
  float scale;
  int zoom_flag;
  int table_flag;
  int web_zoom_flag;
  int no_name_change;
  char color_table_type; /* h or p or g for html, postscript,gif */
  char specific_ann_file[80];
  float zoom_scale;
  int zoom_x, zoom_y;
  web_zoom_flag=FALSE;
  table_flag=FALSE; /* default not to make color table */
  zoom_flag=FALSE; /* default to no zoom */
  strcpy(ps_filename,argv[argc-1]); /* set output file */
  argc=argc-1; /* last 1 is input file and has been used */
  /* replace the period with end of the name */
  fix_title(ps_filename); /* removes suffix */
  strcpy(plt2_filename,ps_filename);
  strcat(plt2_filename,".plt2"); /* fix name of input file */
  strcpy(ann_filename,ps_filename);
  scale=1.0;/* default scale for postscript output*/
  portrait_flag=FALSE; /* Default is not to force portrait mode */
  global_char_colors=FALSE; /* Default is no color annotat of bases*/
  global_extra_colors=FALSE; /* This is true only with -c or -d */
  global_char_dot=FALSE; /* Default is no color annot of dots */
  global_ann=TRUE;/* Default is to use .ann file with -c or -d */
                  /* rather than .ss-count */  
  global_prob_ann=FALSE; /* default is not to use prob annotation with -c or -d switch */
                         /* -q makes this true */
  global_radius_has_been_set=FALSE;
  global_ann_radius_has_been_set=FALSE;
  no_name_change=TRUE; /* Default is to automatically pick name */
  /* of annotation file , -a switches this to FALSE , must use -c or -d*/
  /* Interpret all command line switches */
  i=1;
  while(i<=(argc-1))
     {  
     /* specifiy output file */
     test=strcmp(argv[i],"-o");
     if(test==0)
       {strcpy(ps_filename,argv[i+1]);
        i+=2;
       }
     else /* Check to adjust scale */
        {test=strcmp(argv[i],"-s");
         if(test==0)
           {sscanf(argv[i+1],"%f",&scale);
            if(scale<.05)
              {scale=.05;
               printf("\n The scale was too small, was reset to .05\n");
              }
            if(scale>1.00)
              {scale=1.0;
               printf("\n The scale was too large, was reset to 1.0");
               }
            i+=2;
           } 
    else /* Check to adjust zoom */
        {test=strcmp(argv[i],"-z");
         if(test==0)
           {zoom_flag=TRUE;
            sscanf(argv[i+1],"%f",&zoom_scale);
            sscanf(argv[i+2],"%d",&zoom_x);
            sscanf(argv[i+3],"%d",&zoom_y);
            if(zoom_scale<.01)
              {scale=.01;
               printf("\n The magnification was too small, was reset to .01\n");
              }
            if(zoom_scale>200.00)
              {zoom_scale=200.0;
               printf("\n The magnification was too large, was reset to 200.0");
               }
            if(zoom_x<0)
               zoom_x=0;
            if(zoom_x<0)
               zoom_x=0;
            if(zoom_y<0)
               zoom_y=0;
            if(zoom_y>792)
               zoom_y=792;
            i+=4;
           } else /* Check to adjust zoom */
        {test=strcmp(argv[i],"-y");
         if(test==0)
           {web_zoom_flag=TRUE;
            sscanf(argv[i+1],"%f",&zoom_scale);
            sscanf(argv[i+2],"%d",&zoom_x);
            sscanf(argv[i+3],"%d",&zoom_y);
            if(zoom_scale<.01)
              {scale=.01;
               printf("\n The magnification was too small, was reset to .01\n");
              }
            if(zoom_scale>200.00)
              {zoom_scale=200.0;
               printf("\n The magnification was too large, was reset to 200.0");
               }
            if(zoom_x<0)
               zoom_x=0;
            if(zoom_x<0)
               zoom_x=0;
            if(zoom_y<0)
               zoom_y=0;
            if(zoom_y>792)
               zoom_y=792;
            i+=4;
           }
     else /* check to force portrait mode */
       {test=strcmp(argv[i],"-p");
         if(test==0)
          {portrait_flag=TRUE;
           i++;
           }
 else
    {test=strcmp(argv[i],"-t");
     if(test==0)
      {i++;
        if(i<=(argc-1))
	  {
           if((strlen(argv[i])==1)&&((!strcmp(argv[i],"g"))||
            (!strcmp(argv[i],"h"))))
            {color_table_type=argv[i][0];
             table_flag=TRUE;
            }
           i++;
          }

      }
    else /* check to use colored characters */
        {test=strcmp(argv[i],"-c");
          if(test==0)
           {global_char_colors=TRUE;
            global_extra_colors=TRUE;
            i++;
           }
    else /* check to use colored dots */
        {test=strcmp(argv[i],"-d");
         if(test==0)
          {global_char_dot=TRUE;
           global_extra_colors=TRUE;
           i++;
          }
   else /* check to use .ss-count file */
       {test=strcmp(argv[i],"-n");
        if(test==0)
          {global_ann=FALSE;
           i++;
          }
   else /* check to use probability annotation file */
       {test=strcmp(argv[i],"-q");
        if(test==0)
         {global_prob_ann=TRUE;
          i++;
         }
         
   else
     {/* look for -a flag to specify annotation file */
      test=strcmp(argv[i],"-a");
      if(test==0)
	{strcpy(specific_ann_file,argv[i+1]);
	 no_name_change=FALSE;
	 i+=2;
        }
      else 
       {printf("\n Error !!!!");
        printf("\n Flag %s was not recognized \n",argv[i]);
        i++;
       }
     }}}}}}}} }}
 }
  strcat(ps_filename,".ps");
   if((global_ann!=TRUE)&&(global_extra_colors!=TRUE))
       printf("\n -n flag cannot be used without -c or -d\n");
   if((global_prob_ann==TRUE)&&(global_extra_colors!=TRUE))
       printf("\n -q flag cannot be used without -c or -d\n");
   if((no_name_change==FALSE)&&(global_extra_colors!=TRUE))
       printf("\n -a flag cannot be used without -c or -d\n");
   if((global_prob_ann==TRUE)&&(global_ann==FALSE))
       {printf("\n -q and -n flag cannot be used together ");
        printf("\n -q flag was ignored ");
        global_prob_ann=FALSE;
       }
  if(web_zoom_flag)
    zoom_flag=TRUE;
  make_the_ps(portrait_flag,scale,no_name_change,specific_ann_file,zoom_flag,zoom_x,
                  zoom_y,zoom_scale,web_zoom_flag,table_flag,color_table_type);
}
int main(int argc, char **argv) 
{printf("\n\n\n plt22ps \n\n");
  if((argc<2)||(argv[argc-1][0]=='-'))
    {printf("\n Insufficient arguments:\n");
     printf("\n[   Use:       plt22ps 'name'.plt2                            ]");
     printf("\n[        Where 'name'.plt2 is a plt2 file                     ]");
     printf("\n[        and 'name'.ann exists for -c or -d and is a p-num    ]");
     printf("\n[        (default) or a probability annotation file (see -q   ]");
     printf("\n[        flag below).                                         ]");
     printf("\n[        If -n flag is used, the file 'name'.ss-count is used ]");
     printf("\n[        instead of 'name'.ann                                ]");
     printf("\n[                                                             ]");
   
     printf("\n[                                                             ]");
     printf("\n[                Valid arguments                              ]");
     printf("\n[ -a filename  Specify alternate name for annotation file     ]");
     printf("\n[              with -c or -d, .ann or .ss-count is assumed    ]");
     printf("\n[ -c           Use colored bases for annotation               ]");
     printf("\n[ -d           Use colored dots for annotation                ]");
     printf("\n[ -n           Read ss-count file, default is                 ]");
     printf("\n[                 ann for -c or -d option above.              ]");
     printf("\n[ -o name_out  Specify name of output file.                   ]");
     printf("\n[ -p           Force portrait output.                         ]");
     printf("\n[ -q           read probability annotation file,              ]");
     printf("\n[                  default is 'name'_k.ann, then name.ann     ]");
     printf("\n[                  for -c -d option above                     ]");
     printf("\n[ -s real_number Scales the image size by real_number         ]");
     printf("\n[                .05 to 1.00 is valid                         ]");

     printf("\n[ -t g         Make a gif table of colors,  use with -c,-d    ]");
     printf("\n[ -t h         Make a html table of colors,   use with -c,-d  ]");

     printf("\n[ -y scale x y   Zoom with magnification, scale about         ]");
     printf("\n[                pixel (x,y) from 72 resolution gif           ]");
     printf("\n[ -z scale x y  Zoom with magnification scale about the point ]");
     printf("\n[               (x,y) x,y are postive int, scale is >.01      ]");
     printf("\n[               (0,0) is bottom left, (612,792) is top right  ]");
     printf("\n[See plt22ps.doc for more information                         ]");
     printf("\n");
     exit(1);
    }
  call_file_makers(argc,argv);
  try_exit(1);
  return 0; /* this line is never hit */ 
}

