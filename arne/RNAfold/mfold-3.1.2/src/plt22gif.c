/* Darrin Stewart  and Michael Zuker */
/* Jun 16, 1998 */
/* plt22gif.c */
/* This version was copied from zoom_plt22gif.c */
/* It should add the ability to color letters describing the bases */
/* based on how many other bases they can pair with */
/* This adds the -d and -c options */
/*   */

/* This plt22gif with extra options for zooming */
/* see plt22gif.doc for more information */
/* This version attempts to zoom */
/*  */
/*  */
/*  */
/*  */
/* name is a plt2 file that already exists */
/* name must be the last paramater!*/
/* The suffix .plt2 will be added if necessary */
/* -------------------------------------------------------*/
/* Valid lower case parameters  */
/* See pltt2gif.doc  also   */
/*       */
/* -b enables transparent background for GIF files*/
/*       */
/* -c Use appropriate colors for bases based on .ann or .ss-count */
/*      */
/* -d Substitute appropriately colored dots for the bases */
/*    Position is the position of the previous base, or end of */
/*    connecting line segments */
/*                              */
/* -n Read ss-count file rather than ann file for -c or -d option above */
/*                              */
/* -o name_out                  */
/*     Creates a gif file named name_out.gif */
/*     Default is name.gif  */
/*     The default name is produced from name.plt2 by removing the .plt2 */
/*       */
/*------------------------*/
/* -p forces portrait mode */
/*    The default is to guess between portrait and landscape  by reading */
/*     the file first to decide which works better.*/
/*    */
/* -q with -c or -d use a probability annotation file */
/*    file is name_k.ann or name.ann successively     */
/*    */
/* -r resolution */
/* resolution is an integer value between 20 and 350 */
/* Small values do not look as good, but large ones can produce files */
/* that are too big. */
/* The default value is 72 */
/* Multiply resolution by 8.5 for width, by 11 for height of gif */
/* or vice versa if in landscape mode */
/*  72  for    612 x  792 GIF */
/* 110  for    935 x 1210 GIF */
/* 200  for   1700 x 2200 GIF */
/* 300  for   2550 x 3300 GIF */
/*            */
/* -s      preserves tgd file */
/*         Edit this file and run tgd on it to produce a gif */
/*         Then delete the tgd file */
/*
 -x create a file named name.out.gif2bp
    This lists the locations of the basepairs, (dots, midpoint of line),
    within the gif image.
    The origin is upper left, with x across and y down.
    record consists of:
    x_location y_location base_1 base_2
    where base_1 < base_2 and all are integers */
/*                                      
 -z   s  x  y
      Given a location (x,y) in the original image, create
      a new image with data from location (x,y) in the center
      scaled by s.  (s=1 simply shifts the image)
      ( s=2 will double the size of the basepair dots in the region)
      ( s=.5 will make the dots half as big).
      (using -z 1 hor_size/2 vert_size/2 will produce the original image.)
      The original image has 0,0 as the upper left
      with x increasing across and y increasing down.
      The image is then scaled by s.  
      It is intended that an image is produced without the -z
      The program is then run again based on a (x,y) pixel point specified
      The -r value should be the same on both runs.
      s should be a floating point value greater than zero and less
      than 100.  
      The image size will be the same in both instances.
      x and y are integers >=0 and < image size in the appropriate
      direction.
      */
/* Note:  The file name must be the last parameter */
/* Other parameters will be searched for and processed whenever found */
/* Invalid parameters will probably be ignored and default values used */
/* ----------------------------------------------------------- */

/* For GIF files, this program creates a temporary file name.tgd */
/* The program runs tgd name.temp.tgd to convert name.temp.tgd */
/* into a gif file */
/* If there is an error, name.tgd may need to be deleted. */

/* To compile cc -o plt22gif plt22gif.c -lm */

/*

For the -x option to work:

The plt2 file must contain:
  1.  a)    CM BASE PAIRING LINES WITH BASE PAIRS (occurs exactly once)
      c)     MOV  real_number real_number i j ( occurs multiple times)
      d)     CIA  0.0                     ( as a pair)
        or 
             LI real_number real_number i j (occurs multiple times )
      h)  Some comment CM indicating the end of the basepairs.


The CM comments are essential.  When 1.a) is encountered, each circle or
line will be plotted and the .gif2bp file will be appended with x y i j
where x y are the location in pixels of the center of the segment or circle
representing basepair i j.


The name.gif2bp file contains:
max_dis
x y base_i base_j (occurs number_basepairs times)

where max_dis is the maximum distance in pixels from the center x,y
of a basepair that a user can click and expect to be told that
base i pairs with base j when viewing the gif file name.gif
This value has been multiplied by -1;
The contents after the max_dis are sorted.

Also,  naview must place the basepairs on the appropriate line and move
statements.

 */
typedef int boolean;

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h> 
#include "color_ann.h"

  
/* $$$$ */
#define global_close_width 2.5
 /* multiplty by font width to indicate how close */
/* in pixels a user must click to center of basepair within gif */
 
#define global_line_length 1.2 /* defines the distance from label to base */
/* It is in centimeters */
/* $$$$ */
/* Centimeters are scaled by, SA, or view matrix */
#define TRUE    1
#define FALSE   0
#define global_init_view_x_p 0.83 /* p for portrait, l for landscape */
#define global_init_view_y_p 0.83 /* make them bigger to shrink margins */
#define global_init_view_x_l 1.03 /* make them smaller to enlarge margins */
#define global_init_view_y_l 1.03 /* x and y can be different, but it skews*/
                                  /* the image. */
#define global_character_size 1.2 /* increase to make characters larger */
/* or smaller 1.1 is good */

/* TRUE or FALSE values below */
int global_portrait_flag; /* TRUE or FALSE */
int global_x_flag; /* TRUE or FALSE indicates when to create extra file */
                   /* providing the location of basepairs */
int global_save_tgd; /* non true deletes the tgd file */

int global_base_in_effect; /* The CM BASES # has been read */

int global_png_mode;
int global_jpg_mode;


/* The CM BASE PAIRING LINES has been read */
/* and each dot from CIA 0.0 or line indicates a basepair */


int global_extra_colors; /* true when -c or -d option is used */
int global_char_colors;  /* true when -c is used */
int global_char_dot;  /* true when -d is used */
int global_ann;  /* true with -c or -d unless -n flag is used */
int global_prob_ann; /* true when -q flag is used with -c and or -d */                 

/* Messages for top left corner of output */
char global_title_message[120];

char copyright_message1[]="plt22gif by D. Stewart and M. Zuker";
char copyright_message1_png[]="plt22png by D. Stewart and M. Zuker";
char copyright_message1_jpg[]="plt22jpg by D. Stewart and M. Zuker";
char copyright_message2[]= " Washington University";
char copyright_error_message[120]=" ";

 
char global_current_color[3];
char global_current_font[18];
int  global_font_width;
int  global_font_height;

float global_origin_x,global_origin_y,global_user_tran_x,global_user_tran_y;
float global_view_x,global_view_y;

float global_current_user_x;
float global_current_user_y;

int global_gif_width,global_gif_height; /* width and height of gif in pixels */
int global_first_run;
int global_largest_x;
int global_largest_y;
int global_smallest_x;
int global_smallest_y;
int global_hor_cen; /* offsets for centering */
int global_ver_cen; /* offsets for centering */

int global_distance_to_pair; /* must be within this to find basepair */
int global_left_base,global_right_base; /* bases at end of move line */
/* are stored here for -x option and CIA following the move */

/* for x option below */

float global_ppcm;  /* pixels per centimeter */
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
/*            Ignored for now, all lines are 2 or 1 pixel wide.  */
/*            lines shorter than 11 pixels are 1 pixel wide */
/* CTA display text */
/* CIA circle */
/* LI line */
/* CM comment */
/* CTX display text */
/* MOV move to user cooridinates */
/* COLOR sets color to RED, MAGENTA, or WHITE, which is shown as black */
/*  */
/* Variables that could easily be changed for personal preference are */
/* usually accompanied by $$$. */
/*__________________________________________________________________*/
/* zoom variables */
int global_z_flag,global_z_left,global_z_top,global_z_in_effect;
float global_z_scale;
float global_zoom_ppcm; /* global_z_scale*global_ppcm */
/* global_z_scale is s from comments */
/* global_z_flag is true when -z is specified */
/* global_z_left and global_z_top are based on (x,y), res, and scale  */
/* ______________________________________________________________*/   
/* Graphic  variables  -------------------------------------------*/
                         /* points plotted in each window */
/* _________________________________________________________*/

/* File Variables _______________________________________________*/
/* char rec[180]; */
FILE *xfp; /* for the x file */
FILE *plt2fp;   /* for .plt2 file */
FILE *giffp;  /* for postscript file */
char gif_filename[80];
char plt2_filename[80];
char x_filename[80];
char ann_filename[80];
char global_temporary_name[80];

#include "plt22gif_or_ps.h"
/* _________________________________________________________________*/

void set_max_min(int x,int y)
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




void open_plt2()
{

 printf("\n Trying to open  %s",plt2_filename);   
 if((plt2fp=fopen(plt2_filename,"r"))==NULL)
   {printf("\n Could not open plt2 file: %s \n",plt2_filename);
    exit(1);
   }   
}

/*_______________________________________________________________________*/
void try_exit(int exit_choice)
{ if(exit_choice==1)
   { printf("\n Normal exit from plt22gif \n\n");
     exit(0);
   }
}


int open_temporary(void)
{    /* open specified output file */
   if ((giffp = fopen(global_temporary_name,"w")) == NULL)
                                                 /* open a file */
         {printf ("\n * Could not open file:  %s", global_temporary_name);
          exit(1);
          }
    return (0);
}
int open_x(void)
{    /* open specified output file */
   printf("\n x_filename is %s",x_filename);
   if ((xfp = fopen(x_filename,"w")) == NULL)
         {                                       /* open a file */
          printf ("\n * Could not open file:  %s", x_filename);
          exit(1);
          }
    return (0);
}

/* write basepair in gif2bp file */
void write_pair(int x,int y,int base1,int base2)
  { fprintf(xfp,"%d %d %d %d\n",x,y,base1,base2);
  }
/*___________________________________________________________________________*/
/* create a gif */
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
      return;
     }
 test=strcmp(string,"MAG");
 if(test==0)
     {strcpy(global_current_color,"ma");
      return;
     }
 test=strcmp(string,"GRE");
 if(test==0)
    {strcpy(global_current_color,"gr");
     return;
    }
 test=strcmp(string,"BLU");
 if(test==0)
     {strcpy(global_current_color,"bl");
      return;
     }
 test=strcmp(string,"YEL");
 if(test==0)
     {strcpy(global_current_color,"ye");
      return;
     }
 test=strcmp(string,"CYA");
 if(test==0)
    {strcpy(global_current_color,"cy");
     return;
    }
 strcpy(global_current_color,"bk");
}
void open_files(void)
{ int len;
  if(strlen(gif_filename)==0)
      {printf("\n \n improper filename ");
       exit(1);
      }
  strcpy(global_temporary_name,gif_filename);
  len=strlen(global_temporary_name);
  if(len>4)
    global_temporary_name[len-4]='\0';  
  strcat(global_temporary_name,".tgd");
  open_temporary();
}
void display_time(char *year)
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

void display_title_message(void)

{ char year[8];
  display_time(year);
  /* Try not to make text overwrite image below */
  if(global_png_mode)
   fprintf(giffp,"string im gdFont6x9 %d %d \"%s\" ti\n",
                           6,1,copyright_message1_png);
  else
     {if(global_jpg_mode)
       fprintf(giffp,"string im gdFont6x9 %d %d \"%s\" ti\n",
                           6,1,copyright_message1_jpg);
      else
       fprintf(giffp,"string im gdFont6x9 %d %d \"%s\" ti\n",
                           6,1,copyright_message1);
     }
  fprintf(giffp,"string im gdFont6x9 %d %d \"  %s%s\" ti\n",
                           6,18,year,copyright_message2);
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
  fprintf(giffp,"string im gdFont8x13 %d %d \"%s\" ti\n",
                           6,32,global_title_message);
  if(strlen(copyright_error_message)!=0)
          fprintf(giffp,"string im gdFont8x13 %d %d \"%s\" re\n",
                           6,43,copyright_error_message);
  /* put in copyright symbol */
  fprintf(giffp,"string im gdFont12x24 5 9 \"O\" ti\n");
  fprintf(giffp,"string im gdFont6x13 8 15 \"C\" ti\n"); 
}
void initialize_gif(int clear_flag,int resolution,int portrait_flag,
        int no_name_change,char *specific_ann_file,int table_flag,
        char color_table_type)
{int black[3]; /* adjust black non-letters to make less noticable */
               /* than the letters */
  set_main_gif_colors();
  black[0]=main_color_table_gif[-1*COLOR_TEXT].red;
  black[0]=black[0]+20;
  if(black[0]>255)
      black[0]=main_color_table_gif[-1*COLOR_TEXT].red-20;
  black[1]=main_color_table_gif[-1*COLOR_TEXT].green;
  black[1]=black[1]+20;
  if(black[1]>255)
      black[1]=main_color_table_gif[-1*COLOR_TEXT].green-20; 
  black[2]=main_color_table_gif[-1*COLOR_TEXT].blue;
  black[2]=black[2]+20;
  if(black[2]>255)
      black[2]=main_color_table_gif[-1*COLOR_TEXT].blue-20;
  /* adjust white background to make less white */
  if((main_color_table_gif[-1*COLOR_BACKGROUND].red>253)&&
     (main_color_table_gif[-1*COLOR_BACKGROUND].green>253)&&
     (main_color_table_gif[-1*COLOR_BACKGROUND].blue>253))
    {main_color_table_gif[-1*COLOR_BACKGROUND].red=250;
     main_color_table_gif[-1*COLOR_BACKGROUND].green=250;
     main_color_table_gif[-1*COLOR_BACKGROUND].blue=250;
    }
  open_files(); /* this opens .plt2 temporary*/
  if(portrait_flag!=TRUE)
    {global_gif_width=(int)(11.*resolution+.5)-1;
     global_gif_height=(int)(8.5*resolution+.5)-1;
     if(global_z_flag==TRUE)/* convert from center to top left */
       {global_z_left=global_z_left-5.5*resolution/global_z_scale; 
        global_z_top=global_z_top-4.25*resolution/global_z_scale;
       }
    }
  else
    {global_gif_width=(int)(8.5*resolution)-1;
     global_gif_height=(int)(11.*resolution)-1;
     if(global_z_flag==TRUE)/* convert from center to top left */
       {global_z_left=global_z_left-4.25*resolution/global_z_scale; 
        global_z_top=global_z_top-5.5*resolution/global_z_scale;
       } 
    } 
  /* The items below define brushes for lines. $$$$$$$$$$$$$$$$$$$$$$$$*/
  /* create a brush for black lines */
  fprintf(giffp,"create bkim 2 2\n");
  fprintf(giffp,"colorallocate bkim white %d %d %d\n",
       main_color_table_gif[-1*COLOR_BACKGROUND].red,
       main_color_table_gif[-1*COLOR_BACKGROUND].green,
       main_color_table_gif[-1*COLOR_BACKGROUND].blue);
  fprintf(giffp,"colorallocate bkim bk2 %d %d %d\n",
              black[0],black[1],black[2]);
  fprintf(giffp,"line bkim 0 0 1 0 bk2\n");
  fprintf(giffp,"line bkim 0 1 1 1 bk2\n");
/* create a brush for Red lines */
  fprintf(giffp,"create reim 2 2\n");
  fprintf(giffp,"colorallocate reim white %d %d %d\n",
       main_color_table_gif[-1*COLOR_BACKGROUND].red,
       main_color_table_gif[-1*COLOR_BACKGROUND].green,
       main_color_table_gif[-1*COLOR_BACKGROUND].blue);
  fprintf(giffp,"colorallocate reim re %d %d %d \n",
         main_color_table_gif[-1*COLOR_RED_DOT].red,
         main_color_table_gif[-1*COLOR_RED_DOT].green,
         main_color_table_gif[-1*COLOR_RED_DOT].blue);
  fprintf(giffp,"line reim 0 0 1 0 re\n");
  fprintf(giffp,"line reim 0 1 1 1 re\n");
/* create a brush for Green lines */
  fprintf(giffp,"create grim 2 2\n");
  fprintf(giffp,"colorallocate grim white %d %d %d\n",
       main_color_table_gif[-1*COLOR_BACKGROUND].red,
       main_color_table_gif[-1*COLOR_BACKGROUND].green,
       main_color_table_gif[-1*COLOR_BACKGROUND].blue);
  fprintf(giffp,"colorallocate grim gr 0 204 0\n");
  fprintf(giffp,"line grim 0 0 1 0 gr\n");
  fprintf(giffp,"line grim 0 1 1 1 gr\n");
/* create a brush for yellow lines */
  fprintf(giffp,"create yeim 2 2\n");
  fprintf(giffp,"colorallocate yeim white %d %d %d\n",
       main_color_table_gif[-1*COLOR_BACKGROUND].red,
       main_color_table_gif[-1*COLOR_BACKGROUND].green,
       main_color_table_gif[-1*COLOR_BACKGROUND].blue);
  fprintf(giffp,"colorallocate yeim ye 224 224 0\n");
  fprintf(giffp,"line yeim 0 0 1 0 ye\n");
  fprintf(giffp,"line yeim 0 1 1 1 ye\n");
/* create a brush for Cyn lines */
  fprintf(giffp,"create cyim 2 2\n");
  fprintf(giffp,"colorallocate cyim white %d %d %d\n",
       main_color_table_gif[-1*COLOR_BACKGROUND].red,
       main_color_table_gif[-1*COLOR_BACKGROUND].green,
       main_color_table_gif[-1*COLOR_BACKGROUND].blue);
  fprintf(giffp,"colorallocate cyim cy 0 200 225\n");
  fprintf(giffp,"line cyim 0 0 1 0 cy\n");
  fprintf(giffp,"line cyim 0 1 1 1 cy\n");
  fprintf(giffp,"line cyim 0 1 1 1 cy\n");
/* create a brush for blue lines */
  fprintf(giffp,"create blim 2 2\n");  
  fprintf(giffp,"colorallocate blim white %d %d %d\n",
       main_color_table_gif[-1*COLOR_BACKGROUND].red,
       main_color_table_gif[-1*COLOR_BACKGROUND].green,
       main_color_table_gif[-1*COLOR_BACKGROUND].blue);
  fprintf(giffp,"colorallocate blim bl 0 0 255\n");
  fprintf(giffp,"line blim 0 0 1 0 bl\n");
  fprintf(giffp,"line blim 0 1 1 1 bl\n");
/* create a brush for Magenta lines */
  fprintf(giffp,"create maim 2 2\n");
  fprintf(giffp,"colorallocate maim white %d %d %d\n",
       main_color_table_gif[-1*COLOR_BACKGROUND].red,
       main_color_table_gif[-1*COLOR_BACKGROUND].green,
       main_color_table_gif[-1*COLOR_BACKGROUND].blue);
  fprintf(giffp,"colorallocate maim ma %d %d %d \n",  
       main_color_table_gif[-1*COLOR_BLUE_DOT].red,
       main_color_table_gif[-1*COLOR_BLUE_DOT].green,
       main_color_table_gif[-1*COLOR_BLUE_DOT].blue);
  fprintf(giffp,"line maim 0 0 1 0 ma\n");
  fprintf(giffp,"line maim 0 1 1 1 ma\n");
  /* switch to main image */
  /* Adjust colors below $$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
  /* And above to match  $$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
  fprintf(giffp,"create im %d %d\n",global_gif_width+1,global_gif_height+1);
                     /* declare image size in pixels */ 
  if(global_png_mode)
    printf("\n Creating png ");
  else
    {if(global_jpg_mode)
       printf("\n Creating jpeg ");
     else
       printf("\n Creating gif ");
    }
  printf("  %d by %d \n",global_gif_width+1,global_gif_height+1);
  fprintf(giffp,"colorallocate im white %d %d %d\n",
       main_color_table_gif[-1*COLOR_BACKGROUND].red,
       main_color_table_gif[-1*COLOR_BACKGROUND].green,
       main_color_table_gif[-1*COLOR_BACKGROUND].blue);
  /* first was background color */
  fprintf(giffp,"colorallocate im bk %d %d %d \n",
       main_color_table_gif[-1*COLOR_TEXT].red,
       main_color_table_gif[-1*COLOR_TEXT].green,
       main_color_table_gif[-1*COLOR_TEXT].blue);
     /* declare colors */
  fprintf(giffp,"colorallocate im re %d %d %d\n",  
         main_color_table_gif[-1*COLOR_RED_DOT].red,
         main_color_table_gif[-1*COLOR_RED_DOT].green,
         main_color_table_gif[-1*COLOR_RED_DOT].blue);/*black,red,green */
  fprintf(giffp,"colorallocate im gr 0 204 0\n"); /*yellow,purple,*/
  fprintf(giffp,"colorallocate im ye 224 224 0\n"); /*cyn,blue*/
  fprintf(giffp,"colorallocate im cy 0 200 225\n");
  fprintf(giffp,"colorallocate im bl 0 0 255\n");/* blue */
  fprintf(giffp,"colorallocate im bk2 %d %d %d \n",
                   black[0],black[1],black[2]);
  /* a black for single lines */
  fprintf(giffp,"colorallocate im ti %d %d %d\n",
       main_color_table_gif[-1*COLOR_COPYRIGHT].red,
       main_color_table_gif[-1*COLOR_COPYRIGHT].green,
       main_color_table_gif[-1*COLOR_COPYRIGHT].blue); /* for title */
  fprintf(giffp,"colorallocate im ma %d %d %d \n",
       main_color_table_gif[-1*COLOR_BLUE_DOT].red,
       main_color_table_gif[-1*COLOR_BLUE_DOT].green,
       main_color_table_gif[-1*COLOR_BLUE_DOT].blue);
        /* magenta, perhaps used for blue dot */
  fprintf(giffp,"colorallocate im xx %d %d %d\n",
       main_color_table_gif[-1*COLOR_GRAY_LABEL].red,
       main_color_table_gif[-1*COLOR_GRAY_LABEL].green,
       main_color_table_gif[-1*COLOR_GRAY_LABEL].blue);/* a gray for labels*/
  if((global_char_colors)&&(global_char_dot))
    {fprintf(giffp,"colorallocate im bkl %d %d %d\n",
       main_color_table_gif[-1*COLOR_BLACK_LETTER].red,
       main_color_table_gif[-1*COLOR_BLACK_LETTER].green,
       main_color_table_gif[-1*COLOR_BLACK_LETTER].blue);
       fprintf(giffp,"colorallocate im whl %d %d %d\n",
       main_color_table_gif[-1*COLOR_WHITE_LETTER].red,
       main_color_table_gif[-1*COLOR_WHITE_LETTER].green,
       main_color_table_gif[-1*COLOR_WHITE_LETTER].blue);
       /* for white/black letters on annotated bases */
    }
  if(global_extra_colors)
    {set_extra_gif_colors(giffp,global_prob_ann);
     define_map_ann(ann_filename,global_ann,global_prob_ann,
          no_name_change,specific_ann_file,table_flag,color_table_type,
               global_save_tgd,global_png_mode);
    }
  /* set up title for top left corner of screen */
  /* labels and lines to labels */
  if(clear_flag==TRUE)                      
    fprintf(giffp,"colortransparent im white\n");
  if(!global_jpg_mode) /* interlaced caused problems for xv with jpg */
    fprintf(giffp,"interlace im 1\n"); /* turn on interlace */
  fprintf(giffp,"setbrush im bkim\n");
  /*   */
}

void finish_gif_file(void)
{int error; /* execute tgd to convert file and then remove it */
  char string[100],action[200];
  display_title_message();
  strcpy(string,"gif im ");
  strcat(string,gif_filename); 
  strcat(string,"\n");
  fprintf(giffp,string);
  fclose(giffp);
  if(global_png_mode)
    printf("\n PNG  file: %s\n",gif_filename);
  else
    {if(global_jpg_mode)
       printf("\n JPG  file: %s\n",gif_filename);
     else
       printf("\n GIF  file: %s\n",gif_filename);
    }
  /* convert tgd file to gif */
  /*  strcpy(action,"/usr/people/dstewart/public_html/bin/tgd "); */
  if(global_png_mode)
    strcpy(action,"tgd_png ");
  else
    {if(global_jpg_mode)
       strcpy(action,"tgd_jpg ");
     else
       strcpy(action,"tgd ");
    }
  strcat(action,global_temporary_name);
  /*  printf("\n action is %s\n",action);*/
  error=system(action);
  if(error!=0)
    {printf("\n Error converting the file to GIF format");
     printf("\n error is %d",error);
     printf("\n The program executes tgd to convert file to GIF format.");
     printf("\n Either it was not installed or it was not installed ");
     printf("\n correctly. ");
     printf("\n \n See plt22gif.doc for more information");
     printf("\n ");
     printf("\n %s should be deleted",global_temporary_name);
    }
  else
    {if(global_save_tgd==TRUE)
      {printf("\n The file %s was preserved.",global_temporary_name);
       printf("\n Please delete it when finished!!!!\n");
      }
     else
       {     
        strcpy(action,"rm ");
        strcat(action,global_temporary_name);
        error=system(action);
        if(error!=0)
          printf("\n error deleting %s",global_temporary_name);  
        }
    }
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
  return y_screen;
}
int screen_to_gif_x(float x)
{ int x_screen;
  x_screen=(int)(global_ppcm*x); /* convert centimeters from left edge */
  /* to pixels from left edge */
  x_screen=x_screen+global_hor_cen;
  /* global_hor_cen moves image to center of gif */
  return x_screen;
}
int screen_to_gif_y(float y)
{ int y_screen;
  y_screen=(int)(global_ppcm*y); 
  /* convert centimeters to pixels */
  y_screen=y_screen+global_ver_cen;
  if(global_first_run!=TRUE)
    y_screen=global_gif_height-y_screen; /* omit to get upside down image */ 
 /* Center image vertically with global_ver_cen */
  /* Origin of gif is at top, is backwards from ps */
  return y_screen;
}

int user_to_gif_x(float x)
{/* convert to screen first */
 float x_screen;
 x_screen=user_to_screen_x(x);
 /* convert from screen to gif */
 return screen_to_gif_x(x_screen);
}
float zoom_screen_to_gif_x(float x)
{ /* compute non-zoom cooridinate first */
    float x_f;
    x_f=global_ppcm*x;
    x_f=x_f+global_hor_cen;
    /* now convert to zoom */
    return(int)(global_z_scale*(x_f-(float)global_z_left)+.5);
}
     
int zoom_user_to_gif_x(float x)
{ /* convert to screen first, same as nonzoom */
     float x_screen;
     x_screen=user_to_screen_x(x);
     /* convert from screen to gif */
     return zoom_screen_to_gif_x(x_screen);
}
float zoom_screen_to_gif_y(float y)
{ /* compute non-zoom cooridinate first */
    float y_f;
    y_f=global_ppcm*y;
    y_f=y_f+global_ver_cen;
    if(global_first_run!=TRUE) /* should always be executed */
       y_f=global_gif_height-y_f;
    /* now convert to zoom */
    return (int)(global_z_scale*(y_f-(float)global_z_top)+.5);
}
     
int zoom_user_to_gif_y(float y)
{ /* convert to screen first, same as nonzoom */
     float y_screen;
     y_screen=user_to_screen_y(y);
     /* convert from screen to gif */
     return zoom_screen_to_gif_y(y_screen);
}     
int user_to_gif_y(float y)
{/* convert to screen first */
 float y_screen;
 y_screen=user_to_screen_y(y);
 /* convert from screen to gif */
 return screen_to_gif_y(y_screen);

}
/* functions to create x file */
/*_____________________________________________________________*/

void open_gif2bp(void)
{/* This functin is called when CM BASE PAIRING LINES is reached */
  strcpy(x_filename,gif_filename);
  if(global_jpg_mode)
    {/* fix name for jpg to stay png */
      x_filename[strlen(x_filename)-2]='n';
      x_filename[strlen(x_filename)-3]='p';
    }
  strcat(x_filename,"2bp");
  open_x();
  global_base_in_effect=TRUE;
  printf("\n Base Pair file: %s",x_filename);
  fprintf(xfp,"%d \n",-1*global_distance_to_pair);  
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
            int base1,int base2,int *current_base,int outline_flag)
{/* draw a line in user coodinates from x1,y1 to x2,y2 */
 int gif_x1,gif_x2,gif_y1,gif_y2,gif_x,gif_y,test;
 float line_width;
 /* printf("\n start of LI base is  %d",*current_base); */ 
 if(global_z_in_effect==TRUE)
   {gif_x1=zoom_user_to_gif_x(x1);
    gif_x2=zoom_user_to_gif_x(x2);
    gif_y1=zoom_user_to_gif_y(y1);
    gif_y2=zoom_user_to_gif_y(y2);
    /* omit if not in region, This will fix most cases */
   if(((gif_x1<0)&&(gif_x2<0)) /* left of region */
      ||((gif_x1>global_gif_width)&&(gif_x2>global_gif_width))
        /* right of region,  above region */
      ||((gif_y1<0)&&(gif_y2<0))
      ||((gif_y1>global_gif_height)&&(gif_y2>global_gif_height)))
      {if(outline_flag==TRUE)
      if(*current_base<=global_total_bases_in_file)
         {
            {global_location_base_x[*current_base]=-5.0;
	    /*  printf("\n outline true, setting global_location_base %d ",*current_base); */
            *current_base=*current_base + 1;
            /* The -5 indicates to skip drawing this base later */ 
            }
         }
      return;
      }
   }			    
 else
   {
    gif_x1=user_to_gif_x(x1);
    gif_x2=user_to_gif_x(x2);
    gif_y1=user_to_gif_y(y1);
    gif_y2=user_to_gif_y(y2);  
    if(global_first_run==TRUE)    /* stop here when in */
      {set_max_min(gif_x1,gif_y1); /* check for center and portrait mode */
       set_max_min(gif_x2,gif_y2);
       return;
      }
   }
 if(global_z_in_effect==TRUE)
     line_width=.09*global_view_x*global_zoom_ppcm;
 else
     line_width=.09*global_view_x*global_ppcm; /* lines are .09 cm wide ? */
 /* any line less than 2.0 below will be one pixel wide. others are two */

 if(line_width<2.) /* draw single lines when need to be thin */
   {test=strcmp(global_current_color,"bk");
    if(test==0) /* draw single black lines darker */
       fprintf(giffp,"line im %d %d %d %d bk2 \n",gif_x1,gif_y1,
                          gif_x2,gif_y2);
    else 
       fprintf(giffp,"line im %d %d %d %d %s \n",gif_x1,gif_y1,
                          gif_x2,gif_y2,global_current_color);
   } 
 else /* draw a line 2 pixels wide */
    fprintf(giffp,"line im %d %d %d %d gdBrushed \n",gif_x1,gif_y1,
                          gif_x2,gif_y2);
 if(global_base_in_effect)
   write_pair((gif_x1+gif_x2)/2,(gif_y1+gif_y2)/2,base1,base2);
 if((outline_flag==TRUE)&&(*current_base<=global_total_bases_in_file))
   /* draw colored dots in outline mode for bases */
   {/* stop after global_total_bases_in file since the extra lines */
     /* are basepair lines */
    gif_x=gif_x2;
    gif_y=gif_y2;
    if(global_z_in_effect==TRUE)
        {if((gif_x<0)|| /* circle is left of region */
          (gif_x>=global_gif_width)|| /* right of region */
          (gif_y>=global_gif_height)|| /* below */
          (gif_y<0)) /* above */
           {global_location_base_x[*current_base]=-3.0;
	   /*  printf("\n LI_Extra setting global_location_base_x %d",*current_base); */
            *current_base=*current_base + 1;
            return;
           }
        }
     global_location_base_x[*current_base]=(float)gif_x;
     /*  printf("\n setting global_location_base_x %d",*current_base); */
     global_location_base_y[*current_base]=(float)gif_y;
     *current_base=*current_base + 1;
   }
}

void process_CSZ(float char_size)
{/* sets the character size in centimeters */ 
  /* probabably sets the width */
 float width;/* in pixels */
 float ppcm;
 if(global_z_in_effect)
   {char_size=char_size*global_z_scale;
    ppcm=global_zoom_ppcm;
   }
 else
    ppcm=global_ppcm;
 /* printf("\n char size is %f",char_size); */
 if(char_size<.15)
   width=2.3*char_size*ppcm*global_view_x*global_character_size+1.4;
 else
    if (char_size<.24)
      width=1.6*char_size*ppcm*global_view_x*global_character_size+1.4;
  else 
   width=.7*char_size*ppcm*global_view_x*global_character_size+1.4;
 /* $$$$$$$$$$$$$$$$$$$$$ */
 /* printf("\n width is %f",width); */
 global_distance_to_pair=(int)(global_close_width*width+.5)+5;
 /* set distance to find basepairs */
 /* printf("\n width in pixels: %f",width); */ 
 if(width<6.0)
   {strcpy(global_current_font,"gdFont5x8");
    global_font_width=5;
    global_font_height=8;
    return;
   }
 if(width<7.0)
   {strcpy(global_current_font,"gdFont6x9");
    global_font_width=6;
    global_font_height=9;
    return;
   }
 if(width<8.0)
   {strcpy(global_current_font,"gdFont7x13bold");
    global_font_width=7;
    global_font_height=13;
    return;
   }
 if(width>14)
   {strcpy(global_current_font,"gdFont12x24");
    global_font_width=12;
    global_font_height=24;
    return;
   }
 /* anything else is shown as 8x13bold */
 strcpy(global_current_font,"gdFont8x13bold");
 global_font_width=8;
 global_font_height=13;
}

void process_SA(float scale)
{/* scales the view matrix */
  if(scale!=0)
   {global_view_x=global_view_x/scale; 
    global_view_y=global_view_y/scale; 
   }
 else
   printf("\n attempt to divide by zero");
}
void process_ORI(float x,float y)
{ /* sets the origin of the screen cooridinates */
    global_origin_x=x;
    global_origin_y=y;
}
void process_CTA(float x,float y,char *string)
{ /* print the text centered at x y in screen coordinates */ 
  int len;
  int font_width;
  char current_font[20];
  int gif_x;
  int y_offset;
  if(global_z_in_effect==TRUE)
     return; /* no title for zoomed images */
  /* Make the font no smaller than 7 wide */
  font_width=global_font_width;
  if(font_width<7)
    {font_width=7;
     strcpy(current_font,"gdFont7x13bold");
    }
  else
    strcpy(current_font,global_current_font);
  remove_quotes(string);
  len=strlen(string);
  gif_x=font_width/2;
  gif_x=(global_gif_width-(len-1)*font_width)/2;
  /* The above seems to center the text over its position */
  if(global_portrait_flag==TRUE)
     y_offset=(int)(global_gif_height-1.8*global_ppcm);
  else
     y_offset=(int)(global_gif_height-1.*global_ppcm);
  /* Try not to make text overwrite image above */
  fprintf(giffp,"string im %s %d %d \"%s\" bk\n",current_font,
                           gif_x,y_offset,string);
}
void process_COL(char *color)
{
 fix_color(color);
 /* the following sets the brush to use for each color */
 fprintf(giffp,"setbrush im %sim\n",global_current_color);
}
void process_MOV(float x,float y)
{/* set the current user coordinates */
  /* seems to be used primarily before circles */
 global_current_user_x=x;
 global_current_user_y=y;
}
void process_CIA(float radius)
{ int gif_x,gif_y,radius_in_pix;
  float scale;
  /* draws a circle at current user cooridinates */
  /* radius of zero is converted to look nice */
  /* radius is not implemented */
  if(global_z_in_effect==TRUE)
    {scale=global_z_scale;/* control radius of circle below */
     radius_in_pix=(int)(global_ppcm*.4*global_view_x*scale);
     gif_x=zoom_user_to_gif_x(global_current_user_x);
     gif_y=zoom_user_to_gif_y(global_current_user_y);
     if((gif_x+radius_in_pix<0)|| /* circle is left of region */
        (gif_x-radius_in_pix>global_gif_width)|| /* right of region */
        (gif_y-radius_in_pix>global_gif_height)|| /* below */
        (gif_y+radius_in_pix<0)) /* above */
         return;
    }
  else
    {scale=1.; /* keep this line the same as the one of control above */
     radius_in_pix=(int)(global_ppcm*.4*global_view_x*scale);
     gif_x=user_to_gif_x(global_current_user_x);
     gif_y=user_to_gif_y(global_current_user_y);
    }
  /* radius of zero is drawn at .4 cm */
  /* basepair circles result from zero radius */ 
  if(radius!=0.0)
     {/* printf("\n weird radius is %f",radius); */
      radius_in_pix=(int)(global_ppcm*radius*global_view_x*scale);
     }
  else
    {/* $$$$$$ Line well above controls size of circles  .4 is good */
     if(global_base_in_effect==TRUE)
        write_pair(gif_x,gif_y,global_left_base,global_right_base);
    }
    fprintf(giffp,"arc im %d %d %d %d %d %d %s \n",gif_x,gif_y,
                          radius_in_pix,radius_in_pix,
                          0,360,global_current_color);
    if((gif_x>0)&&(gif_y>0)&&(gif_x<global_gif_width)
         &&(gif_y<global_gif_height))
    fprintf(giffp,"filltoborder im %d %d %s %s\n",gif_x,gif_y,
                            global_current_color,global_current_color); 
}

void process_CTX(float x,float y,char *string,int *current_base)
{/* draw text centered at user coordinates x,y */
  /* convert to screen coordinates*/
  /* call process_CTA */
 int len,gif_x,gif_y;
 float screen_x,screen_y;
 screen_x=user_to_screen_x(x);
 screen_y=user_to_screen_y(y);
 /* printf("\n screen_x is %f, screen_y is %f",screen_x,screen_y); */
 remove_quotes(string);
 len=strlen(string)/2;
 /* set location of each base when extra colors is set */
 /* on first run, location is not important, but a running base count */
 /* is used to later set the outline flag */

 /* convert to gif cooridinates and center vertically over position */
 if(global_z_in_effect==TRUE)
   {
    gif_x=zoom_screen_to_gif_x(screen_x);
    gif_y=zoom_screen_to_gif_y(screen_y);
    if((global_char_colors==TRUE)||(global_extra_colors!=TRUE))
     {gif_y=gif_y-7*global_font_height/16;
      if(((gif_x-10)<0)|| /* left of region */
       ((gif_x+10)>global_gif_width)|| /* right of region */
       ((gif_y-10)<0)|| /* above region */
       ((gif_y+10)>global_gif_height)) /* below region */
         {
          if(global_extra_colors==TRUE)
            {if(*current_base<= global_total_bases_in_file) 
               {
                global_location_base_x[*current_base]=(float)-4.0;
	        *current_base=*current_base+1; /* this base in not on screen */
		}
             /* The -4.0 indicates to skip drawing it later */
             /* Otherwise its true position is stored below */
            }
          return; /* Skip drawing text outside viewable region */
         }
     }
   }
 else
   {
    gif_x=screen_to_gif_x(screen_x);
    gif_y=screen_to_gif_y(screen_y);
    if((global_char_colors==TRUE)||(global_extra_colors!=TRUE))
        gif_y=gif_y-7*global_font_height/16; /* center vertically */
   }
 /* Center horizontally over the indicated position */
 if(((global_char_colors==TRUE)||(global_extra_colors!=TRUE))||
         (*current_base>global_total_bases_in_file) )
    gif_x=gif_x-len*global_font_width-global_font_width/2;
 /* on the first run update extreme values and keep a count of how */
 /* many bases are used */
 /* The count indicates whether it is an outline structure or not */
 /* If not, no base annotation will be done */
 if(global_first_run==TRUE)
   {set_max_min(gif_x,gif_y);
    if(global_extra_colors==TRUE)
      {*current_base=*current_base+1;
      /* printf("\n setting global_location_base_x %d",*current_base);*/ 
      }
    return;
   }
  /* The above seems to center the text over its position */
  /* tgd uses the bottom left of the character as a location */
  /* The coordinates seem to be different */

 /* following if line added for loop labels */
 /* It might be better without it for testing that annotation file matches */
 /* This way extra CTX lines are assumed to be loop labels and ignored */
 /* only annotation options make it to here */
 if(*current_base<=global_total_bases_in_file)
   {global_char_base[*current_base]=string[0];
     /* storing which base is at this position for later drawing */
     /* printf("\n setting base %d to %c",*current_base,string[0]); */ 
    if(global_char_dot==TRUE)
      {/* plan to draw a dot at computed position */ 
        if(global_z_in_effect==TRUE)
          { if((gif_x<0)|| /* circle is left of region */
               (gif_x>global_gif_width)|| /* right of region */
               (gif_y>global_gif_height)|| /* below */
               (gif_y<0)) /* above */
             {global_location_base_x[*current_base]=(float)-2.0;
	     *current_base=*current_base + 1;
             /* The -2 indicates not to draw this base */
              return;
	     }
	  }
      }   
   if(global_extra_colors==TRUE)
     {global_location_base_x[*current_base]=(float)gif_x;
     /*   printf("\n setting current _base_ x %d",*current_base); */ 
      global_location_base_y[*current_base]=(float)gif_y;
     }
     else
       fprintf(giffp,"string im %s %d %d \"%s\" bk\n",global_current_font,
                           gif_x,gif_y,string);
    *current_base=*current_base + 1;         
   }
   else
     {
       fprintf(giffp,"string im %s %d %d \"%s\" bk\n",global_current_font,
                           gif_x,gif_y,string);
     }
}

void process_TEX(float angle,float x,float y,char *string)
{/* plot string with left character at x,y in user coordinates*/
  /* given user coordinates */
 float x1,y1; /* location to place label */
 int gif_sx,gif_x1,gif_x2;
 int gif_sy,gif_y1,gif_y2;
 /* gif_s is x any y location for string */
 int len,i;
 char number[15];
 /* place label global_line_length away based on angle*/
 x1=x+global_line_length*(float)cos((double)(angle*-3.14159/180.));
            /* convert to radians */
 y1=y+global_line_length*(float)sin((double)(angle*-3.14159/180.));
 if(global_z_in_effect==TRUE)
   {gif_sx=zoom_user_to_gif_x(x1);
    gif_sy=zoom_user_to_gif_y(y1);
   }
 else
   {gif_sx=user_to_gif_x(x1);
    gif_sy=user_to_gif_y(y1);
   }
 remove_quotes(string);
 len=strlen(string);
 for(i=2;i<len;i++)
    {number[i-2]=string[i]; /* copy nonspace characters */
    }
 number[len-2]='\0';
 len=strlen(number);
 if(global_first_run==TRUE)
    {set_max_min(gif_sx-global_font_width*len/2,
                     gif_sy-global_font_height/2);
     return;
    }
 x1=x+.65*global_line_length*(float)cos((double)(angle*-3.14159/180.));
           /* convert to radians */
 y1=y+.65*global_line_length*(float)sin((double)(angle*-3.14159/180.));
 if(global_z_in_effect==TRUE)
    {gif_x1=zoom_user_to_gif_x(x1);
     gif_y1=zoom_user_to_gif_y(y1);
    }
 else
   {gif_x1=user_to_gif_x(x1);
     gif_y1=user_to_gif_y(y1);

   }
 x1=x+.15*global_line_length*(float)cos((double)(angle*-3.14159/180.));
          /* convert to radians */
 y1=y+.15*global_line_length*(float)sin((double)(angle*-3.14159/180.));
 if(global_z_in_effect==TRUE)
   {
    gif_x2=zoom_user_to_gif_x(x1);
    gif_y2=zoom_user_to_gif_y(y1);
   }
else
    {
    gif_x2=user_to_gif_x(x1);
    gif_y2=user_to_gif_y(y1);
   }
 fprintf(giffp,"string im %s %d %d \"%s\" xx\n",global_current_font,
         gif_sx-global_font_width*len/2,gif_sy-global_font_height/2,number);
 fprintf(giffp,"line im %d %d %d %d xx\n",gif_x1,gif_y1,gif_x2,gif_y2); 
} 

/* __________________________________________________________________*/
/*___________________________________________________________________*/
void uncode_line(char *rec,int *current_text_base,int outline_flag)
{ int len;
  int test;
  float float1,float2,float3,float4,float5,float6;
  char junk[101];
  char *string_pointer;
  char char_string[225];
  char char_string2[225];
  float angle; 
  int base1,base2; /*The number of bases */
  len=strlen(rec);

  /* printf("\n rec is %s",rec);*/  
  if(len<5)
    {/*printf("\n short line of length %d  was not processed. It was:  %s\n",
           len,rec);*/
     return;
    }
/* only TEX,CTX,ORI,SA,OD,LI are needed on the first run */
  test=strncmp("OD ",rec,3);
  if(test==0)
    {sscanf(rec,"%s%f%f",junk,&float1,&float2);
     process_OD(float1,float2);
     return;
    }
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
  test=strncmp("LI ",rec,3);
  if(test==0)
    {if(global_x_flag!=TRUE)
       {sscanf(rec,"%s%f%f%f%f%f%f",junk,&float1,&float2,&float3,
                                       &float4,&float5,&float6);
       process_LI(float1,float2,float4,float5,0,0,current_text_base,
            outline_flag);
       }
     else
       {sscanf(rec,"%s%f%f%f%f%f%f%d%d",junk,&float1,&float2,&float3,
                                 &float4,&float5,&float6,&base1,&base2);
        process_LI(float1,float2,float4,float5,base1,base2,
             current_text_base,outline_flag);
       }
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
  if(global_first_run==TRUE)
     return;
  /* commands below are not needed for first run */
  /* centering is based on lines and circles */
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
  test=strncmp("MOV",rec,3);
  if(test==0)
    {if(global_base_in_effect)
       sscanf(rec,"%s%f%f%d%d",
             junk,&float1,&float2,&global_left_base,&global_right_base);
     else
        sscanf(rec,"%s%f%f",junk,&float1,&float2);
     process_MOV(float1,float2);
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
      /* printf("\n current_record is %s",rec);  */
     if(global_x_flag==FALSE)
         return;
     global_base_in_effect=FALSE;
     test=strncmp(rec,"CM BASE PAIRING LINES WITH BASE PAIRS",37);
     if(test==0)
        {open_gif2bp();
        } 
     return;
    }
  printf("\n unknown command---->  %s\n",rec);
  return;
}
void process_plt2lines(int outline_flag,int *bases_used)
{ char rec[101];

  while(fgets(rec,99,plt2fp)!=NULL)
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
 char rec[120];
 global_smallest_x=10000; /* current_text_base is not needed here */
 global_smallest_y=10000;
 global_largest_x=-10000;
 global_largest_y=-10000; /* set these so that first read resets them; */
 open_plt2();
 *current_text_base=1;
 while(fgets(rec,99,plt2fp)!=NULL)
      {uncode_line(rec,current_text_base,FALSE);
      }
 rewind(plt2fp);
 /* printf("\n x dif was %f",largest_x-smallest_x); */
 /* printf("\n y dif was %f",largest_y-smallest_y); */
 dif=((float)(global_largest_x-global_smallest_x))/
              ((float)(global_largest_y-global_smallest_y));
 /* printf("\n Portrait test of width/height= %.3f   indicates ",dif);*/
 if(dif<1.2)/* 1.1 creates a tendency to be portrait $$$$$ */
    {/*printf(" Portrait mode\n");
     printf(" Values < 1.1 are  portrait\n",dif);*/ 
     return TRUE;
    }
 else
    {/*printf(" Landscape mode.\n");
     printf(" Values < 1.1 are  portrait\n",dif);*/ 
     return FALSE;
    }
}
/*_____________________________________________________________________________*/

void initialize_zoom(void)
{ global_zoom_ppcm=global_ppcm*global_z_scale;
  global_z_in_effect=TRUE;
}

void create_colored_dots_characters(int bases_used)
     /* This is done last in order to place them on top of the */
     /* previously drawn lines */
{int i;
 int gif_x;
 int gif_y;
 int radius_in_pix;
 char color[10];
 float scale;
 if(global_char_dot==TRUE)
   {    /* set radius */
	/* draw dots */
      if(global_z_in_effect==TRUE)
          scale=global_z_scale;/* control radius of circle below */
        else
          scale=1.; 
        radius_in_pix=(int)(global_ppcm*.75*global_view_x*scale);
      for(i=1;i<=bases_used;i++)
       { if(global_location_base_x[i]>=0.0)
	   {
            gif_x= (int)(global_location_base_x[i]+.5);
            gif_y= (int)(global_location_base_y[i]+.5);
            if(global_char_colors==TRUE)
	      {gif_x=gif_x+(9*global_font_width)/16;
               gif_y=gif_y+global_font_height/2;
              }

	    /*  printf("\n making dot for base %d of color %d",i,
            ann_to_color[i]); */
            
            fprintf(giffp,"arc im %d %d %d %d %d %d %s \n",gif_x,gif_y,
                   radius_in_pix,radius_in_pix,
                   0,360,color_table_char[ann_to_color[i]]);
	   fprintf(giffp,"filltoborder im %d %d %s %s\n",gif_x,gif_y,
            color_table_char[ann_to_color[i]],
            color_table_char[ann_to_color[i]]);
            if(global_char_colors==TRUE)
	      { if(ann_to_color[i]<WHITE_BLACK_SWITCH)
                  strcpy(color,"bkl");
                else
                  strcpy(color,"whl");
                fprintf(giffp,"string im %s %d %d \"%c\" %s\n",
                global_current_font,gif_x-(9*global_font_width)/16,
                   gif_y-(global_font_height)/2,
                 global_char_base[i],
                 color);
 /*   printf("\n placing base %d using character %c",i,global_char_base[i]);*/ 
              }
	   }
       } 
   }
 else if(global_char_colors==TRUE)
   {  for(i=1;i<=bases_used;i++)
	  {if((global_location_base_x[i])>=0.0)
           {gif_x=(int)(global_location_base_x[i]+.5);
            gif_y=(int)(global_location_base_y[i]+.5);
            fprintf(giffp,"string im %s %d %d \"%c\" %s\n",
            global_current_font,gif_x,gif_y,global_char_base[i],
            color_table_char[ann_to_color[i]]);
             }
          }
   }
}
void make_the_gif(int clear_flag,int resolution,int portrait_flag,
       int no_name_change,char *specific_ann_file,int table_flag,
       char color_table_type)
{ /* The initialization of global var. is only a precaution */
  /* against a poor ct file in most cases */
  int error;
  int bases_used; /* needed for making dots for bases in outline mode */
  char action[200];
  int image_width,image_height;
  int outline_flag; /* true indicates we are processing an outline only */
  global_portrait_flag=TRUE;
  global_view_x=global_init_view_x_p;
  global_view_y=global_init_view_y_p;
  global_origin_x=0.;       /* a page is about 21 cm wide */
  global_origin_y=0.;       /* or 28 cm wide */
  global_user_tran_x=0.;
  global_user_tran_y=0.;
  global_ppcm=resolution/2.54; /* pixels per real centimeter */ 
  global_hor_cen=0;/* for centering */
  global_ver_cen=0;
  global_first_run=TRUE;
  global_base_in_effect=FALSE;
  global_z_in_effect=FALSE;
  global_gif_width=(int)(8.5*resolution+.5)-1;
  global_gif_height=(int)(11.*resolution+.5)-1;
  bases_used=1;

  outline_flag=FALSE;
  global_portrait_flag=set_portrait_flag(&bases_used);/* also sets center */

  if((bases_used==1)&&(global_char_dot==TRUE))
    {/* printf("\n outline_flag set to true" ); */
        outline_flag=TRUE;
    }
  if((bases_used==1)&&(global_char_colors==TRUE))
    {printf("\n -----Error   No bases can be colored in -------");
     printf("\n ---- an outline structure.              ------- \n");
     global_char_colors=FALSE;
     global_extra_colors=FALSE;
     global_char_dot=FALSE;
     outline_flag=FALSE;
     strcat(copyright_error_message,"Cannot color bases in an outline structure, try the Dot option ");
    }
  if(portrait_flag==TRUE)
    global_portrait_flag=TRUE;
  if(global_portrait_flag!=TRUE)
    {global_gif_width=(int)(11*resolution+.5)-1;
     global_gif_height=(int)(8.5*resolution+.5)-1;
    }
  image_height=global_largest_y-global_smallest_y;
  image_width=global_largest_x-global_smallest_x;
  if(global_portrait_flag==TRUE)
     {global_hor_cen=(global_gif_width-image_width)/2;  /* 72*8.5=612 */
      global_ver_cen=(global_gif_height-image_height)/2; /* 72*11=792 */
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
    {global_hor_cen=(global_gif_width-
        image_width/global_init_view_x_p*global_init_view_x_l)/2;
     global_ver_cen=(global_gif_height-
        image_height/global_init_view_y_p*global_init_view_y_l)/2;
     global_hor_cen=global_hor_cen-global_smallest_x/
         global_init_view_x_p*global_init_view_x_l;
     global_ver_cen=global_ver_cen-global_smallest_y/
         global_init_view_y_p*global_init_view_y_l;
     global_view_x=global_init_view_x_l;
      /* The page is wider for landscape mode */
     global_view_y=global_init_view_y_l;/* use 1.05 later */
    }
  /* printf("hor cen is %d , ver cen is %d",global_hor_cen,global_ver_cen); */
  global_first_run=FALSE;
  global_base_in_effect=FALSE;
  initialize_gif(clear_flag,resolution,global_portrait_flag,no_name_change,specific_ann_file,table_flag,color_table_type);
  /*  temp file */
  /* Sets width in pixels, defines colors, brush and writes them */
  /* to temp file */
  /* width of page specified in naview affects this. */
  /* When naview is correct, this should be 1, but .9 looks better */
  strcpy(global_current_font,"gdFont10x20");
  global_font_width=10; /* in pixels */
  global_font_height=20;
  strcpy(global_current_color,"bk");
  if(global_z_flag==TRUE)
    initialize_zoom();
  bases_used=1;
  process_plt2lines(outline_flag,&bases_used);
  bases_used--;
  if(bases_used>1)
    create_colored_dots_characters(bases_used);
  fclose(plt2fp);
  finish_gif_file(); /* Write end of gif file, convert it to gif */
  /* and close gif file */
  if(global_x_flag==TRUE)
      {fclose(xfp); /* close file */
       strcpy(action,"sort -n -o ");
       strcat(action,x_filename); /* sort the file */
       strcat(action," ");
       strcat(action,x_filename); 
       strcat(action,"\n");
       error=system(action);
       if(error!=0)
           printf("\n Error sorting %s from plt22gif action was %s ",
            x_filename,action);
      }
}  

/*________________________________________________________________*/
 
/* ______________________________________________________________*/


void call_file_makers(int argc,char **argv)
{ int resolution;
  int i,test;
  char color_table_type; /* h or p or g for html, postscript,gif */
  char specific_ann_file[80];
  int no_name_change;
  int clear_flag,portrait_flag,table_flag;
  strcpy(gif_filename,argv[argc-1]); /* set output file */
  argc=argc-1; /* last 1 is input file and has been used */
  /* determine clear background option */
  clear_flag=FALSE;/* Default is not to use a clear background in the GIF */
  /* Adjust names */
  portrait_flag=FALSE; /* Default is not force portrait mode */
  global_save_tgd=FALSE; /* Default is not to save tgd file */
  fix_title(gif_filename); /* removes plt2 suffix */
  strcpy(plt2_filename,gif_filename);
  strcpy(ann_filename,gif_filename);
  strcat(plt2_filename,".plt2"); /* fix name of input file */
  global_char_colors=FALSE; /* True when -c is used */
  global_char_dot=FALSE;/* True when -d is used */
  global_png_mode=FALSE;
  global_jpg_mode=FALSE;
  global_extra_colors=FALSE; /* True when -c or -d has been used */
  global_ann=TRUE; /* Default is to use .ann not .ss-count when -c or -d is used */
                  /* -n will make this false */
  global_prob_ann=FALSE; /* default is not to use prob annotation with -c or -d switch */
                         /* -q makes this true */
  no_name_change=TRUE; /* Default is to use default names for -c or -d options */
  resolution=72; /* Default resolution for Gif file */
  global_x_flag=FALSE; /* Default is not make .gif2bp file for locating basepairs in .Gif */
  global_z_flag=FALSE; /* Default is not too zoom */
  table_flag=FALSE; /* Default is not to make table of colors */
  i=1;
  while(i<=(argc-1))
  {  /* look for -b */ 
  test=strcmp(argv[i],"-b");
   if(test==0)
         {clear_flag=TRUE; /* create a clear background gif */
          i++;
         }
   else   
     {test=strcmp(argv[i],"-o");
     if(test==0)
       {strcpy(gif_filename,argv[i+1]);
        i+=2;
       }
  else
    { /* set portrait mode for top of output */
     test=strcmp(argv[i],"-p");
     if(test==0)
      {portrait_flag=TRUE;
       i++;
      }
  else
     {/* set save tgd file mode  */
      test=strcmp(argv[i],"-s");
      if(test==0)
        {global_save_tgd=TRUE;
         i++;
        }
  else
	{ /* look for -c */ 
          test=strcmp(argv[i],"-c");
          if(test==0)
           {global_char_colors=TRUE;
            global_extra_colors=TRUE;
            i++;
	   } 
  else
	{ /* look for -png */ 
          test=strcmp(argv[i],"-png");
          if(test==0)
           {global_png_mode=TRUE;
            i++;
	   }
  else
	{ /* look for -jpg */ 
          test=strcmp(argv[i],"-jpg");
          if(test==0)
           {global_jpg_mode=TRUE;
            i++;
	   }
  else
      {/* Look for -d flag  */
       test=strcmp(argv[i],"-d");
       if(test==0)
         {global_char_dot=TRUE;
          global_extra_colors=TRUE;
          i++;
         }
  else
       {test=strcmp(argv[i],"-n");
        if(test==0)
         {global_ann=FALSE;
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
  else /* set resolution for gif file */
     {test=strcmp(argv[i],"-r");
      if(test==0)
	 {resolution=atoi(argv[i+1]); /* Resolution must follow -r */
	  if(resolution<20)
	    resolution=20;
	  else if(resolution>400)
	    resolution=400;
          i+=2;
	 }
  else
    {/* set xtra flag , makes the .gif2bp file for locating basepairs */
     test=strcmp(argv[i],"-x");
     if(test==0)
	{global_x_flag=TRUE;
         i++;
     
        }
  else /* Check for zoom */
    {test=strcmp(argv[i],"-z");
     if(test==0)
      {global_z_flag=TRUE;
       sscanf(argv[i+1],"%f",&global_z_scale);
       if(global_z_scale<=0.0)
	{printf("\n scale %f not allowed, must be >0.0",global_z_scale);
	 exit(1);
	}
       global_z_left=atoi(argv[i+2]);
       global_z_top=atoi(argv[i+3]);
       i+=4;
      }
  else
    {test=strcmp(argv[i],"-q");
     if(test==0)
      {global_prob_ann=TRUE;
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
  else
   {printf("\n Error !!!!");
    printf("\n Flag %s was not recognized \n",argv[i]);
    i++;
   }
    }}}}}}}}}}}}}}
  }
 if(global_png_mode)
   strcat(gif_filename,".png");
 else
   {if(global_jpg_mode)
     strcat(gif_filename,".jpg");
    else
      strcat(gif_filename,".gif");
   }
   /* Attatch proper suffix to default name or specified name */
   if((global_ann!=TRUE)&&(global_extra_colors!=TRUE))
       printf("\n -n flag cannot be used without -c or -d\n");
   if((global_prob_ann==TRUE)&&(global_extra_colors!=TRUE))
       printf("\n -q flag cannot be used without -c or -d\n");
   if((table_flag==TRUE)&&(global_extra_colors!=TRUE))
       printf("\n -t flag cannot be used without -c or -d\n");
   if((no_name_change==FALSE)&&(global_extra_colors!=TRUE))
       printf("\n -a flag cannot be used without -c or -d\n");
   if((global_prob_ann==TRUE)&&(global_ann==FALSE))
       {printf("\n -q and -n flag cannot be used together ");
        printf("\n -q flag was ignored ");
        global_prob_ann=FALSE;
       }
   make_the_gif(clear_flag,resolution,portrait_flag,no_name_change,specific_ann_file,table_flag,color_table_type);
}

int main(int argc, char **argv) 
{int i;
  printf("\n\n\n plt22gif \n\n");
  printf("\n Command line arguments were ...");
  for(i=0;i<argc;i++)
    {printf(" %s",argv[i]);
    }
  printf("\n");
  if((argc<2)||(argv[argc-1][0]=='-'))

    {printf("\n Insufficient arguments:\n");
     printf("\n[   Use:       plt22gif 'name'.plt2                           ]");
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
     printf("\n[ -b           Enable clear background for GIF file output.   ]");
     printf("\n[ -c           Use colored bases for annotation               ]");
     printf("\n[ -d           Use colored dots for annotation                ]");
     printf("\n[ -jpg         Create JPEG instead of GIF                      ]");
     printf("\n[ -n           read ss-count file, default is                 ]");
     printf("\n[                   .ann for -c or -d option above            ]");
     printf("\n[ -o name_out  Specify name of output file.                   ]");
     printf("\n[ -p           Force portrait mode                            ]");
     printf("\n[ -png         Create PNG instead of GIF                      ]");
     printf("\n[ -q           read probability annotation file,              ]");
     printf("\n[                  default is 'name'_k.ann, then name.ann     ]");
     printf("\n[                  for -c -d option above                     ]");
     printf("\n[ -r res       Specify resolution of gif file(20 to 400)      ]");
     printf("\n[ -s           Save the tgd file                              ]");
     printf("\n[ -t g         Make a gif table of colors,  use with -c,-d    ]");
     printf("\n[ -t h         Make a html table of colors,   use with -c,-d  ]");
     printf("\n[ -x           Create name.gif2bp of basepair locations       ]");
     printf("\n[ -z s x y     zoom with scale s > 0.0 and                    ]");
     printf("\n[              integers x,y at the center                     ]");
     printf("\n[                                                             ]");
     printf("\n[See plt22gif.doc for more information                        ]");
     printf("\n");
     exit(1);
    }
  call_file_makers(argc,argv);
  try_exit(1);
  return 0; /* this line is never hit */ 
}



