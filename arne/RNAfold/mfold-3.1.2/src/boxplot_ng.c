/* Darrin Stewart and Michael Zuker */
/* May 17 2000 */
/* boxplot_ng.c version 2.0 */
/* execute with boxplot_ng $1 */
/*      */

/* Colors are now supplied by boxplot_col */
/* Look in boxplot_setcolor.inc for path to boxplot.col */
/* boxplot.col supplies the colors */
/* $1 is a plot file that already exists */
/* $1 must be the last paramater!*/
/* -------------------------------------------------------*/
/* Valid lower case parameters  */
/* See boxplot_ng.doc  also   */
/*                              */
/* -o name.out                  */
/*     Creates a postscript of gif file named name.out. */
/*     Default is name.gif or name.ps depending on -g.  */
/*     The default name is produced from name.plot by removing the .plot */
  
 
/* -t "This is the Title"   */
/* The title is displayed at the top of the gif or plot file */
/* The quotation marks are required */
/*                                  */
/* -g   indicates creation of  GIF file, default is a postscript file */
/* tgd must exist for this to work. */
/*     */
/* -go Turn off grid lines */
/*     */
/* -z left right top bottom */
/* left right top and bottom are integers specifying the boundaries of */
/* a  zoom region.  Some care is made by the program to make the region */
/* valid.  Regions below the main diagonal may be adjusted.*/
/* default value is no zoom. Show region defined by helices within plot file */
/*                                                     */
/* -f filter */
/* filter is an integer */
/* helices with length greater than or equal to filter will be displayed */
/* default value is 1 */
/*    */
/* -r resolution */
/*    use values from 50 to 400 , 72 works well*/
/*                                          */
/* -c colors values from 4 to 7 are allowed */
/*       */
/* -b enables transparent background for GIF files*/
/*       */
/* -p   PLot file is read as a probability file */

/*     */
/* -i increment */
/* increment is a floating point number with one digit right of the decimal */
/* point.  Positive values will be added to the optimal energy and */
/* energies within the range will be displayed */
/* negative values will be used as an absolute cutoff */
/* valid values are 5.1 3.0   or -113.5  */
/* values refer to Kcal/mole */ 
/* default value is the worst energy in the file */

/*----------*/
/* -c colors */
/* colors is an integer from 4 to 7 */
/* default value is 4 */
/* This option is not valid for a probability plot file */
/*------------------------*/
/* -r resolution */
/* resolution is an integer value of 72,110,200, or 300 */
/* The default value is 72 */
/* Multiply resolution by 8.5 for width, by 11 for height of gif */
/*  72  for    612 x  792 GIF */
/* 110 or    935 x 1210 GIF */
/* 200 or   1700 x 2200 GIF */
/* 300 or   2550 x 3300 GIF */
/* Resolutions from 50 to 400 will work */
/*                                      */
/* Note:  The file name must be the last parameter */
/* Other parameters will be searched for and processed whenever found */
/* Invalid parameters will probably be ignored and default values used */
/* ----------------------------------------------------------- */

/* The plot file may contain  energies or probabilities */
/* try box2plot_ng example or example2.plot */
/* For GIF files, this program creates a temporary file name.tgd */
/* The program runs tgd name.tgd to convert name.tgd */
/* into a gif file */
/* If there is an error, name.tgd may need to be deleted. */

/* change Jan 5, fixed points plotted with .ps file  and added */
/* outline is drawn 2 pixels wide for .gif file when resolution > 70 */
/* changes May 5, 1999 make dot sizes and colors within postscript file */
/* easy to change */
/* June 24, 1999 change, added degree, copyright, and delta to ps */


/* Change Sep 3, 1999 
   Removed unistd.h so program will work on Windows NT 
   */

/* Jan 3, 2000, adjust title to avoid long titles causing errors */

/* May 9, 2000 Extra digit added to probability plot data */
/* It will now read .0001  Before it could only read data */
/* as small as .001 */

/* May 16, 2000 Added -mi for mutual information */
/* Reorganized source code */

/* #include <unistd.h>*/  /* for reading wc -l of plot file */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h> /* for sqrt */
#include <time.h> /* for time within postscript file */
#include <limits.h> /* for INT_MAX */


#define TRUE    1
#define FALSE   0
#define MAIN_POST    1 /* create postcript for main or zoom */
#define ZOOM_POST    2 /* used in general_post, create_post */

#include "boxplot_setcolor.inc"   

char copyright_message1[]="boxplot_ng by D. Stewart and M. Zuker";
char copyright_message2[]= " Washington University";

int   global_helix_array_size;

int global_length;  /* len of RNA sequence, guessed at based on input file */

/* Variables that could easily be changed for personal preference are */
/* usually accompanied by $$$.  */
/*__________________________________________________________________*/


/* ______________________________________________________________*/   
/* Graphic  variables  -------------------------------------------*/
struct helix {
              int row;
              int column;
              int diagonal;
              int color;
              int energy;
              int length;
             };

int global_dumb_switch;
float global_color_increment; /* range of energy for each color */
int global_energy_cutoff; /* only probabilities <= */
			  /* global_energy_cutoff will be shown. */
                  /* note that probabilities have a negative sign */
                 /* -30 indicates probabilities >=.3 will be shown */
float global_energy_cutoff_stored;
int global_energy_cutoff_set;
int global_chain_len=1;/* chains smaller will not be displayed */

int global_points_plotted,global_points_plotted_zoom; /* indicates number of */
                         /* points plotted in each window */
/* _________________________________________________________*/
/* Zoom Variables __________________________________________*/

int display_l,display_r,display_t,display_b,display_w,display_h;
/* left,right,top,bottom,width and height for zoom window display */
/* They are in necleotide coordinates */
int global_display_from_arg;

/* File Variables _______________________________________________*/
char rec[90];
FILE *fp;   /* for .plot file */
char global_post_filename[120];
char global_plot_filename[120]; /* with .plot */
char global_output_name[120];
char global_name_of_file_img[120];
char global_name_of_file_ps[120];

int global_contains_prob;
int global_png_mode;
int global_jpg_mode;
int global_opt_prob_flag;
int global_mi_flag;
int global_clear_flag;
int global_save_tgd_flag;
int global_make_gifdat_flag;
int global_make_grid_flag;
int global_make_label_flag;
int global_label_i;
int global_label_j;
int global_post_adjust;
int global_type_of_file;
int global_resolution;
int global_adjust_col_flag=FALSE;

char global_title[150];

int global_optimal_energy;
int global_worst_energy;
int global_post_file; /* 1 for main, 2 for zoom */
                      /* used from text display */
struct helix *global_diag;

int *global_diag_count;
int *global_diag_start;/* 10000 is the maximum size */
int global_helices_in_plot_file; /* number of entries in file */

/* _________________________________________________________________*/


#include "boxplot_input.inc"
#include "boxplot_ps.inc"
#include "boxplot_img.inc"


/*____________________________________________________________________*/
/* Initialization Routines for nongraphics */

void initialize_data2(void)
{ global_energy_cutoff=global_worst_energy;
  if(global_contains_prob!=TRUE)
     global_color_increment=((float)(global_energy_cutoff -
                    global_optimal_energy))/(float)(global_number_of_colors-1);
  initialize_colors();

}

 
void initialize_data1(int argc, char **argv)
{ 
  open_file(argc,argv);
  initialize_len(FALSE,0);

}


/*____________________________________________________________________*/
 
void try_exit(int exit_choice)
{ if(exit_choice==1)
   {printf("\n Normal exit from boxplot_ng \n\n");
    free(global_diag);
   }
}



void set_energy_cutoff1(float cutoff_level)  /* for probabilities */
{  
   global_energy_cutoff=(int)(cutoff_level*-10000);
}

/* This function is called from the menu and sets the filter */
void set_chain_len(int chain_len_choice)
{
       global_chain_len=chain_len_choice;

}     
/*________________________________________________________________*/



void check_parameters(int argc,char **argv)
{ int number_of_colors,filter;
  int i,test;
  int title_set;
  char title2[140];
  float cutoff;
  strcpy(global_plot_filename,argv[argc-1]); /* set output file */
  argc=argc-1; /* last 1 is input file and has been used */
  global_type_of_file=2; /* default type is postscript */
  global_make_gifdat_flag=FALSE;/* default is no data file for gif zooming */ 
  global_make_grid_flag=TRUE; /* Default is to make grid lines */
  global_make_label_flag=FALSE; /* Default is no labels */
  global_label_i=0; /* default label positions */
  global_label_j=0;
  global_save_tgd_flag=FALSE; /* default is not to save the tgd file */
  global_png_mode=FALSE;
  global_jpg_mode=FALSE;
  global_clear_flag=FALSE;/* default is not a clear background for GIF */
  global_contains_prob=FALSE;
  global_energy_cutoff_set=FALSE;
  global_display_from_arg=FALSE;
  global_mi_flag=FALSE;
  /* Set Default Title */
  strcpy(global_title,global_plot_filename);
  /* replace the period with end of the title */
  fix_title(global_title); /* set title to default as input file without .plot*/
  /* set default output file */
  strcpy(global_plot_filename,global_title);
  strcpy(global_output_name,global_title);
  strcat(global_plot_filename,".plot");
  title_set=FALSE;
  filter=1; /* default is to show all helices */
  number_of_colors=4; /* default is 4 colors for dotplot*/
  /* set zoom region */
  global_resolution=72;
  i=1; /* check each argument */
  while(i<=(argc-1))
    {/* printf("\n argument %d is %s ",i,argv[i]);*/
      /* look for -p */  
     {test=strcmp(argv[i],"-m");
       if(test==0)
       {
        global_post_adjust=TRUE;
        i++;
       }
     else
     {test=strcmp(argv[i],"-p");
      if(test==0)
       {global_contains_prob=TRUE;
         i++; 
       }
     else 
      {test=strcmp(argv[i],"-g");
       if((test==0)&&(!global_dumb_switch))
       {/* printf("\n recognized -g flag"); */
        if(global_type_of_file!=3)
          global_type_of_file=1; /* create a gif file */
        i++;
       }
      else
       {test=strcmp(argv[i],"-pg"); /* look for -pg */
        if(test==0)
         {global_type_of_file=3; /* create a gif and a postscript file */
          i++;
         }
     else
       {test=strcmp(argv[i],"-mi"); /* look for -mi */
        if(test==0)
         {global_mi_flag=TRUE; /* Treat input data as mutual information */
          global_contains_prob=TRUE; /* treat as probability too */
          i++;
         }
        else 
        {test=strcmp(argv[i],"-png"); /* look for -png */
        if(test==0)
         {if(global_type_of_file!=3)
            global_type_of_file=1;
          global_png_mode=TRUE; /* create png, not gif */
          i++;
         } 
        else 
        {test=strcmp(argv[i],"-jpg"); /* look for -jpg */
        if(test==0)
         {if(global_type_of_file!=3)
            global_type_of_file=1;
          global_jpg_mode=TRUE; /* create jpg, not gif */
          i++;
         }
        else
	  {test=strcmp(argv[i],"-b");
           if(test==0) /* determine whether to use clear background */
              {global_clear_flag=TRUE;
	       i++;
              }
           else
             { test=strcmp(argv[i],"-d");/* determine whether to make .gifdat file */
                if(test==0)
                  {global_make_gifdat_flag=TRUE; /* create .gifdat file */
                   i++;
                  }
                else
                  {test=strcmp(argv[i],"-go");
                    if(test==0)
                      {global_make_grid_flag=FALSE; /* No grid lines  */
                       i++;
                      }
                    else
                      {test=strcmp(argv[i],"-l");
                        if(test==0)
                          {global_make_label_flag=TRUE;/* Label extra point and make energy file */
                            global_label_i=atoi(argv[i+1]);
                            global_label_j=atoi(argv[i+2]);
                             i+=3;
                          }
                        else
                           {test=strcmp(argv[i],"-s");
                            if(test==0)
                             { global_save_tgd_flag=TRUE; /* Save the tgd file */
                               i++;
                             }
                            else
                              {test=strcmp(argv[i],"-o");
                               if(test==0)
                                 {strcpy(global_output_name,argv[i+1]);
                                  i+=2;
                                  }
                               else
                                 {test=strcmp(argv[i],"-t");
				   if(test==0)
				     {strcpy(global_title,argv[i+1]);
                                           /* Title must follow -t */
				      global_title[75]='\0';
				      title_set=TRUE;
                                      i+=2;
				     }
                                   else
                                     {test=strcmp(argv[i],"-f"); /* Set filter */
				      if(test==0)
					{filter=atoi(argv[i+1]); /* filter must follow -f */
					  if(filter<1)
					    filter=1;
                                          i+=2;
                                        }
                                      else
					{test=strcmp(argv[i],"-c"); /* Set colors */
					  if(test==0)
					    {number_of_colors=atoi(argv[i+1]);
					     if(number_of_colors>8)
					       number_of_colors=8;
					     else
					       if(number_of_colors<2)
						 number_of_colors=4;
                                             i+=2;
                                            }
                                          else
					    {test=strcmp(argv[i],"-z");
					     if(test==0)
					       {if((i+4)>(argc-1))
						 {printf("\n Error!  improper use of -z \n");
                                                  global_display_from_arg=FALSE;
                                                  i++;
                                                 }
                                                else
						 {
                                                 global_display_from_arg=TRUE;
                                                 display_l=atoi(argv[i+1]);
					         display_r=atoi(argv[i+2]);
						 display_t=atoi(argv[i+3]);
						 display_b=atoi(argv[i+4]);
                                                 i+=5;
						 }
                                               }
                                         else {test=strcmp(argv[i],"-i");
					        if(test==0) /* set increment */
						  {sscanf(argv[i+1],"%f",&cutoff);
                                                   /* copy float into cutoff */
                                                   global_energy_cutoff_set=TRUE;
                                                   global_energy_cutoff_stored=cutoff;
                                                  i+=2;
						  }
                                                else
                                                  {test=strcmp(argv[i],"-r");
						    if(test==0) /* set gif resolution */
						      {global_resolution=atoi(argv[i+1]);
                                                       /* Resolution must follow -r */
						       if(global_resolution<20)
							 global_resolution=20;
						       else if(global_resolution>400)
							 global_resolution=400;
                                                       i+=2;
						      }
                                         else
					   {printf("\n\n !!!!!!!!!!!!!!Warning!!!!! \n");
                                            printf("   flag %s was not used!!",argv[i]);
                                            printf("\n !!!!!!!!!!!!!!!!!!!!!!!!!\n");
                                            i++;
                                           }
						  }
					 }		    }}}}}}}}}}}}}}}}}
   }
   /* adjust title */
  if((global_contains_prob)&&(global_number_of_colors<4))
    global_number_of_colors=4;
  if(title_set!=TRUE)
    {
     if(global_contains_prob==TRUE)
        {if(global_mi_flag)
          strcpy(title2,"Mutual Information Dotplot for ");
         else
          strcpy(title2,"Probability Dotplot for ");
        }
     else
        strcpy(title2,"Energy Dotplot for ");
     strcat(title2,global_title);
     strcpy(global_title,title2);
    }
  if(global_contains_prob)
    {if(filter>1)
      {printf("\n Error, The filter option is not available with probability");
       printf("\n or Mutual Information.");
       printf("\n All helices are length 1, i.e., only individual base pairs");
       printf("\n exist.\n\n");
       filter=1;
      }
    }
  set_chain_len(filter);
  global_number_of_colors=number_of_colors;
  if(global_contains_prob)
    {if(global_mi_flag)
       printf("\n Treating data as Mutual Information");
     else
       {if(global_opt_prob_flag)
          printf("\n Optimal Structure detected in Probability Data");
        else
          printf("\n Treating data as Probability");
       }
    }
 

}



void finish_making_post(char *text_window_input,char *title,

              int make_grid_flag,int make_label_flag,int label_i,
              int label_j,int adjust_post,int mi_flag,int opt_prob_flag)
{  int error;
   printf("\n Postscript file is     %s",
              text_window_input);
   strcpy(global_post_filename,text_window_input);
   if(strlen(global_post_filename)==0)
     {printf ("\n !!!! ERROR !!! Postscript filename has zero length \n");
       return;
     }
   error=open_post(global_post_filename);
   if(error)
     {printf("\n error creating postscript file ");
      return;
     }
   general_post(title,make_grid_flag,make_label_flag, label_i,label_j,
    adjust_post,0,NULL,NULL);
}


void make_output(void)
{int temp;
  float cutoff;
  if(global_energy_cutoff_set)
    {cutoff=global_energy_cutoff_stored;
     if(global_contains_prob==TRUE)
       {if(cutoff>3.)
           cutoff=3.0; 
	else if(cutoff<0.0)
	       cutoff=0.0;  
        set_energy_cutoff1(cutoff);
       }
     else
       {finish_setting_energy_cutoff(cutoff);
       }
    }
 if(global_display_from_arg==TRUE)
   {if((display_l>global_length)|| (display_l<1))
      display_l=1;
    if((display_r>global_length)|| (display_r<1))
      display_r=global_length;
    if((display_t>global_length)|| (display_t<1))
      display_t=1;
    if((display_b>global_length)|| (display_b<1))
      display_b=global_length;
    /* keep left smaller than right */
    if(display_r<display_l)
       {temp=display_r;
	display_r=display_l;
	display_l=temp;
       }
    /* keep top smaller than bottom */
    if(display_b<display_t)
       {temp=display_b;
	display_b=display_t;
	display_t=temp;
       }
   }
 else
   {display_l=1;
    display_t=1;
    display_r=global_length;
    display_b=global_length;
   } 
  display_w=display_r-display_l+1;
  display_h=display_b-display_t+1;
  if(global_type_of_file!=1)
     {strcpy(global_name_of_file_ps,global_output_name);
      strcat(global_name_of_file_ps,".ps");
      if((display_l==1)&&(display_t==1)&&
            (display_r==global_length)&&(display_b==global_length)&&
            (global_contains_prob!=TRUE))
       global_post_file=MAIN_POST;
     else
       global_post_file=ZOOM_POST;
     finish_making_post(global_name_of_file_ps,global_title,
               global_make_grid_flag,
               global_make_label_flag,
           global_label_i,global_label_j,global_post_adjust,
           global_mi_flag,global_opt_prob_flag);
     }
  if(global_type_of_file!=2)
     {strcpy(global_name_of_file_img,global_output_name);
      if(global_png_mode)
        strcat(global_name_of_file_img,".png");
      else
        {if((!global_jpg_mode)&&(!global_dumb_switch))
          strcat(global_name_of_file_img,".gif");
         else
           {strcat(global_name_of_file_img,".jpg");
            global_jpg_mode=TRUE;
           }
        }

      finish_making_gif(global_save_tgd_flag,global_clear_flag,
        global_resolution,
         global_name_of_file_img,global_title,global_make_gifdat_flag,
             global_make_label_flag,global_label_i,
             global_label_j,global_make_grid_flag,global_png_mode,
             global_jpg_mode,global_mi_flag,global_opt_prob_flag,0,
             NULL,NULL);
     }

}


int main(int argc, char **argv) 
{global_dumb_switch=FALSE;
  printf("\n\n\n boxplot_ng Version 2.0 \n\n");
  if((argc<2)||(argv[argc-1][0]=='-'))
    {printf("\n Insufficient arguments:\n");
     printf("\n[   Use:       boxplot_ng name                              ]");
     printf("\n[        Where name.plot is a plot file                     ]");
     printf("\n[                                                           ]");
     printf("\n[                Valid arguments                            ]");
     printf("\n[              (Default is to create a postscript file)     ]");
     printf("\n[ -b           Enable clear background for jpg,png file     ]");
     printf("\n                output.                                     ]");
     printf("\n[ -c colors    Specifies 2 to 8 colors for dots in plot.    ]");
     printf("\n[ -d           Create .gifdat file for WEB zooming          ]");
     printf("\n[ -f filter    Display only helices of length >= filter.    ]");
     if(!global_dumb_switch)
       printf("\n[ -g           Create only a GIF file for output.           ]");
     printf("\n[ -go          Grid lines off, default is on                ]");
     printf("\n[ -i increment Specify energy increment.                    ]");
     printf("\n[ -l i j       label i,j and create a file .gifeng listing  ]");
     printf("\n[               the energy at i,j                           ]");
     printf("\n[ -jpg         Create jpg output                            ]");
     printf("\n[ -m           Magnify dots in postscript output            ]");
     printf("\n[ -mi          Treate input data as Mutual Information      ]");
     printf("\n[ -o name_out  Specify name of output file.                 ]");
     printf("\n[ -p           Plot file is read as probability file        ]");
     printf("\n[              (default is energy)                          ]");
     printf("\n[ -pg          Create postscript and a png jpg file         ]");
     printf("\n                      Use with -png or -jpg                 ]");
     printf("\n[ -png         Create png output                            ]");
     printf("\n[ -r res       Specify resolution of png/jpg file(50 to 300)]");
     printf("\n[ -s           Save the tgd file for editing                ]");
     printf("\n[ -t \"TITLE\" ");
     printf("  Specify TITLE as title                       ]");
     printf("\n[ -z l r t b   Specify a zoom region: left,right,top,bottom ]");
     printf("\n[                                                           ]");
     printf("\n[    This program reads boxplot-3.0.col to determine colors.]");
     printf("\n[See boxplot_ng.doc for more information                    ]");
     printf("\n");
     exit(1);
    }
  initialize_data1(argc,argv);
  make_output();
  try_exit(1);
  return 0; 
}
