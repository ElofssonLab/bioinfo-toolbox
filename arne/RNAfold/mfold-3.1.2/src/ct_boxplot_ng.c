/* Darrin Stewart  and Michael Zuker*/
/* Feb 4, 1999  */
/* ct_boxplot_ng.c  Version 2.0*/

/* execute with ct_boxplot_ng name.ct                   */
/* where name.ct is a ct file */
/* try ct_boxplot  test_1.ct */

/* This program creates a temporary file ct_boxplot_out.tgd */
/* as it creates the gif file                          */
/*  */
/* Arguments */

/*  -a Structure,   Display this structure in the lower triangle
    -b              Enable clear background for Gif output
    -d              Look for ct files in this directory, ends with /
    -g              Create gif output, default is postscript
    -l i j          Label the point i,j in the GIF or zoomed postscript image.
    -m  value       Magnify dot size by this floating point value.
    -o name_out 
                    Use name_out.gif or name_out.ps for output

   -pg              Create both gif and postscript output
   -r resolution    specify gif resolution from 20 to 300
   -s               Save the tgd file for editing.


    -t "TITLE"
    Specify a title for output
    Otherwise the first line of the first ct file is used. 
   -y s i j         Zoom with scale s about row i and column j
   -w create .gifdat for web zooming 
   -z left right top bottom        Specifies a zoom region */
/*  .gif or .ps is substituted for .ct on the input file */
/*  */ 

/* gif comments */
/* For GIF files, this program creates a temporary file name.tgd */
/* The program runs tgd name.tgd to convert name.tgd */
/* into a gif file */
/* If there is an error, temp.tgd may need to be deleted. */
/* resolution is an integer value between 20 and 300 */
/* Multiply resolution by 8.5 for width, by 11 for height of gif */
/*  72  for    612 x  792 GIF */
/* 110  for    935 x 1210 GIF */
/* 200  for   1700 x 2200 GIF */
/* 300  for   2550 x 3300 GIF */




/* Made changes October 18, 1999 to allow larger name in 
   first line of ct file and to center the title on the postscript file */




#include <stdio.h>
#include <string.h>
#include <stdlib.h>   
#include <math.h>
#include <time.h>
#include <limits.h>
#include "ct_boxplot_general.inc"
#include "ct_boxplot_gif.inc"
#include "ct_boxplot_ps.inc"
#include "ct_boxplot_read_ct.inc"
#include "ct_boxplot_setcolor.inc"

int global_png_mode;

int global_jpg_mode;



struct helix global_diag[MAXIMUM_HELICES]; 
char global_sequence_name[160];
int global_sequence_name_set;
float global_color[COLOR_TABLE_SIZE][3];

int global_len;  /* rightmost base that is part of a helix */

int *global_gray3_content;  
int global_lower_structure=-2; /* default to show partial and full overlap */
/* when right mouse button is clicked, transform coordinates into i,j */


/* ______________________________________________________________*/

/* Variables that could easily be changed for personal preference are */
/* usually accompanied by $$$.  Colors in particular                */

/* ______________________________________________________________ */

int global_last_structure; /* number of last structure */
/* It is the number of ct files specified +1 */
/* colors 2 to global_last structure */
/* are used for drawing dots */

/* Value to multiple by size of dots in all output */


					      /* for zoom */
/* ______________________________________________________________*/   
/* Graphic  variables  -------------------------------------------*/


/* _________________________________________________________*/



int global_zoom_labels_count ; /* the number of labels used */
int global_zoom_labels_row[50];
int global_zoom_labels_col[50];

char *global_zoom_labels_string[50]; /* value to place at x,y */

/* File Variables _______________________________________________*/
char rec[90];
FILE *fp;   /* for .ct file */
char global_post_filename[120];
int global_ct_files; /* number of ct files used for input */
char filename_ct[MAXIMUM_CT_FILES][80]; /* allow at most 10 filenames for input */
char file_data_ct[MAXIMUM_CT_FILES][80]; /* first line of ct file */

char filename[120];
char global_gif_filename[120]; /* this depends on -o or input file */

int global_diag_count[2*MAXIMUM_SIZE];
int global_diag_start[2*MAXIMUM_SIZE];
int global_helices_in_plot_file; /* number of entries in file */

/* _________________________________________________________________*/


/* given a row and column, return the color of that position */

void display_all_helices(void)
{int i;
 for(i=0;i<global_helices_in_plot_file;i++)
   {printf("\n helix %d row =%d column=%d color=%d length=%d diag is %d or %d",
        i,global_diag[i].row,global_diag[i].column,
        global_diag[i].color,global_diag[i].length,
        global_diag[i].row+global_diag[i].column-1,global_diag[i].diagonal);
      }
 for(i=1;i<=2*global_len;i++)
   {printf("\n diag %d starts at %d",i,global_diag_start[i]);
   }
}

/*____________________________________________________________________*/
/* Initialization Routines for nongraphics */


void initialize_labels(void)
{int i;
 global_zoom_labels_count=0; /* initialize to no labels */
 for(i=0;i<50;i++) /* allocate memory for string */
   {global_zoom_labels_string[i]=(char *)malloc(15*sizeof(char));
    if(global_zoom_labels_string[i]==NULL)
      {printf("\n Error------- insufficient memory for zoom labels");
       exit(1);
      }
   }
}
void process_input_files(void)
{global_len=0;
 global_helices_in_plot_file=0;
 read_of_all_ct_files(fp,global_diag_count,global_diag_start,global_ct_files,
                  &global_helices_in_plot_file,&global_len,
                  filename_ct,global_diag,file_data_ct,global_sequence_name,
                  global_sequence_name_set);
 initialize_colors(global_color,global_ct_files);
}



void display_bad_arguments(int end_flag)
{
    printf("\n Bad  arguments:\n");
    printf("\n[ Usage:          ct_boxplot name.ct                       ]");
    printf("\n[             or  ct_boxplot -o name_out name.ct           ]");
    printf("\n[              (.ct suffix is optional)                      ]");
    printf("\n[                                                            ]");
    printf("\n[                Valid arguments                             ]");
    printf("\n[ -a Struc     Structure to display in lower triangle        ]");
    printf("\n[                (-3<=Struc<=number of ct files)             ]");
    printf("\n[                Default is Full and Partial overlap         ]");
    printf("\n[                -3 is Full and partial, -1 is Full overlap, ]");
    printf("\n[                 0 is Partial overlap, -2 is none           ]");
    printf("\n[                For 3 sequences, -11 produces Partial       ]");
    printf("\n[                    Overlap interpreted for content.        ]");
    printf("\n[                    -13 for Partial Overlap interpreted and ]");
    printf("\n[                     Full overlap.                          ]");
    printf("\n[ -b           Enable clear background for GIF output        ]");
    printf("\n[ -d direct/   Look for ct files in directory direct         ]");
    printf("\n[ -f           Force gray as multicolor                      ]");
    printf("\n[ -g           Create a Gif file, default is ps              ]");
    printf("\n[ -i           For 3 sequences, gray is shown as 3 colors    ]");
    printf("\n[                to indicate content                         ]");
    printf("\n[ -j           Draw postscript dots at least 1/72 inch       ]");
    printf("\n[ -jpg         Create jpg instead of gif                     ]");
    printf("\n[ -l i j       Label the point i,j in the GIF or zoomed ps   ]");
    printf("\n[ -m value     Magnify dot size by value. (.2 to 100.)       ]");
    printf("\n[ -o name_out  Specify name of output file as                ]");
    printf("\n[              name_out.ps or name_out.gif                   ]");
    printf("\n[ -pg          Create both GIF and postscript output         ]");
    printf("\n[ -png         Create png instead of gif                     ]");
    printf("\n[ -r  RES      Specify Gif resolution from 20 to 300         ]");
    printf("\n[ -s           Save the tgd file for editing the gif         ]"); 
    printf("\n[ -t \"TITLE\"   Specify TITLE as title, (Quotes Required)     ]");
    printf("\n[ -w           Create .gifdat for webzooming                ]");
    printf("\n[ -y s i j     Zoom with scale s about row i, column j       ]");
    printf("\n[ -z l r t b   Specify a zoom region: left,right,top,bottom  ]");
    printf("\n \n ");
    if(end_flag==TRUE)
      {printf("  -End of session  . .. \n");
       exit(1);
      }
   }
 
void initialize_data(int argc, char **argv)
{ int i,test;
  int gif_resolution;
  int zoom_flag;
  int save_tgd_flag;
  int grid_flag;
  int type_of_output_file;
  int lower_structure;
  int post_file;
  int clear_gif_flag;
  int dot_zoom_flag,dot_zoom_row,dot_zoom_column;
  int dot_zoom_show;
  int gifdat_flag;
  int big_dot_flag;
  int nongray_for_3_flag;
  int forced_stripe_flag;
  char ct_directory[80];
  char temp_file_name[130];
  float dot_zoom_scale;
  int temp;
  int ct_directory_flag;
  char number_rec[25];
  float dot_magnifier;
  int display_l,display_r,display_t,display_b;
   /* For output,2 is postscript, 1 is gif, 3 is both */
  char argument[80];
  big_dot_flag=FALSE;
  nongray_for_3_flag=FALSE;
  forced_stripe_flag=FALSE;  
  gifdat_flag=FALSE; /* do not create gifdat flag for web zooming */
  dot_zoom_flag=FALSE; /* zoom about a dot, default is FALSE */
  grid_flag=TRUE; /* Default is to have a grid */
  clear_gif_flag=FALSE;/* Default is not to use a clear gif */
  save_tgd_flag=FALSE; /* Default is not to save the tgd file */
  global_sequence_name_set=FALSE;
  ct_directory_flag=FALSE;
  zoom_flag=FALSE; /* Default is to display the full image */
  /* This zooms based on left, right, top,bottom */
  gif_resolution=72; /* Default resolution for gif */ 
  i=1;
  lower_structure=-12; /* default is to show both full and partial*/
                       /* overlap in lower triangle */
  /* show interpreted if 3 sequences */
  dot_magnifier=1.0; /* default is no magnification */
  type_of_output_file=2; /* default is postscript */
  global_ct_files=0;
  strcpy(global_sequence_name,"Name not in first ct file?");
  /* set default filenames */
  strcpy(global_post_filename,"ct_boxplot_out");
  strcpy(global_gif_filename,"ct_boxplot_out");
  initialize_labels();
  global_png_mode=FALSE;
  global_jpg_mode=FALSE;
  if(argc<2)
    {display_bad_arguments(TRUE);
    }
  while(i<=(argc-1))
      { /* printf("\n argument %d is %s",i,argv[i]); */
       test=strcmp(argv[i],"-s");
       if(test==0)
         {save_tgd_flag=TRUE; /* Save the tgd file */
          i++;
         }
       else
         {test=strcmp(argv[i],"-d");
          if(test==0)
            {
             if(i<(argc-1))
	      {/* make sure there is an argument to use */
               strcpy(ct_directory,argv[i+1]);
               ct_directory_flag=TRUE;
              }
             i+=2;
            }
       else
         {test=strcmp(argv[i],"-o");
          if(test==0)
            {
             if(i<(argc-1))
	      {/* make sure there is an argument to use */
               strcpy(global_post_filename,argv[i+1]);
               strcpy(global_gif_filename,argv[i+1]);
              }
             i+=2;
            }
       else
	 {test=strcmp(argv[i],"-b");
           if(test==0)
             {clear_gif_flag=TRUE;
              i++;
             }
       else
	 {test=strcmp(argv[i],"-png");
           if(test==0)
             {global_png_mode=TRUE;
              type_of_output_file=1;
              printf("\n caught png mode");
              i++;
             } 
         {test=strcmp(argv[i],"-jpg");
           if(test==0)
             {global_jpg_mode=TRUE;
              type_of_output_file=1;
              printf("\n caught jpg mode");
              i++;
             }
        else
	 {test=strcmp(argv[i],"-i");
           if(test==0)
             {nongray_for_3_flag=TRUE;
              i++;
             }
        else
	 {test=strcmp(argv[i],"-j");
           if(test==0)
             {big_dot_flag=TRUE;
              i++;
             }
         else
	 {test=strcmp(argv[i],"-f");
           if(test==0)
             {forced_stripe_flag=TRUE;
              i++;
             }
       else
	 {test=strcmp(argv[i],"-w");
           if(test==0)
             {gifdat_flag=TRUE;
              i++;
             }
      else
	 {test=strcmp(argv[i],"-a");
           if((test==0)&&(i<(argc-1)))
             {lower_structure=atoi(argv[i+1]);
              if(lower_structure>MAXIMUM_CT_FILES)
                 lower_structure=MAXIMUM_CT_FILES;
              lower_structure=lower_structure+1;
              i+=2;
             }
       else
	 {test=strcmp(argv[i],"-r");
           if((test==0)&&(i<(argc-1)))
             {gif_resolution=atoi(argv[i+1]);
              if(gif_resolution<20)
                 gif_resolution=20;
              if(gif_resolution>300)
                 gif_resolution=300;
              i+=2;
             }
       else {test=strcmp(argv[i],"-m");
	      if((test==0)&&(i<(argc-1)))/* set magnification */
	       {sscanf(argv[i+1],"%f",&dot_magnifier);
		 if(dot_magnifier>100.)
                        dot_magnifier=100.0; 
		 else if(dot_magnifier<0.2)
			 dot_magnifier=0.2;
                i+=2;  
               }
      else
	 {test=strcmp(argv[i],"-g");
           if(test==0)
             {type_of_output_file=1;
              i++;
             }
      else
	 {test=strcmp(argv[i],"-pg");
           if(test==0)
             {type_of_output_file=3;
              i++;
             }
	   else
            {test=strcmp(argv[i],"-y");
             if((test==0)&&(i<(argc-3)))
               {dot_zoom_flag=TRUE;
                sscanf(argv[i+1],"%f",&dot_zoom_scale);
                if(dot_zoom_scale<.01)
                 dot_zoom_scale=.1;
                if(dot_zoom_scale>200.)
                  dot_zoom_scale=200.;
                dot_zoom_row=atoi(argv[i+2]);
                dot_zoom_column=atoi(argv[i+3]);
                if(dot_zoom_row<1)
                 dot_zoom_row=1;
                if(dot_zoom_column<1)
		  dot_zoom_column=1;
                 i+=4;
               }
      else
	 {test=strcmp(argv[i],"-z");
           if((test==0)&&(i<(argc-4)))
             {zoom_flag=TRUE;
              display_l=atoi(argv[i+1]);
	      display_r=atoi(argv[i+2]);  
	      display_t=atoi(argv[i+3]);
	      display_b=atoi(argv[i+4]);
              i+=5;
             }
       else
         {test=strcmp(argv[i],"-l");
           if((test==0)&&(i<(argc-2)))
             {global_zoom_labels_count=1;
              global_zoom_labels_row[0]=atoi(argv[i+1]);
	      global_zoom_labels_col[0]=atoi(argv[i+2]);
              i+=3;
             }
       else
         {test=strcmp(argv[i],"-t");
           if((test==0)&&(i<(argc-1)))
             {global_sequence_name_set=TRUE;
              strcpy(global_sequence_name,argv[i+1]);
              i+=2;
             }
                
       else
          {strcpy(argument,argv[i]);
           if(argument[0]=='-')
               {printf("\n Error!!   %s was not recognized. --------",
                          argument);
                printf("\n Filenames cannot begin with a hyphen");
                printf("\n Follow each flag with enough parameters\n");
                i++;
                display_bad_arguments(FALSE);
                }
       else
          {strcpy(filename_ct[global_ct_files],argv[i]);
           if(global_ct_files<MAXIMUM_CT_FILES)
               global_ct_files++;
           else
              {printf("\n\n  !!!!Error, too many ct files ");
               printf("\n Only %d will be used !!!!!!!!\n\n",MAXIMUM_CT_FILES);
              }
            i++;
           }
	   }
	 }
         }}
	 }
	    }}
         }
         }
         }}}}
         }}
         }}
	 }}
      }
  for(i=0;i<global_ct_files;i++)
    {fix_title(filename_ct[i]);
     strcat(filename_ct[i],".ct");
     if(ct_directory_flag)
       {
        strcpy(temp_file_name,ct_directory);
        strcat(temp_file_name,filename_ct[i]);
        strcpy(filename_ct[i],temp_file_name);
       }
     printf("\n filename %d is %s",i,filename_ct[i]); 
    }
 if(lower_structure<-2)
   {if(global_ct_files!=3)
      lower_structure=-2;
    else
      {if(lower_structure<-10)
        lower_structure=-12;
       else
        lower_structure=-10;
      }
   }
 if(nongray_for_3_flag)
   {if(global_ct_files!=3)
      {printf("\n Error, Must have 3 ct files to use -i option\n");
       nongray_for_3_flag=FALSE;
      }
   }
 if(ct_directory_flag)
       {/* printf("fixing filenames");*/
        strcpy(temp_file_name,ct_directory);
        strcat(temp_file_name,global_post_filename);
        strcpy(global_post_filename,temp_file_name);
        strcpy(temp_file_name,ct_directory);
        strcat(temp_file_name,global_gif_filename);
        strcpy(global_gif_filename,temp_file_name);
        printf("\n gif is %s  postscript is %s",global_gif_filename,
               global_post_filename);
       }
  global_last_structure=global_ct_files+1;
  strcat(global_post_filename,".ps");
  if(global_png_mode)
    strcat(global_gif_filename,".png");
  else
    {if(global_jpg_mode)
       strcat(global_gif_filename,".jpg");
     else
       strcat(global_gif_filename,".gif");
    }
  process_input_files();
  if((nongray_for_3_flag)||lower_structure<=-10)
    {global_gray3_content=store_gray_content(global_gray3_content,
          global_helices_in_plot_file,
          global_len,global_diag,global_diag_start);
    }
  if(dot_zoom_flag)
    {dot_zoom_show=(int)(global_len/dot_zoom_scale);
     if(dot_zoom_show<3)
        dot_zoom_show=3;
     zoom_flag=TRUE;
     display_l=dot_zoom_column-dot_zoom_show/2;
     display_r=dot_zoom_column+dot_zoom_show/2;
     display_t=dot_zoom_row-dot_zoom_show/2;
     display_b=dot_zoom_row+dot_zoom_show/2;
     /* keep zoom region in bounds without altering zoom scale */
     if(display_l<1)
         {display_r=display_r-display_l;
          display_l=1;
         }
     else
       {if(display_r>global_len)
          {display_l=display_l-(display_r-global_len);
           display_r=global_len;
          }
       }
     if(display_t<1)
       {display_b=display_b-display_t;
        display_t=1;
       }
     else
       {if(display_t>global_len)
	 {display_b=display_b-(display_t-global_len);
          display_t=global_len;
         }
       }


    }
  if(zoom_flag)
    {if((display_l>global_len)|| (display_l<1))
       display_l=1;
     if((display_r>global_len)|| (display_r<1))
       display_r=global_len;
     if((display_t>global_len)|| (display_t<1))
       display_t=1;
     if((display_b>global_len)|| (display_b<1))
       display_b=global_len;
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
     display_r=global_len;
     display_t=1;
     display_b=global_len;
    }
 if(global_zoom_labels_count==1) /* Force label in bounds and convert */
                                 /* it to string */
   {if(global_zoom_labels_row[0]<display_t)
       global_zoom_labels_row[0]=display_t;
    else
      if(global_zoom_labels_row[0]>display_b)
       global_zoom_labels_row[0]=display_b;
    if(global_zoom_labels_col[0]<display_l)
       global_zoom_labels_col[0]=display_l;
    else
      if(global_zoom_labels_col[0]>display_r)
       global_zoom_labels_col[0]=display_r;
    strcpy(global_zoom_labels_string[0],"\\(");
    sprintf(number_rec,"%d",global_zoom_labels_row[0]);
    strcat(global_zoom_labels_string[0],number_rec);
    strcat(global_zoom_labels_string[0],",");
    sprintf(number_rec,"%d",global_zoom_labels_col[0]);
    strcat(global_zoom_labels_string[0],number_rec);
    strcat(global_zoom_labels_string[0],"\\)");
   }
 if(type_of_output_file!=2)
   finish_making_gif(TRUE,global_gif_filename,gif_resolution,
      global_zoom_labels_count,global_zoom_labels_row,
       global_zoom_labels_col,grid_flag,global_last_structure,
       global_diag,global_color,save_tgd_flag,display_l,
       display_r,display_t,
          display_b,lower_structure,global_sequence_name,
          dot_magnifier,global_diag_start,global_diag_count,
          file_data_ct,clear_gif_flag,gifdat_flag,global_len,
          forced_stripe_flag,
          nongray_for_3_flag,global_gray3_content,global_png_mode,
          global_jpg_mode);
  if(type_of_output_file!=1)
    {if(!zoom_flag)
      post_file=MAIN_POST; /* create post for main */
     else
      post_file=ZOOM_POST; /* create post for zoom */
     finish_making_post(big_dot_flag,post_file,global_post_filename,
               global_zoom_labels_count,display_l,
               display_r,display_t,display_b,global_color,
               global_zoom_labels_row,global_zoom_labels_col,
               global_zoom_labels_string,grid_flag,
               global_last_structure,global_diag_start,
               global_diag_count,global_diag,
               global_len,lower_structure,global_sequence_name,
               file_data_ct,dot_magnifier,forced_stripe_flag,
               nongray_for_3_flag,global_gray3_content);
    }

}


/*_______________________________________________________________________*/


/*____________________________________________________________________*/
/*_______________________________________________________________*/

 
void try_exit(int exit_choice)
{ if(exit_choice==1)
   { printf("\n Normal exit from ct_boxplot_ng \n\n");
     exit(0);
   }
}


/*________________________________________________________________*/
/* Plotting functions _____________________*/


int main(int argc, char **argv) 
{ printf("\n ct_boxplot_ng Version 2.0");
  initialize_data(argc,argv);
  try_exit(1);
  return 0;     /* the program never hits this line */
}





