/* ct_boxplot.h  Feb 4, 1999*/
/* Darrin Stewart and Michael Zuker */

#define MAXIMUM_HELICES 200000 /* This may need to be adjusted */
/* It is important that the number of helices in the plot file does */
/* not exceed this number */

#define MAXIMUM_CT_FILES 200 /* Only this many ct files can be read in */
/* Changing the above number is not advised */

#define MAXIMUM_SIZE 20000 /* maximum sequence length */
#define TRUE    1
#define FALSE   0

#define MAIN_POST    1 /* create postcript for main or zoom */
#define ZOOM_POST    2 /* used in general_post, create_post */

#define COLOR_TABLE_SIZE MAXIMUM_CT_FILES+6
#define COL_PART_OVERLAP  1
#define COL_COMPLETE_OVERLAP 0
#define COL_TEXT      MAXIMUM_CT_FILES+2
#define COL_BACKGROUND    MAXIMUM_CT_FILES+3
#define COL_DOTS    MAXIMUM_CT_FILES+4
#define COL_GRID    MAXIMUM_CT_FILES+5

#define STRIPE_SCREEN 4.0  /* 4.2 is Good for the screen */
#define STRIPE_PS     2.8  /* 5.0 is good for postscript */
#define STRIPE_GIF    2.8  /* 6.0 is good for the GIF */

/* When the dot switches from black or gray for overlap to multiple */
/* Colors depends on whether there are stripe_* pixels,points for each */
/* Color. Use smaller values to make switching occur with smaller dots */

struct helix {
              int row;
              int column;
              int diagonal;
              int color;         
              int length;
             };

void read_of_all_ct_files(FILE *,int *,int *,int,int *,int *,
                 char filename_ct[MAXIMUM_CT_FILES][80],struct helix *,
                 char file_data_ct[MAXIMUM_CT_FILES][80],char *,int);
 

void fix_title(char *);

int step_fun(int,int);

int start_fun(int,int);

void finish_making_post(int,int, char *,int,int,int,int,int, 
         float color_table[][3],int *,int *,char *string[50],
         int,int,int*,int*,struct helix *,int,int,char *,
         char file_data_ct[][80],float ,int,int,int *);

void finish_making_gif(int,char *,int,int,
                       int*, int*, int,
                       int, struct helix *,float table[][3],
                       int,int,int,int,int,int,char *,float,
                       int *,int *, char file_data_ct[][80],int,int,
                       int,int,int,int *,int,int);

char *num_string_int(int,char *);

void initialize_colors(float color_table[][3],int global_ct_files);

