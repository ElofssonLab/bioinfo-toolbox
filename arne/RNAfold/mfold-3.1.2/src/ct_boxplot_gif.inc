/* Darrin Stewart and Michael Zuker */
/* Feb 4, 1999  */
/* ct_boxplot_gif.c */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "copyright_gif.inc"


FILE *giffp;
FILE *gifdatfp;

const char *gif_color[ COLOR_TABLE_SIZE]={"c0","c1","c2","c3","c4","c5","c6",
             "c7","c8","c9","c10","c11","c12","c13","c14","c15","c16",
             "c17","c18","c19","c20",
 "c21","c22","c23","c24","c25","c26","c27","c28","c29","c30",
 "c31","c32","c33","c34","c35","c36","c37","c38","c39","c40",
 "c41","c42","c43","c44","c45","c46","c47","c48","c49","c50",
 "c51","c52","c53","c54","c55","c56","c57","c58","c59","c60",
 "c61","c62","c63","c64","c65","c66","c67","c68","c69","c70",
 "c71","c72","c73","c74","c75","c76","c77","c78","c79","c80",
 "c81","c82","c83","c84","c85","c86","c87","c88","c89","c90",
 "c91","c92","c93","c94","c95","c96","c97","c98","c99","c100",
  "c101","c102","c103","c104","c105","c106",
  "c107","c108","c109","c110","c111","c112","c113","c114","c115","c116",
  "c117","c118","c119","c120",
 "c121","c122","c123","c124","c125","c126","c127","c128","c129","c130",
 "c131","c132","c133","c134","c135","c136","c137","c138","c139","c140",
 "c141","c142","c143","c144","c145","c146","c147","c148","c149","c150",
 "c151","c152","c153","c154","c155","c156","c157","c158","c159","c160",
 "c161","c162","c163","c164","c165","c166","c167","c168","c169","c170",
 "c171","c172","c173","c174","c175","c176","c177","c178","c179","c180",
 "c181","c182","c183","c184","c185","c186","c187","c188","c189","c190",
 "c191","c192","c193","c194","c195","c196","c197","c198","c199","c100",
  "c201",
					  "ctext","cbak","cd","cgrid"
           };

/* create gif output for ct_boxplot */

int open_gif(char *gif_tempname)
{    /* open specified output file */
   if ((giffp = fopen(gif_tempname, "w")) == NULL) /* open a file */
         {printf ("\n * Could not open file:  %s", gif_tempname);
	  printf("\n This file is a temporary file used to create a ");
          printf("\n gif file.\n");
          return(1);
          }
    return (0);
}/*___________________________________________________________________________*/
/* create a gif */
/*___________________________________________________________________________*/
/* convert a row to a position in pixels down from the top */
int gif_row_to_y(int dis_t,int y_start,float y_pic_per_dot,int row)
{ 
   return y_start+(int)(y_pic_per_dot*((float)(row-dis_t))+.500001);
  
}
/* convert a column to position in pixels right from the left edge */
int gif_col_to_x(int dis_l,int x_start,float x_pic_per_dot,int column)
{
   return x_start+(int)(x_pic_per_dot*((float)(column-dis_l))+.50001);
 
}/* draw a single dot */
void gif_plot_dot(int x_start,int y_start, float x_pic_per_dot,
                 float y_pic_per_dot,int dis_l,int dis_t,
                 char *color_string,int color,int row,int column,
                 float dot_magnifier,
                  int stripe_flag,int last_structure)
{  int x_pos;/* positions for center of dot */
   int y_pos;
   int ix_dot_size,iy_dot_size;
   int upper_left_x,upper_left_y,lower_right_x,lower_right_y;
   float adjuster;
   float scale_adjust_x,scale_adjust_y;/* make large dots smaller */
   /* make small dots larger */
   int lx,rx,ty,by; /* positions for corners of dot */
   x_pos=gif_col_to_x(dis_l,x_start,x_pic_per_dot,column);
   y_pos=gif_row_to_y(dis_t,y_start,y_pic_per_dot,row);
   /* draw a rectangle here */
     {scale_adjust_x=dot_magnifier;
      scale_adjust_y=dot_magnifier; /* make small dots bigger */
      ix_dot_size=(int)(x_pic_per_dot*scale_adjust_x+.7);
      iy_dot_size=(int)(y_pic_per_dot*scale_adjust_y+.7);
      if(ix_dot_size<2)
        ix_dot_size=2;
      if(iy_dot_size<2)
        iy_dot_size=2;
      lx=(int)(x_pos-ix_dot_size/2+.5);
      rx=(int)(lx+ix_dot_size-1+.5);
      ty=(int)(y_pos-iy_dot_size/2+.5);
      by=(int)(ty+iy_dot_size-1+.5);
      if(!stripe_flag)
        fprintf(giffp,"filledrectangle im %d %d %d %d %s\n",
         lx,ty,rx,by,color_string);
      else
	{adjuster=(color-1)*1./(float)(last_structure-1);
         upper_left_x=(int)(x_pos-adjuster*ix_dot_size/2+.5);
         upper_left_y=(int)(y_pos-adjuster*iy_dot_size/2+.5);
         lower_right_x=(int)(x_pos+adjuster*ix_dot_size/2+.5-1.);
         lower_right_y=(int)(y_pos+adjuster*iy_dot_size/2+.5-1.);
         fprintf(giffp,"filledpolygon im %d %d %d %d %d %d %d %d 4 %s\n",
             upper_left_x,upper_left_y,rx,ty,
             lower_right_x,lower_right_y,lx,by,color_string);    

        }
     } 
}

/* create the tic marks */
void gif_single_hor_tic(int x_start,int j_col,float x_pic_per_dot,int ty,
                        int by,int dis_l,int grid_flag)
{ int x_pos;
  x_pos=gif_col_to_x(dis_l,x_start,x_pic_per_dot,j_col);
  if(grid_flag==TRUE)
    fprintf(giffp,"line im %d %d %d %d %s\n", x_pos,ty,x_pos,by-2,
	    gif_color[COL_GRID]);
  fprintf(giffp,"line im %d %d %d %d %s\n",x_pos,ty,x_pos,ty-15,
          gif_color[ COL_TEXT]);
  fprintf(giffp,"string im gdFont10x20 %d %d \"%d\" %s\n",
                         x_pos-10,ty-36,j_col,gif_color[COL_TEXT]);
}

void gif_plot_hor_tics(int x_start,float x_pic_per_dot,int dis_l,
                       int dis_r,int ty,int by,int grid_flag)
{ 
  int j_col,step;
  int start;
  step=step_fun(dis_l,dis_r);
  start=start_fun(dis_l,step);
  gif_single_hor_tic(x_start,dis_l,x_pic_per_dot,ty,by,dis_l,FALSE);
  gif_single_hor_tic(x_start,dis_r,x_pic_per_dot,ty,by,dis_l,FALSE);
  for(j_col=start;j_col<(dis_r-2*step/3);
                        j_col=j_col+step)
      {gif_single_hor_tic(x_start,j_col,x_pic_per_dot,ty,by,dis_l,
              grid_flag);
      }
}
void gif_single_ver_tic(int y_start,int i_row,float y_pic_per_dot,int lx,int rx,
                        int dis_t,int grid_flag)
{ int y_pos; 
  y_pos=gif_row_to_y(dis_t,y_start,y_pic_per_dot,i_row);
  if(grid_flag==TRUE)
    fprintf(giffp,"line im %d %d %d %d %s\n",rx-2,y_pos,lx+2,y_pos,
	    gif_color[COL_GRID]);
  fprintf(giffp,"line im %d %d %d %d %s\n",rx,y_pos,rx+15,y_pos,
            gif_color[COL_TEXT]);
  fprintf(giffp,"string im gdFont10x20 %d %d \"%d\" %s\n",rx+18,
                          y_pos-12,i_row,gif_color[COL_TEXT]);
}
void gif_plot_ver_tics(int y_start,float y_pic_per_dot,int dis_t,
                       int dis_b,int lx,int rx,int grid_flag)
{ int i_row,step,start;
  step=step_fun(dis_t,dis_b);
  start=start_fun(dis_t,step);
  gif_single_ver_tic(y_start,dis_t,y_pic_per_dot,lx,rx,dis_t,FALSE);
  gif_single_ver_tic(y_start,dis_b,y_pic_per_dot,lx,rx,dis_t,FALSE);
  for(i_row=start;i_row<(dis_b-step/2);i_row=i_row+step)
   {gif_single_ver_tic(y_start,i_row,y_pic_per_dot,lx,rx,dis_t,
                 grid_flag);  
   }
}

void  gif_tics(int x_start,int y_start,float x_pic_per_dot,float y_pic_per_dot,
               int dis_l,int dis_r,int dis_t,int dis_b,int lx,
               int rx,int ty,int by,int grid_flag)
{  gif_plot_ver_tics(y_start,y_pic_per_dot,dis_t,dis_b,lx,rx,grid_flag);
   gif_plot_hor_tics(x_start,x_pic_per_dot,dis_l,dis_r,ty,by,grid_flag);
}
void  gif_points(int counter[MAXIMUM_CT_FILES+1],int stripe_flag,
                int x_start,int y_start,float x_pic_per_dot,
                float y_pic_per_dot,int dis_l,int dis_r,
                  int dis_t,int dis_b,int last_structure,
                  struct helix *diag,
                  int lower_structure,int *diag_start,int *diag_count,
                  float dot_magnifier,int nongray_for_3_flag,
                  int *gray3_content,int gray3[3])
{
 int row,column,n_column,new_row;
 int location,new_color;
 int color,starting_color,ending_color;
 int diagonal,start_diag,end_diag,position,count;
 int old_color;
 char color_string[4];
 int gray3_column;
 start_diag=dis_t+dis_l-1;
 end_diag=dis_b+dis_r-1;
 /* single structures */
  gray3[0]=0;
  gray3[1]=0;
  gray3[2]=0;
 for(color=last_structure;color>=2;color--)
   {strcpy(color_string,gif_color[color]);
    for(diagonal=start_diag;diagonal<=end_diag;diagonal++)
      {position=diag_start[diagonal];
       for(count=0;((count<diag_count[diagonal])&&
           (diag[position+count].row<=dis_b) ) ;count++)
	 {location=position+count;
          new_color=diag[location].color;
              if(new_color==color)
	      { new_row=diag[location].row;
                column=diag[location].column;
                for(row=new_row;
                      row<(new_row+diag[location].length);row++)
	           {if((row>=dis_t)&&(row<=dis_b)
                          &&(column>=dis_l)&&(column<= dis_r))
                          { gif_plot_dot(x_start,y_start,x_pic_per_dot,
                                         y_pic_per_dot,dis_l,dis_t,
                                         color_string,color,row,column,
                                          dot_magnifier,stripe_flag,
                                          last_structure);
                            counter[color]++;  
                           }
                             column--;
		   }
              }                        

       }    
     }
   }

 /* black or gray dots above main diagonal */

 if((last_structure>3)&&(!stripe_flag))
    starting_color=1;
   else
    starting_color=0;
   for(color=starting_color;color>=0;color--)
   {strcpy(color_string,gif_color[color]);
    if((nongray_for_3_flag)&&(color==1))
       strcpy(color_string,gif_color[5]);
    for(diagonal=start_diag;diagonal<=end_diag;diagonal++)
      {position=diag_start[diagonal];
       for(count=0;((count<diag_count[diagonal])&&
           (diag[position+count].row<=dis_b) ) ;count++)
	 {location=position+count;
          new_color=diag[location].color;
              if(new_color==color)
	      { new_row=diag[location].row;
                column=diag[location].column;
                for(row=new_row;
                      row<(new_row+diag[location].length);row++)
	           {if((row>=dis_t)&&(row<=dis_b)
                          &&(column>=dis_l)&&(column<= dis_r))
                          { gif_plot_dot(x_start,y_start,x_pic_per_dot,
                                         y_pic_per_dot,dis_l,dis_t,
                                         color_string,color,row,column,
                                         dot_magnifier,FALSE,last_structure); 
                           counter[color]++;
                           }
                             column--;
		   }
		 
            }                      

       }    
     }
   }
 /* plot partial overlap for 3 sequences for gray interpreted */

 if(nongray_for_3_flag)
   {for(row=dis_t;row<=dis_b;row++)
     {column=gray3_content[row];
      if(column!=0)
       {if(column<0)
          n_column=-1*column;
        else
          n_column=column;
        if((n_column>=dis_l)&&(n_column<=dis_r))
	   { 
            if(column>0)
             {/* for seq 2 nd 3 */
              strcpy(color_string,gif_color[7]);
              gray3[2]++;
             }
            else
             {/* for seq 1 & 3 */
              strcpy(color_string,gif_color[6]); 
              gray3[1]++;           
             }
            gif_plot_dot(x_start,y_start,x_pic_per_dot,
                                         y_pic_per_dot,dis_l,dis_t,
                                         color_string,1,row,n_column,
                                          dot_magnifier,FALSE,
                                          last_structure);
           }
	       
       }
     }
   gray3[0]=counter[1]-gray3[1]-gray3[2];
 
   }
 /* subdiagonal */
 /* points below the main diagonal */
 /* draw dots on lower diagonal */

 if(lower_structure!=-1) /* -1 is no dots below */
  {if(lower_structure!=-10)
   {if(lower_structure==-2) /* -2 is black and gray dots below */
    {starting_color=1;
     ending_color=0;
    }
  else
   {if(lower_structure==-12)
     {starting_color=0;
      ending_color=0;
     }
   else
    {starting_color=lower_structure; /* a single structure color dot below */
     ending_color=lower_structure;
    }
   }
    for(color=starting_color;color>=ending_color;color--)
     {/* printf("\n current color is %d",color); */
      strcpy(color_string,gif_color[color]);
     for(diagonal=start_diag;diagonal<=end_diag;diagonal++)
      {position=diag_start[diagonal];
       for(count=0;((count<diag_count[diagonal])&&
           (diag[position+count].row<=dis_r) ) ;count++)
	 {location=position+count;
          new_color=diag[location].color;
              if(new_color==color)
	      { new_row=diag[location].row;
                column=diag[location].column;
                if(column>=dis_t)
		  {
                   for(row=new_row;
                      row<(new_row+diag[location].length);row++)
	           {if((column>=dis_t)&&(column<=dis_b)
                          &&(row>=dis_l)&&(row<= dis_r))
                          { gif_plot_dot(x_start,y_start,x_pic_per_dot,
                                         y_pic_per_dot,dis_l,dis_t,
                                         color_string,color,column,row,
                                          dot_magnifier,FALSE,last_structure);
                           }
                             column--;
		   }
                  }
              }                        

       }    
     }
   }
   }
 if(lower_structure<=-10)
  {color=1; /* take care of gray interpreted in lower triangle */
   for(diagonal=start_diag;diagonal<=end_diag;diagonal++)
      {position=diag_start[diagonal];
       for(count=0;((count<diag_count[diagonal])&&
           (diag[position+count].row<=dis_r) ) ;count++)
	 {location=position+count;
          new_color=diag[location].color;
              if(new_color==color)
	      { new_row=diag[location].row;
                column=diag[location].column;
                if(column>=dis_t)
		  {
                   for(row=new_row;
                      row<(new_row+diag[location].length);row++)
	           {if((column>=dis_t)&&(column<=dis_b)
                          &&(row>=dis_l)&&(row<= dis_r))
                          {gray3_column=gray3_content[row]; 
                           if(gray3_column==0)
                             {strcpy(color_string,gif_color[5]);
                              old_color=5;
                             }
                           else
                             {if(gray3_column<0)
                               {strcpy(color_string,gif_color[6]);
                                old_color=6;
                               }
                             else
                              {
                               strcpy(color_string,gif_color[7]);
                               old_color=7;
                              }
                             }
                             if(!nongray_for_3_flag)
                                gray3[old_color-5]++;
                            gif_plot_dot(x_start,y_start,x_pic_per_dot,
                                         y_pic_per_dot,dis_l,dis_t,
                                         color_string,color,column,row,
                                          dot_magnifier,FALSE,last_structure);
                           }
                             column--;
		   }
                  }
              }                        

       }    
     }
  }
 }
}

/* make diagonal line */
void gif_make_diag(int display_l,int display_r,int display_t,int display_b,
 int  x_start,int start_y,float x_pic_per_dot,float y_pic_per_dot)
{int ty,rx;
 int lx,by;
 int right,top,bottom,left;
  if((display_r>=display_t)&&(display_b>=display_l))
    {
     top=display_t;
     bottom=display_b;
     right=display_r;
     left=display_l;
      if(display_r>=display_b)
        {bottom=display_b;
         right=display_b;
        }
      else
        {right=display_r;
         bottom=display_r;
        }
      if(display_l>=display_t)
       {left=display_l;
        top=display_l;
       }
      else
       {top=display_t;
        left=top;
       }   

    ty=gif_row_to_y(display_t,start_y,y_pic_per_dot,top)-
              (int)(y_pic_per_dot/2.+.51);
    lx=gif_row_to_y(display_l,x_start,x_pic_per_dot,left)-
              (int)(x_pic_per_dot/2.+.51);
    rx=gif_col_to_x(display_l,x_start,x_pic_per_dot,right)+
               (int)(x_pic_per_dot/2.+.51);
    by=gif_row_to_y(display_t,start_y,y_pic_per_dot,bottom)+
               (int)(y_pic_per_dot/2.+.51);
    fprintf(giffp,"line im %d %d %d %d gdBrushed \n",rx-1,by,lx-1,ty);
      /* The above draws a wide diagonal line */
   }
} 

void gif_label_points(int x_start,int y_start,float x_pic_per_dot,
                float y_pic_per_dot,int dis_l,
                  int dis_t,int zoom_labels_count,
                  int *zoom_labels_row,int *zoom_labels_col
                  )
{int i,hor_pos,ver_pos,label_hor_pos,lx,ty,rx,by;
 int hor_shift,ver_shift;
 char string[80];
 char temp_string[50];
 hor_shift=(int)(x_pic_per_dot/3);
 ver_shift=(int)(y_pic_per_dot/3);
 for(i=0;i<zoom_labels_count;i++)
   {hor_pos=gif_col_to_x(dis_l,x_start,x_pic_per_dot,
                  zoom_labels_col[i]);
    ver_pos=gif_row_to_y(dis_t,y_start,y_pic_per_dot,
                  zoom_labels_row[i]);
    label_hor_pos=gif_col_to_x(dis_l,x_start,x_pic_per_dot,
                  zoom_labels_col[i]+1);
    lx=hor_pos-hor_shift;
    rx=hor_pos+hor_shift;
    ty=ver_pos-ver_shift;
    by=ver_pos+ver_shift;
    strcpy(string,num_string_int(zoom_labels_row[i],temp_string));
    strcat(string,",");
    strcat(string,num_string_int(zoom_labels_col[i],temp_string));
    fprintf(giffp,"filledrectangle im %d %d %d %d %s\n",lx,ty,rx,by,
            gif_color[COL_DOTS]);
    fprintf(giffp,"string im gdFontSmall %d %d \"(%s) \"%s\n",
          label_hor_pos+1,ver_pos-5,
          string, gif_color[COL_DOTS]);
   }
			
}

void gif_create_colors(float color_table[][3],int last_structure,
         int nongray_for_3_flag,int lower_structure,int jpg_mode)
{int color[COLOR_TABLE_SIZE][3];
 int i,j;
  /* convert colors to integers */
  if((nongray_for_3_flag)||(lower_structure<=-10))
       last_structure=last_structure+3;
  for(j=0;j<COLOR_TABLE_SIZE;j++)
    {for(i=0;i<=2;i++)
       color[j][i]=(int)(255.*color_table[j][i]);
    }
  fprintf(giffp,"create bkim 2 2\n");

  fprintf(giffp,"colorallocate bkim %s %d %d %d\n",gif_color[COL_BACKGROUND],
    color[COL_BACKGROUND][0],color[COL_BACKGROUND][1],
          color[COL_BACKGROUND][2]);
  fprintf(giffp,"colorallocate bkim %s %d %d %d\n",gif_color[COL_TEXT],
    color[0][0],color[0][1],color[0][2]);
  fprintf(giffp,"line bkim 0 0 1 0 %s\n",gif_color[COL_TEXT]);
  fprintf(giffp,"line bkim 0 1 1 1 %s\n",gif_color[COL_TEXT]); 

  /* end of brush definition */
  fprintf(giffp,"colorallocate im %s %d %d %d\n",gif_color[COL_BACKGROUND],
       color[COL_BACKGROUND][0],color[COL_BACKGROUND][1],
       color[COL_BACKGROUND][2]);
  /* first is background color */
   for(j=0;j<=last_structure;j++)
     fprintf(giffp,"colorallocate im %s %d %d %d\n",gif_color[j],
          color[j][0], color[j][1],color[j][2]);
   for(j=202;j<COLOR_TABLE_SIZE;j++)
       fprintf(giffp,"colorallocate im %s %d %d %d\n",gif_color[j],
          color[j][0], color[j][1],color[j][2]);
               /* Gray for grid lines */ 
  fprintf(giffp,"setbrush im bkim\n");/* Sets the brush  for wide lines*/
  if(!jpg_mode)
  fprintf(giffp,"interlace im 1\n"); /* $$$$$ turn on interlace */
  }
/* Make the gif */


int open_gifdat(char *filename)
{    /* open specified output file */
       /*  filename[strlen(filename) - 1] = NULL; */
   if ((gifdatfp = fopen(filename, "w")) == NULL) /* open a file */
         {printf ("\n * Could not open gifdat file:  %s", filename);
          return(1);
          }
    return (0);
}
void make_gifdat(char *filename,int resolution,int x_start,int y_start,
      float x_pic_per_dot,float y_pic_per_dot,int dis_l,int dis_t,
      int max_len)
{int error;
 int len;
 char gifdat_name[120];
 strcpy(gifdat_name,filename);
 len=strlen(gifdat_name);
/* strip off .gif from name and replace with .gifdat */
 if(len>4)
   {gifdat_name[len-4]='\0';
   }
 strcat(gifdat_name,".gifdat");
 printf("\n The name for the GIFDAT file is %s \n",gifdat_name);
 error=open_gifdat(gifdat_name);
 if(error)
   {printf("Error creating .gifdat file name was %s",gifdat_name);
    return;
   }
 fprintf(gifdatfp,"%d %d %d %f %f %d %d %d \n",resolution,x_start,y_start,
                x_pic_per_dot,y_pic_per_dot,max_len,dis_l,dis_t);
 fclose(gifdatfp); 
}


void finish_making_gif(int stripe_flag,char *gif_filename,int resolution,
           int zoom_labels_count, 
             int *zoom_labels_row,int *zoom_labels_col,int grid_flag,
             int last_structure,struct helix *diag,float color_table[][3],
              int save_tgd_flag,
               int display_l,int display_r,int display_t,int display_b,
                int lower_structure,char *sequence_name,float dot_magnifier,
                int *diag_start,int *diag_count, char file_data_ct[][80],
                int clear_flag,int gifdat_flag, int max_len,
                int forced_stripe_flag,
                int nongray_for_3_flag,int *gray3_content,int png_mode,
                int jpg_mode)
{ char action[80]; /* used for system call */
  int error,error_occurred; /* resolution is 72, 110, or 200 ? */
  char gif_tempname[80]; /* this is name.tgd */
  int dis_l,dis_r,dis_t,dis_b;/* local versions of display_l,display_r */
  char string[80];
  int name_len;
  int center_position;
  int hor_pos,ver_pos;
  int color,offset;
  int other_stop_point;
  int counter[MAXIMUM_CT_FILES+1];
  int stop_point;
  /* use 8.5 wide by 11 inches wide gif ? */
  int gif_width; /* width of image in pixels */
  int gif_height; /* height of image in pixels */
  int lx,ty,rx,by; /* top,bottom,right,left corn of rec */
  int i;
  int gray3[3];
  float x_pic_per_dot,y_pic_per_dot;
  int x_start,y_start;
  for(i=0;i<=last_structure;i++)
     counter[i]=0;
  dis_l=display_l;
  dis_r=display_r;
  dis_t=display_t;
  dis_b=display_b;
  gif_width=(int)(8.5*resolution)-1;
  gif_height=11*resolution-1; 
  strcpy(gif_tempname,gif_filename);
  name_len=strlen(gif_tempname);
  if(name_len>4)
     gif_tempname[name_len-4]='\0';
  strcat(gif_tempname,".tgd");
  error=open_gif(gif_tempname);
  if(error)
    {if(png_mode)
      printf("\n error creating PNG file ");
     else
      {if(jpg_mode)
         printf("\n error creating GIF file "); 
       else
         printf("\n error creating GIF file ");
      }
     return;
    }                                 /* hor by ver or x by y */
  fprintf(giffp,"create im %d %d\n",gif_width+1,gif_height+1);
                     /* declare image size in pixels */
  if(png_mode)
    printf("\n Creating png ");
  else
    {if(jpg_mode)
        printf("\n creating jpg ");
     else
       printf("\n Creating gif ");
    } 
  printf("\n  %d by %d \n",gif_width+1,gif_height+1);
  gif_create_colors(color_table,last_structure,nongray_for_3_flag,
       lower_structure,jpg_mode);
   /* declare gif colors */
  if(clear_flag==TRUE)
    fprintf(giffp,"colortransparent im %s\n",gif_color[COL_BACKGROUND]);
  /* make box */ /* 7 by 7 inches */
  ty=(int)(1.57*resolution); /* come down 3 inches */
  by=(int)(gif_height-2.37*resolution);
  lx=(int)(.47*resolution); /* left side of box */ 
  rx=gif_width-(.95*resolution); /* 6.5 inches wide drawing region*/
  /* an extra .1 inch was used to make the line outside the region */
  /* Using a wide line too */
  display_copyright_message_gif(giffp,"ct_boxplot");
  fprintf(giffp,"line im %d %d %d %d gdBrushed\n",lx,ty,rx,ty);
  fprintf(giffp,"line im %d %d %d %d gdBrushed\n",lx,by,rx,by);
  fprintf(giffp,"line im %d %d %d %d gdBrushed\n",lx,ty,lx,by);
  fprintf(giffp,"line im %d %d %d %d gdBrushed\n",rx,ty,rx,by);
  center_position=4*resolution-(22+strlen(sequence_name))*6;
  if(center_position<0)
     center_position=0;
  fprintf(giffp,
           "string im gdFont12x24 %d %d \"Structure dot plot for %s\"  %s\n",
              center_position,(int)(resolution/2),
              sequence_name,gif_color[COL_TEXT]);
  x_pic_per_dot=resolution*7.0/(dis_r-dis_l+1);
   /* width of basepair*/
  y_pic_per_dot=resolution*7.0/(dis_b-dis_t+1);/* height of basepair*/
  x_start=(int)(.5*resolution+x_pic_per_dot/2+.5);
   /* distance from left edge for first basepair*/
  y_start=(int)(1.59*resolution+y_pic_per_dot/2+.5);
  if(gifdat_flag==TRUE)
     {make_gifdat(gif_filename,resolution,x_start,y_start,x_pic_per_dot,
         y_pic_per_dot,dis_l,dis_t,max_len);
     }
     /* distance from top for first basepair*/
  /* both above are in pixels */
  /* create the tic marks */
  /*                      */
  gif_tics(x_start,y_start,x_pic_per_dot,y_pic_per_dot,dis_l,
                    dis_r,dis_t,dis_b,lx,rx,ty,by,grid_flag);
 /* make the diagonal line */
  gif_make_diag(dis_l,dis_r,dis_t,dis_b,x_start,
              y_start,x_pic_per_dot,y_pic_per_dot);
  /* create the points */
 if(forced_stripe_flag)
    stripe_flag=TRUE;
 else
  {
    if((nongray_for_3_flag)||(last_structure<4))
       stripe_flag=FALSE;
    if(stripe_flag==TRUE)
     {
      if( !( (x_pic_per_dot>(STRIPE_GIF*(last_structure-1)))&&
          (y_pic_per_dot>(STRIPE_GIF*(last_structure-1))) ) )
       stripe_flag=FALSE;
     }
  }
  gif_points(counter,stripe_flag,x_start,y_start,x_pic_per_dot,y_pic_per_dot,
         dis_l,dis_r,dis_t,dis_b,last_structure,diag,
         lower_structure,diag_start,diag_count,dot_magnifier,
         nongray_for_3_flag,gray3_content,gray3);
  /* Display colors */
  hor_pos=resolution/2;
  ver_pos=gif_height-2.2*resolution;
  fprintf(giffp,"string im gdFontLarge %d %d \"Overlap\" %s\n",
       hor_pos,ver_pos,gif_color[COL_TEXT]);
  fprintf(giffp,"string im gdFontLarge %d %d \"Full:   %4d \" %s\n",
       hor_pos,ver_pos+18,counter[0],gif_color[COL_COMPLETE_OVERLAP]);
  if(!stripe_flag)
   {
    if(last_structure>3)
       {fprintf(giffp,"string im gdFontLarge %d %d \"Partial: %3d \" %s\n",
        hor_pos,ver_pos+36,counter[1],gif_color[COL_PART_OVERLAP]);
       } 
   } 
   stop_point=last_structure/2+1;
   if(last_structure>16)
    { stop_point=9;
      other_stop_point=15; 
      hor_pos=resolution*6-10;
      fprintf(giffp,"string im gdFontLarge %d %d \"%d others \" c16\n",
          hor_pos,(int)(gif_height-.47*resolution+2),last_structure-15);
     }
  else
    other_stop_point=last_structure;
   hor_pos=resolution*3-50; 
   for(color=2,offset=0;color<=stop_point;
             color++,offset=offset+18)
         {fprintf(giffp,"string im gdFontLarge %d %d \"%s %3d\" %s\n",
          hor_pos,ver_pos+offset,file_data_ct[color-2] ,counter[color],gif_color[color]);
          }
   hor_pos=resolution*6-30; 
   for(color=(stop_point+1),offset=0;color<=other_stop_point;
             color++,offset=offset+18)
         {fprintf(giffp,"string im gdFontLarge %d %d \"%s %3d\" %s\n",
          hor_pos,ver_pos+offset,file_data_ct[color-2] ,counter[color],gif_color[color]);
         }
   if((nongray_for_3_flag)||(lower_structure<=-10))
    {fprintf(giffp,"string im gdFontLarge %d %d \"1 & 2    %4d\" %s\n",hor_pos,
            ver_pos+38,gray3[0],gif_color[5]);
    fprintf(giffp,"string im gdFontLarge %d %d \"1 & 3    %4d\" %s\n",hor_pos,
            ver_pos+56,gray3[1],gif_color[6]);
    fprintf(giffp,"string im gdFontLarge %d %d \"2 & 3    %4d\" %s\n",hor_pos,
            ver_pos+74,gray3[2],gif_color[7]);
    }
   /* Display labels */
   if(zoom_labels_count>0)
     gif_label_points(x_start,y_start,x_pic_per_dot,y_pic_per_dot,dis_l,
             dis_t,zoom_labels_count, zoom_labels_row,
                   zoom_labels_col);
   strcpy(string,"gif im ");
  strcat(string,gif_filename); 
  strcat(string,"\n");
  fprintf(giffp,string);
  fclose(giffp);
  if(png_mode)
    printf(" The name for the PNG file is       %s\n",gif_filename);
  else
    {if(jpg_mode)
        printf(" The name for the JPG file is       %s\n",gif_filename);
     else
      printf(" The name for the GIF file is       %s\n",gif_filename); 
    }
  /* convert tgd file to gif */
  if(png_mode)
    strcpy(action,"tgd_png ");
  else
    {if(jpg_mode)
      strcpy(action,"tgd_jpg ");
     else
      strcpy(action,"tgd ");
    }
  strcat(action,gif_tempname);
  error_occurred=system(action);
  if(error_occurred!=0)
    {printf("\n Error converting file to GIF");
     printf("\n Perhaps tgd does not exist");
     printf("\n %s should be deleted",gif_tempname);
    }
  else
    {if(save_tgd_flag!=TRUE)
      {strcpy(action,"rm ");
       strcat(action,gif_tempname);
       system(action);
      }
     else
        {printf("\n Edit %s",gif_tempname);
         printf("\n and run tgd on %s",gif_tempname);
         printf("\n to customize the gif.");
         printf("\n Please delete %s when finished.\n",gif_tempname);
         }
    }

}  


