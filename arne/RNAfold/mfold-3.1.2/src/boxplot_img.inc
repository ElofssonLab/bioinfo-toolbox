/* boxplot_ps.inc */
/* by Darrin Stewart and Michael Zuker */
/* May 17, 2000 */



FILE *giffp; /* for gif file */

void clean_title_gif(char *title)
{ int i;
  int length;
  length=strlen(title);
  /* replace double with single quotes */
  for(i=0;i<length;i++)
     {if(title[i]==34)
        title[i]=39;
     }
}
/*___________________________________________________________________________*/
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
 
}


/* draw a single dot */
void gif_plot_dot(int x_start,int y_start, float x_pic_per_dot,
                 float y_pic_per_dot,int dis_l,int dis_t,
                 char *color,int row,int column,int energy)
{  int x_pos;/* positions for center of dot */
   int y_pos;
   int ix_dot_size,iy_dot_size;
   float scale_adjust_x,scale_adjust_y;/* make large dots smaller */
   /* make small dots larger */
   int lx,rx,ty,by; /* positions for corners of dot */
   x_pos=gif_col_to_x(dis_l,x_start,x_pic_per_dot,column);
   y_pos=gif_row_to_y(dis_t,y_start,y_pic_per_dot,row);
   if((x_pic_per_dot<1.2)||(y_pic_per_dot<1.2))
     {/* draw a dot here */
      fprintf(giffp,"setpixel im %d %d %s \n",x_pos,y_pos,color);
     }
   else /* draw a rectangle here 94% of the proper width */
     {if(global_contains_prob==TRUE) /* adjust for value of prob */
       {if(energy<=-12000)
         energy=-12000;
        scale_adjust_x=sqrt((-1.0*(float)energy)/10000.);
        scale_adjust_y=scale_adjust_x;
       }
      else
       {if(x_pic_per_dot>3.)
         scale_adjust_x=.85;
        else
          scale_adjust_x=1;
        if(y_pic_per_dot>3.)
           scale_adjust_y=.80; /* make big dots smaller */
        else
          scale_adjust_y=1.1; /* make small dots bigger */
       } 
      ix_dot_size=(int)(x_pic_per_dot*scale_adjust_x+.7);
      if(ix_dot_size<1)
        ix_dot_size=1;
      iy_dot_size=(int)(y_pic_per_dot*scale_adjust_y+.7);
      if(iy_dot_size<1)
        iy_dot_size=1;
      lx=x_pos-ix_dot_size/2;
      rx=lx+ix_dot_size-1;
      ty=y_pos-iy_dot_size/2;
      by=ty+iy_dot_size-1;
      fprintf(giffp,"filledrectangle im %d %d %d %d %s\n",lx,ty,rx,by,color);
     } 
}


/* create the tic marks */
void gif_single_hor_tic(int x_start,int j_col,float x_pic_per_dot,int ty,
                        int by,int dis_l,int make_grid_lines)
{ int x_pos;
  x_pos=gif_col_to_x(dis_l,x_start,x_pic_per_dot,j_col);
  if(make_grid_lines==TRUE)
    fprintf(giffp,"line im %d %d %d %d gl\n", x_pos,ty,x_pos,by-2);
  fprintf(giffp,"line im %d %d %d %d b0\n",x_pos,ty,x_pos,ty-15);
  fprintf(giffp,"string im gdFont10x20 %d %d \"%d\" b0\n",
                         x_pos-10,ty-36,j_col);
}


void gif_plot_hor_tics(int x_start,float x_pic_per_dot,int dis_l,
                       int dis_r,int ty,int by,int make_grid_lines)
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
                            make_grid_lines);
      }
}

void gif_single_ver_tic(int y_start,int i_row,float y_pic_per_dot,int lx,
                   int rx,int dis_t,int make_grid_lines)
{ int y_pos; 
  y_pos=gif_row_to_y(dis_t,y_start,y_pic_per_dot,i_row);
  if(make_grid_lines==TRUE)
    fprintf(giffp,"line im %d %d %d %d gl\n",rx-2,y_pos,lx+2,y_pos );
  fprintf(giffp,"line im %d %d %d %d b0\n",rx,y_pos,rx+15,y_pos);
  fprintf(giffp,"string im gdFont10x20 %d %d \"%d\" b0\n",rx+18,
                          y_pos-12,i_row);
}


void gif_plot_ver_tics(int y_start,float y_pic_per_dot,int dis_t,
                       int dis_b,int lx,int rx,int make_grid_lines)
{ int i_row,step,start;
  step=step_fun(dis_t,dis_b);
  start=start_fun(dis_t,step);
  gif_single_ver_tic(y_start,dis_t,y_pic_per_dot,lx,rx,dis_t,FALSE);
  gif_single_ver_tic(y_start,dis_b,y_pic_per_dot,lx,rx,dis_t,FALSE);
  for(i_row=start;i_row<(dis_b-step/2);i_row=i_row+step)
   {gif_single_ver_tic(y_start,i_row,y_pic_per_dot,lx,rx,dis_t,make_grid_lines);  
   }
}




void  gif_tics(int x_start,int y_start,float x_pic_per_dot,float y_pic_per_dot,
               int dis_l,int dis_r,int dis_t,int dis_b,int lx,
               int rx,int ty,int by,int make_grid_lines)
{  gif_plot_ver_tics(y_start,y_pic_per_dot,dis_t,dis_b,lx,rx,make_grid_lines);
   gif_plot_hor_tics(x_start,x_pic_per_dot,dis_l,dis_r,ty,by,make_grid_lines);
}


  /* make the basepairs */
int  gif_points(int x_start,int y_start,float x_pic_per_dot,
                float y_pic_per_dot,int dis_l,int dis_r,
                  int dis_t,int dis_b)
{
 int row,column,new_row;
 int diag,position,location,new_color;
 int color;
 char color_string[4];
 int count;
 int energy;
 int start_diag,end_diag;
 int points_plotted;
 start_diag=dis_t+dis_l-1;
 end_diag=dis_b+dis_r-1;
 points_plotted=0; 
 for(color=global_number_of_colors;color>=1;color--)
   {
    strcpy(color_string,gif_color[color]);
    for(diag=start_diag;diag<end_diag;diag++)
       {position=global_diag_start[diag];
        for(count=0;count<global_diag_count[diag];count++)
          {location=position+count;
           if((!global_opt_prob_flag)
                      ||(global_diag[location].energy>=-20000))
		    {
          if(global_diag[location].row<=dis_b)
	  {
           if((global_diag[location].energy==global_optimal_energy)
              ||((global_diag[location].energy<=global_energy_cutoff)&&
             (global_diag[location].length>=global_chain_len)))
	     {new_color=global_diag[location].color;
              if(new_color==color)
	      { new_row=global_diag[location].row;
                column=global_diag[location].column;
                for(row=new_row;
                      row<(new_row+global_diag[location].length);row++)
	           {if((row>=dis_t)&&(row<=dis_b)
                          &&(column>=dis_l)&&(column<= dis_r))
                          { gif_plot_dot(x_start,y_start,x_pic_per_dot,
                                 y_pic_per_dot,dis_l,dis_t,color_string,
                            row,column,global_diag[location].energy);
                            points_plotted++;
                           }
                             column--;
		   }
              }                        

	     }
	  }}
	  }
       }
    }
 /* plot dots below main diagonal */
  color=1;
  strcpy(color_string,"OP");
  for(diag=start_diag;diag<end_diag;diag++)
       {position=global_diag_start[diag];
        for(count=0;count<global_diag_count[diag];count++)
          {
            location=position+count;
            if(global_diag[position+count].row<=dis_r)
	      {
            if (global_diag[location].color==1)
             {location=position+count;
	      new_row=global_diag[location].row;
                column=global_diag[location].column;
                if(column>=dis_t) 
		  {
                for(row=new_row;
                      row<(new_row+global_diag[location].length);row++)
	           {if((column>=dis_t)&&(column<=dis_b)
                          &&(row>=dis_l)&&(row<= dis_r))
                          {energy=global_diag[location].energy;
                           if(global_contains_prob)
                             {if(global_opt_prob_flag)
                               {if(energy<=-20001)
                                 gif_plot_dot(x_start,y_start,x_pic_per_dot,
                                   y_pic_per_dot,dis_l,dis_t,color_string,
                                  column,row,-10000);
			       }
                             else
                               gif_plot_dot(x_start,y_start,x_pic_per_dot,
                                   y_pic_per_dot,dis_l,dis_t,color_string,
                                  column,row,energy);
			     }
                          else
                            gif_plot_dot(x_start,y_start,x_pic_per_dot,
                                   y_pic_per_dot,dis_l,dis_t,color_string,
                                  column,row,energy);
			  }

                             column--;
		   }
                  }
        
	     }}
  
	  }
       }
 return points_plotted;
}







/* make diagonal line */
void gif_make_diag(int dis_l,int dis_r,int dis_t,int dis_b,
 int  x_start,int start_y,float x_pic_per_dot,float y_pic_per_dot,
 int resolution)
{int ty,rx;
 int lx,by;
 int right,left,top,bottom;
 if((!((dis_r>=dis_t)&&(dis_b>=dis_l)))) /* cases all below, all above */
   return; /* There is no diagonal line */
 /* find right end of line */
 if(dis_r>dis_b)
   { /* diagonal crosses bottom border */
     bottom=dis_b;
     right=bottom;
   }
 else
   { /* diagonal crosses right border */
     right=dis_r;
     bottom=right;
   }
 /* find left end of line */
 if(dis_l>dis_t)
   {/* diagonal crosses left side */
    left=dis_l;
    top=dis_l;
   }
 else
   {/* diagonal crosses top side */
    top=dis_t;
    left=top;
   }
 ty=gif_row_to_y(dis_t,start_y,y_pic_per_dot,top)-
              (int)(y_pic_per_dot/2.+.5);
 lx=gif_row_to_y(dis_l,x_start,x_pic_per_dot,left)-
              (int)(x_pic_per_dot/2.+.51);
 rx=gif_col_to_x(dis_l,x_start,x_pic_per_dot,right)+
               (int)(x_pic_per_dot/2.+.51);
 by=gif_row_to_y(dis_t,start_y,y_pic_per_dot,bottom)+
               (int)(y_pic_per_dot/2.+.5);
 if(resolution<70) /* draw a narrow or wide line */
       fprintf(giffp,"line im %d %d %d %d b0\n",lx-1,ty-1,rx,by);
 else
       fprintf(giffp,"line im %d %d %d %d gdBrushed \n",lx-1,ty-1,rx,by);
} 


void gif_create_colors(char *file_type,int grid_flag)
{int color[TOTAL_COLORS][3];
 int i,j;
  /* convert colors to integers */
  for(j=0;j<TOTAL_COLORS;j++)
    {for(i=0;i<=2;i++)
       color[j][i]=(int)(255.*global_color[j][i]);
    }
  if(global_number_of_colors>=7)
    {for(j=1;j<=global_number_of_colors;j++)
      {for(i=0;i<=2;i++)
        color[j][i]=(int)(255.*global_prob_color[j][i]);
      }
    }
  fprintf(giffp,"create bkim 2 2\n");
  fprintf(giffp,"colorallocate bkim bak %d %d %d\n",
    color[COLOR_BACKGROUND][0],
    color[COLOR_BACKGROUND][1],color[COLOR_BACKGROUND][2]);
  fprintf(giffp,"colorallocate bkim b0 %d %d %d\n",
    color[0][0],color[0][1],color[0][2]);
  fprintf(giffp,"line bkim 0 0 1 0 b0\n");
  fprintf(giffp,"line bkim 0 1 1 1 b0\n");
  /* end of brush definition */
  fprintf(giffp,"colorallocate im bak %d %d %d\n",
       color[COLOR_BACKGROUND][0],color[COLOR_BACKGROUND][1],
       color[COLOR_BACKGROUND][2]);
  /* first is background color */
  fprintf(giffp,"colorallocate im b0 %d %d %d\n",
                       color[0][0],color[0][1],color[0][2]);
  fprintf(giffp,"colorallocate im bk %d %d %d\n",
            color[1][0],color[1][1],color[1][2]);
  /* special optimal color */
 fprintf(giffp,"colorallocate im OP %d %d %d\n",
             (int)(255.*global_color[1][0]), 
             (int)(255.*global_color[1][0]), 
             (int)(255.*global_color[1][0]));
  /* declare colors */
  fprintf(giffp,"colorallocate im re %d %d %d\n",
                      color[2][0],color[2][1],color[2][2]);
  fprintf(giffp,"colorallocate im gr %d %d %d\n",
                      color[3][0],color[3][1],color[3][2]);
  fprintf(giffp,"colorallocate im ye %d %d %d\n",
                      color[4][0],color[4][1],color[4][2]);
  fprintf(giffp,"colorallocate im pu %d %d %d\n",
                   color[5][0],color[5][1],color[5][2]);
  fprintf(giffp,"colorallocate im br %d %d %d\n",
                   color[6][0],color[6][1],color[6][2]);
  fprintf(giffp,"colorallocate im bl %d %d %d\n",
                   color[7][0],color[7][1],color[7][2]);
  fprintf(giffp,"colorallocate im fl %d %d %d\n",
                   color[8][0],color[8][1],color[8][2]);
  fprintf(giffp,"colorallocate im la %d %d %d\n",
                  color[COLOR_LABEL][0],color[COLOR_LABEL][1],
                  color[COLOR_LABEL][2]);
   fprintf(giffp,"colorallocate im co %d %d %d\n",
                  color[COLOR_COPY][0],color[COLOR_COPY][1],
                  color[COLOR_COPY][2]);

  if(grid_flag==TRUE)
         fprintf(giffp,"colorallocate im gl %d %d %d\n",
                  color[COLOR_GRID][0],color[COLOR_GRID][1],
                  color[COLOR_GRID][2]);
               /* Gray for grid lines */
  fprintf(giffp,"setbrush im bkim\n");/* Sets the brush  for wide lines*/
  if(strcmp(file_type,"jpg"))
   {
    fprintf(giffp,"interlace im 1\n"); /* $$$$$ turn on interlace */
   }

}



void gif_label_points(int label_i,int label_j,int x_start,int y_start,
      float x_pic_per_dot,float y_pic_per_dot,int dis_l,int dis_t)
{int hor_pos,ver_pos,label_hor_pos,lx,ty,rx,by;
 int hor_shift,ver_shift;
 char string[100];
 /* printf("\n within label subroutine i is %d, j is %d",label_i,label_j); */
 hor_shift=(int)(x_pic_per_dot/3);
 ver_shift=(int)(y_pic_per_dot/3);
 hor_pos=gif_col_to_x(dis_l,x_start,x_pic_per_dot,
                  label_j);
 ver_pos=gif_row_to_y(dis_t,y_start,y_pic_per_dot,
                  label_i);
 label_hor_pos=gif_col_to_x(dis_l,x_start,x_pic_per_dot,
                  label_j+1);
 lx=hor_pos-hor_shift;
 rx=hor_pos+hor_shift;
 ty=ver_pos-ver_shift;
 by=ver_pos+ver_shift;
 strcpy(string,num_string_int(label_i));
 strcat(string,",");
 strcat(string,num_string_int(label_j)); 
 fprintf(giffp,"filledrectangle im %d %d %d %d la\n",lx,ty,rx,by);
 fprintf(giffp,"string im gdFontSmall %d %d \"(%s) \"la\n",
          label_hor_pos+1,ver_pos-5,string);
}





void display_title_message_gif(void)

{ char year[8];
  display_time_year(year);
  /* Try not to make text overwrite image below */
  
  fprintf(giffp,"string im gdFont6x9 %d %d \"%s\" %s\n",
                           6,1,copyright_message1,gif_color[COLOR_COPY]);
  fprintf(giffp,"string im gdFont6x9 %d %d \"  %s%s\" %s\n",
                           6,18,year,copyright_message2,gif_color[COLOR_COPY]);
  /* put in copyright symbol */
  fprintf(giffp,"string im gdFont12x24 5 9 \"O\" %s\n",gif_color[COLOR_COPY]);
  fprintf(giffp,"string im gdFont6x13 8 15 \"C\" %s\n",gif_color[COLOR_COPY]);
}



int open_temp(char *filename)
{    /* open specified output file */
       /*  filename[strlen(filename) - 1] = NULL; */
   if ((giffp = fopen(filename, "w")) == NULL) /* open a file */
         {printf ("\n * Could not open file:  %s", filename);
          return(1);
          }
    return (0);
}   


int open_gif(char *filename)
{    /* open specified output file */
       /*  filename[strlen(filename) - 1] = NULL; */
   if ((giffp = fopen(filename, "w")) == NULL) /* open a file */
         {printf ("\n * Could not open file:  %s", filename);
          return(1);
          }
    return (0);
}   


void make_gifdat(char *filename,int resolution,int x_start,int y_start,
      float x_pic_per_dot,float y_pic_per_dot,int dis_l,int dis_t)
{int error;
 int len;
 len=strlen(filename);
/* strip off .tgd from name and replace with .gifdat */
 if(len>4)
   {filename[len-4]='\0';
   }
 strcat(filename,".gifdat");
 printf("\n GIFDAT file is %s",filename);
 error=open_gif(filename);
 if(error)
   {printf("Error creating .gifdat file name was %s",filename);
    return;
   }
 fprintf(giffp,"%d %d %d %f %f %d %d %d \n",resolution,x_start,y_start,
                x_pic_per_dot,y_pic_per_dot,global_length,dis_l,dis_t);
 fclose(giffp); 
}



void finish_making_gif(int save_tgd_flag,int clear_flag,int resolution,
                       char *name_for_gif,char *title,
                       int make_gifdat_flag,
                     int make_label_flag,int label_i,int label_j,
                     int make_grid_flag,int png_mode,int jpg_mode,
                     int mi_flag,int opt_prob_flag,int label_array_count,
                     int *label_row_array,int *label_col_array)
{ char action[90]; /* needed for system call */
  int error; /* resoulution is 72, 110, or 200 ? */
  int error_occurred;
  int hor_pos;
  int name_len;
  int color;
  int i;
  char col_str[3];
  char file_type[6];
  char temp_filename[120];
  char gifeng_filename[120];
  int offset,offset_adjust;
  int dis_l,dis_r,dis_t,dis_b;/* local versions of display_l,display_r */
  char string[120],gif_filename[120];
  int x_float; 
  /* use 8.5 wide by 11 inches wide gif ? */
  int gif_width; /* width of image in pixels */
  int gif_height; /* height of image in pixels */
  int lx,ty,rx,by; /* top,bottom,right,left corn of rec */
  float x_pic_per_dot,y_pic_per_dot;
  int x_start,y_start;
  int points_plotted;
  dis_l=display_l;
  dis_r=display_r;
  dis_t=display_t;
  dis_b=display_b;
  /* correct message for text window */
  strcpy(gif_filename,name_for_gif);
  if(strlen(gif_filename)==0)
      {printf("\n Error, filename has zero length.");
       return; /* abort making file */
      }
  gif_width=(int)(8.5*resolution)-1;
  gif_height=11*resolution-1;
  strcpy(temp_filename,gif_filename);
  name_len=strlen(temp_filename);
  if(name_len>4)
     temp_filename[name_len-4]='\0';
  if(make_label_flag==TRUE)
    strcpy(gifeng_filename,temp_filename);
  strcat(temp_filename,".tgd");
 if(png_mode)
      {strcpy(file_type,"png");
      }
     else
       {if(jpg_mode)
         strcpy(file_type,"jpg");
        else
         strcpy(file_type,"gif");
       }
  error=open_gif(temp_filename);
  if(error)
    {printf("\n error for %s file ",file_type);
     }
  /* hor by ver or x by y */
  fprintf(giffp,"create im %d %d\n",gif_width+1,gif_height+1);
                     /* declare image size in pixels */ 
  printf("\n Creating %s %d by %d \n",file_type,gif_width+1,gif_height+1);
  gif_create_colors(file_type,make_grid_flag);
  if(clear_flag==TRUE)
    fprintf(giffp,"colortransparent im bak\n");
  /* make box  7 by 7 inches */
  ty=(int)(1.97*resolution); /* come down 3 inches */
  by=(int)(gif_height-1.97*resolution);
  lx=(int)(.47*resolution); /* left side of box */ 
  rx=gif_width-(.95*resolution); /* 6.5 inches wide drawing region*/
  /* an extra .1 inch was used to make the line outside the region */
  if(resolution<70)
    fprintf(giffp,"polygon im %d %d %d %d %d %d %d %d 4 b0\n",
                   lx,ty,rx,ty,rx,by,lx,by);
  else  /* Make a wide line */
    {fprintf(giffp,"line im %d %d %d %d gdBrushed\n",lx,ty,rx,ty);
     fprintf(giffp,"line im %d %d %d %d gdBrushed\n",lx,by,rx,by);
     fprintf(giffp,"line im %d %d %d %d gdBrushed\n",lx,ty,lx,by);
     fprintf(giffp,"line im %d %d %d %d gdBrushed\n",rx,ty,rx,by);
    }
  hor_pos=(int)(4.25*resolution-(strlen(title)+2)/2*12);
  if(hor_pos<1)
     hor_pos=1;
  clean_title_gif(title);
  fprintf(giffp,"string im gdFont12x24 %d %d \"%s\"  b0\n",
              hor_pos,(int)(resolution/2),title);

  x_pic_per_dot=resolution*7.0/(dis_r-dis_l+1);
   /* width of basepair*/
  y_pic_per_dot=resolution*7.0/(dis_b-dis_t+1);/* height of basepair*/
  x_start=(int)(.5*resolution+x_pic_per_dot/2+.5);
   /* distance from left edge for first basepair*/
  y_start=(int)(2.0*resolution+y_pic_per_dot/2+.5);
     /* distance from top for first basepair*/
  /* both above are in pixels */
  /* create the tic marks */
  /*                      */
  gif_tics(x_start,y_start,x_pic_per_dot,y_pic_per_dot,dis_l,
                    dis_r,dis_t,dis_b,lx,rx,ty,by,make_grid_flag);
 /* make the diagonal line */
  gif_make_diag(dis_l,dis_r,dis_t,dis_b,x_start,
              y_start,x_pic_per_dot,y_pic_per_dot,resolution);
  /* both above are in pixels */
  /* create the points */
  points_plotted=gif_points(x_start,y_start,x_pic_per_dot,y_pic_per_dot,
         dis_l,dis_r,dis_t,dis_b);
  fprintf(giffp,
         "string im gdFont10x20 %d %d \"Lower Triangle Shows\" OP\n",
                     lx,by+15);
 if(global_contains_prob)
    {if(opt_prob_flag)
      {fprintf(giffp,
         "string im gdFont10x20 %d %d \"Optimal Structure\" OP\n",
                     lx,by+35);
      }
     else
       {fprintf(giffp,
         "string im gdFont10x20 %d %d \"Highest Color Range\" OP\n",
                     lx,by+35);
          
       }
    }
  else
    {fprintf(giffp,
       "string im gdFont10x20 %d %d \"Optimal Energy\" OP\n",
                     lx,by+35); 
    }

  fprintf(giffp,
    "string im gdFont10x20 %d %d \"Upper Triangle\" b0\n",
              lx,by+68);
  fprintf(giffp,
    "string im gdFont10x20 %d %d \"Basepairs Plotted: %d\" b0\n",
              lx,by+91,points_plotted);
  if(global_contains_prob!=TRUE)
    {strcpy(string,
    num_string_fancy_float(global_energy_cutoff-global_optimal_energy));
     strcat(string," kcal/mole");
    if(global_chain_len>1)
      {strcat(string,"     Filter: ");
       strcat(string,num_string_int(global_chain_len));
      }
    fprintf(giffp,
       "string im gdFont10x20 %d %d \"deltaG in Plot File =%s \" b0\n",
       resolution/2+5,resolution,string);
    }
  else
    {strcpy(string,num_string_fancy_float(global_energy_cutoff));
     if(mi_flag)
      strcat(string," <= M.I. <= ");
     else
      strcat(string," <= Probability <= ");
     strcat(string,
     num_string_fancy_float(global_optimal_energy));
     fprintf(giffp,"string im gdFont10x20 %d %d \" %s \" b0\n",
           resolution/2+5,resolution,string);
    }
  /* display optimal energy */
  if(global_contains_prob!=TRUE)
     {strcpy(string,num_string_fancy_int(global_optimal_energy));
      fprintf(giffp,
         "string im gdFont10x20 %d %d \"Optimal energy: %s \"bk\n",
          (int)(4.2*resolution)+14,by+3,string);
     }
  if(global_contains_prob!=TRUE)
    strcpy(string,
    num_string_fancy_float(global_energy_cutoff-global_optimal_energy));
  else
    strcpy(string,
    num_string_fancy_float(global_energy_cutoff));
  if(global_contains_prob!=TRUE)
     strcat(string," kcal/mole");
  if(global_chain_len>1)
    {strcat(string,"         Filter: ");
     strcat(string,num_string_int(global_chain_len));
    }
  /* print the color scheme */
 if(global_contains_prob!=TRUE)
    {x_float=(float)global_optimal_energy;
     if((global_number_of_colors<6)||(resolution>90))
       offset_adjust=23;
     else
       offset_adjust=17;
     for(color=2,offset=22;color<=global_number_of_colors;
       color++,offset+=offset_adjust)
       {strcpy(col_str,gif_color[color]);
        strcpy(string,num_string_fancy_float(x_float+(float)(color-2)*
                        global_color_increment));
        strcat(string,"< energy <=");
        strcat(string,num_string_fancy_float(x_float+(float)(color-1)*
                     global_color_increment));
        fprintf(giffp,"string im gdFont10x20 %d %d \"%s \" %s\n",
          (int)(4.2*resolution),by+offset,
          string,col_str); 
       }
    }
 else
   {
    strcpy(string,num_string_float(global_prob[global_number_of_colors-4][0]));
    strcat(string," < Prob.");
    fprintf(giffp,"string im gdFont10x20 %d %d \"%s \" bk\n",
          (int)(4.2*resolution),by+6,
          string);
    strcpy(col_str,"re");
    if((global_number_of_colors<6)||(resolution>90))
       offset_adjust=23;
     else
       offset_adjust=17;
    for(color=2,offset=6+offset_adjust;color<=global_number_of_colors-1;
                color++,offset+=offset_adjust)
	{strcpy(string,
            num_string_float(global_prob[global_number_of_colors-4][color-1]));
	 strcat(string," < Prob. <=");
	 strcat(string,
           num_string_float(global_prob[global_number_of_colors-4][color-2])); 
         fprintf(giffp,"string im gdFont10x20 %d %d \"%s \" %s\n",
            (int)(4.2*resolution),by+offset,
            string,col_str);
          strcpy(col_str,gif_color[color+1]); 
	  }
    strcpy(string,"         Prob. <=");
    strcat(string,num_string_float(
     global_prob[global_number_of_colors-4][global_number_of_colors-2]));
    
    fprintf(giffp,"string im gdFont10x20 %d %d \"%s \" %s\n",
          (int)(4.2*resolution),by+offset,
          string,col_str);
 
   }
  if(make_label_flag==TRUE)
    {gif_label_points(label_i,label_j,x_start,y_start,
      x_pic_per_dot,y_pic_per_dot,dis_l,dis_t);
     }
  if(label_array_count>0)
    {for(i=0;i<label_array_count;i++)
       gif_label_points(label_row_array[i],label_col_array[i],x_start,y_start,
      x_pic_per_dot,y_pic_per_dot,dis_l,dis_t);
    }
  display_title_message_gif();
  strcpy(string,"gif im ");
  strcat(string,gif_filename); 
  strcat(string,"\n");
  fprintf(giffp,string);
  fclose(giffp);
  printf(" %s file is            %s\n",file_type,gif_filename);
 /* convert tgd file to gif */
  if(png_mode)
    strcpy(action,"tgd_png ");
  else
    {if(jpg_mode)
      strcpy(action,"tgd_jpg ");
     else
       strcpy(action,"tgd ");
    }
  strcat(action,temp_filename);
  error_occurred=system(action);
  if(error_occurred!=0)
    {printf("\n Error converting the file to %s format",file_type);
     printf("\n The program executes tgd_png or tgd_jpg to convert");
     printf(" the file to %s format.",file_type);
     printf("\n Either it was not installed or it was not installed ");
     printf("\n correctly. ");
     printf("\n See boxplot_ng.doc or boxplot.doc for more information");
     printf("\n ");
     printf("\n %s should be deleted",temp_filename);
    }
  else 
    {if(save_tgd_flag!=TRUE)
        {strcpy(action,"rm ");
         strcat(action,temp_filename);
         system(action);
        }
     else
	{ printf("\n You may edit %s",temp_filename);
          printf("\n Run tgd on %s",temp_filename);
          printf("\n Please delete %s when finished\n",temp_filename);
         }
    } 
   if(make_gifdat_flag==TRUE)
     {make_gifdat(temp_filename,resolution,x_start,y_start,x_pic_per_dot,
         y_pic_per_dot,dis_l,dis_t);
     }
   if(make_label_flag==TRUE)
     make_energy_file(gifeng_filename,label_i,label_j);
}  
