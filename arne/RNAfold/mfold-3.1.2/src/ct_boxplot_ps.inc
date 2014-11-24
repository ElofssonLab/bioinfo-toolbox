/* Darrin Stewart and Michael Zuker */
/* May 3, 1999  Version 2.0*/
/* ct_boxplot_ps.c */
/* create postscript output for ct_boxplot */

#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "copyright_ps.inc"


FILE *postfp;

/* Set the color */
void fixcolor_plot_postscript(int color)
{ /*  postscript colors */
  if((color>=0)&&(color<COLOR_TABLE_SIZE))
      fprintf(postfp,"\n  color_%d\n",color);
}

void fixcolor_plot_postscript_define(int last_structure,int lower_structure,
                   float color_table[][3],int nongray_for_3_flag) 
{ /*  postscript colors */
  int color;
  fprintf(postfp,"\n%s Set colors","%%");
  fprintf(postfp,"\n%s 0 is full overlap, 1 is partial overlap","%%");
  fprintf(postfp,"\n%s 2 is structure 1, 3 is structure 2, 4 is structure 3,etc ","%%");
  fprintf(postfp,"\n%s %d is text, ","%%",MAXIMUM_CT_FILES+2);
  fprintf(postfp,"\n%s %d is background, %d is labels, %d is grid lines","%%",
               MAXIMUM_CT_FILES+3, MAXIMUM_CT_FILES+4, MAXIMUM_CT_FILES+5);  
  if((nongray_for_3_flag)||(lower_structure<=-10))
     fprintf(postfp,"\n%s 5,6,7 are for interpreted partial overlap","%%");
  fprintf(postfp,"\n%s original colors from ct_boxplot.col\n\n","%%");
  for(color=0;color<=last_structure;color++)
    {fprintf(postfp,"\n/color_%d {  %.3f %.3f %.3f setrgbcolor} def ",color,
        color_table[color][0],color_table[color][1],color_table[color][2]);
    } 
  if((nongray_for_3_flag)||(lower_structure<=-10))
    for(color=(last_structure+1);color<=(last_structure+3);color++)
    {fprintf(postfp,"\n/color_%d {  %.3f %.3f %.3f setrgbcolor} def ",color,
        color_table[color][0],color_table[color][1],color_table[color][2]);
    }
  for(color=MAXIMUM_CT_FILES+2;color<COLOR_TABLE_SIZE;color++)
  {fprintf(postfp,"\n/color_%d {  %.3f %.3f %.3f setrgbcolor} def ",color,
        color_table[color][0],color_table[color][1],color_table[color][2]);
    } 
  fprintf(postfp,"\n");
}



/* open the output file */
int open_post(char *post_filename)
{    /* open specified output file */
   if ((postfp = fopen(post_filename, "w")) == NULL) /* open a file */
         {printf ("\n * Could not open file:  %s", post_filename);
          return(1);
          }
    return (0);
}

void single_hor_tic_zoom_post(int j_col,float scale,int grid_flag,
             int display_l)
{ fprintf(postfp,"\n %f 1. scale",scale);
  fprintf(postfp,"\n %d 540 moveto",j_col-display_l+1);
  fprintf(postfp,"\n 0 8 rlineto");
  if(grid_flag==TRUE)
   {fprintf(postfp,"\n stroke");
    fixcolor_plot_postscript(COL_GRID);
    fprintf(postfp,"\n %d 539 moveto",j_col-display_l+1);
    fprintf(postfp,"\n 0 -538 rlineto");
    fprintf(postfp,"\n stroke");
    fixcolor_plot_postscript(COL_TEXT);
   }
  fprintf(postfp,"\n %d 550 moveto",j_col-display_l+1);
  fprintf(postfp,"\n %f 1. scale",1./scale);
  fprintf(postfp,"\n -4 0 rmoveto");
  fprintf(postfp,"\n (%d) show",j_col);
}
void plot_hor_tics_zoom_post(float scale,int grid_flag,int display_l,
                    int display_r)
{ 
  int j_col,step;
  int start;
  step=step_fun(display_l,display_r);
  start=start_fun(display_l,step);
  single_hor_tic_zoom_post(display_r,scale,FALSE,display_l);
  single_hor_tic_zoom_post(display_l,scale,FALSE,display_l);
  for(j_col=start;j_col<(display_r-2*step/3);
                        j_col=j_col+step)
      {single_hor_tic_zoom_post(j_col,scale,grid_flag,display_l);
      }
}
void single_ver_tic_zoom_post(int i_row,float scale,int grid_flag,
               int display_b)
{ fprintf(postfp,"\n 1. %f scale",scale);
  fprintf(postfp,"\n 540.5 %d moveto",display_b-i_row+1);
  fprintf(postfp,"\n 8 0 rlineto");
 if(grid_flag==TRUE)
   {fprintf(postfp,"\n stroke");
    fixcolor_plot_postscript(COL_GRID);
    fprintf(postfp,"\n 539.5 %d moveto",display_b-i_row+1);
    fprintf(postfp,"\n -538 0 rlineto");
    fprintf(postfp,"\n stroke");
    fixcolor_plot_postscript(COL_TEXT);
   }
  fprintf(postfp,"\n 550.5 %d moveto",display_b-i_row+1);
  fprintf(postfp,"\n 1. %f scale",1./scale);
  fprintf(postfp,"\n (%d) show",i_row);
}
void plot_ver_tics_zoom_post(int display_t,int display_b,float scale,
                      int grid_flag)
{ int i_row,step,start;
 step=step_fun(display_t,display_b);
 start=start_fun(display_t,step);
 single_ver_tic_zoom_post(display_t,scale,FALSE,display_b);
 single_ver_tic_zoom_post(display_b,scale,FALSE,display_b);
 for(i_row=start;i_row<(display_b-step/2);i_row=i_row+step)
   {single_ver_tic_zoom_post(i_row,scale,grid_flag,display_b);  
   }
}

void draw_diagonal_line_post(int display_l,
         int display_r,int display_t,int display_b)
{int left,top,bottom,right;
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
     fixcolor_plot_postscript(COL_TEXT);
     fprintf(postfp,"\n %.4f  %.4f moveto", 
             left-display_l+0.0, display_b-top+2.0);
     fprintf(postfp,"\n %.4f %.4f  lineto stroke",
        right-display_l+2.0,display_b-bottom+0.0);
    }
}

void display_time(void)
{time_t now;
 now=time(NULL);
 fprintf(postfp,"\n%s %s","%%",ctime(&now));
}


void draw_postscript_dots(int counter[MAXIMUM_CT_FILES+1],
               int stripe_flag,
               int last_structure,
               int *diag_start,int *diag_count,struct helix *diag,
               int display_l,int display_r,int display_t,int display_b,
               int lower_structure,int nongray_for_3_flag,
               int *gray3_content,int gray3[3])
{int row,col,new_row,count;
 int diagonal,position,location,new_color;
 int color;
 int starting_color,ending_color;
 int start_diag,end_diag;
 float adjuster;
 int n_col,column;
 int old_color;
  gray3[0]=0;
  gray3[1]=0;
  gray3[2]=0;
 start_diag=display_t+display_l-1;
 end_diag=display_b+display_r-1;
 /* single structure colors above main diagonal */
 for(color=last_structure;color>=2;color--)
 {if(stripe_flag==TRUE)
    {adjuster=(color-1.)/(last_structure-1.)*.5;
     }
  fixcolor_plot_postscript(color);
  for(diagonal=start_diag;diagonal<=end_diag;diagonal++)
    {position=diag_start[diagonal];
     
     for(count=0;( (count<diag_count[diagonal]) 
   &&(diag[position+count].row<=display_b) ) ;count++)
       { location=position+count;
	  new_color=diag[location].color;
             if(new_color==color)
	       {new_row=diag[location].row;
                col=diag[location].column;
                for(row=new_row;
                    row<(new_row+diag[location].length);row++)
		      {if((row>=display_t)&&(row<=display_b)
                          &&(col>=display_l)&&(col<= display_r))
                        {if(stripe_flag)
                          fprintf(postfp,"\n%d %d %.5f qb%d",
                               row-display_t+1,col-display_l+1,
                                 adjuster,color);
                         else
                          fprintf(postfp,"\n%d %d b%d",
                               row-display_t+1,col-display_l+1,color);
                        counter[color]++;  
                        }
                          col--;
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
   {fixcolor_plot_postscript(color); 
     if((color==1)&&(nongray_for_3_flag))
        fixcolor_plot_postscript(5); /* for overlap of 1 and 2 */
    for(diagonal=start_diag;diagonal<=end_diag;diagonal++)
     {position=diag_start[diagonal];
      for(count=0;( (count<diag_count[diagonal]) 
        &&(diag[position+count].row<=display_b) ) ;count++)
       {location=position+count;
	  new_color=diag[location].color;
             if(new_color==color)
	       {new_row=diag[location].row;
                col=diag[location].column;
                for(row=new_row;
                    row<(new_row+diag[location].length);row++)
		      {if((row>=display_t)&&(row<=display_b)
                          &&(col>=display_l)&&(col<= display_r))
                        {fprintf(postfp,"\n%d %d b%d",
                               row-display_t+1,col-display_l+1,color);
                         counter[color]++;
                        }
                          col--;
                      }                       

	       }
          
  
        }
    }
   }
 /* plot gray interpreted above main diagonal*/
 if(nongray_for_3_flag)
   {
     for(row=display_t;row<=display_b;row++)
     {
      col=gray3_content[row];
      if(col!=0)
       {if(col<0)
          n_col=-1*col;
        else
          n_col=col;
        if((n_col>=display_l)&&(n_col<=display_r)&&(row<=n_col))
	   { 
            if(col>0)
             {fixcolor_plot_postscript(7); /* for seq 2 nd 3 */
              gray3[2]++;
             }
            else
             {fixcolor_plot_postscript(6); /* for seq 1 & 3 */
              gray3[1]++;
             }
            fprintf(postfp,"\n%d %d b1",
                               row-display_t+1,n_col-display_l+1);
           
           }
	       
       }
     }
   gray3[0]=counter[1]-gray3[1]-gray3[2];
   }
 /* points below the main diagonal */
 /* draw dots on lower structure */
 if(lower_structure!=-10)
  {
 if(lower_structure!=-1) /* -1 is no dots below */
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
   {fixcolor_plot_postscript(color);
    for(diagonal=start_diag;diagonal<=end_diag;diagonal++)
     {position=diag_start[diagonal];
      for(count=0;( (count<diag_count[diagonal]) 
        &&(diag[position+count].row<=display_r) ) ;count++)
       {location=position+count;
	  new_color=diag[location].color;
             if(new_color==color)
	       {new_row=diag[location].row;
                col=diag[location].column;
                if(col>=display_t)
                 {
                for(row=new_row;
                    row<(new_row+diag[location].length);row++)
		      {if((col>=display_t)&&(col<=display_b)
                          &&(row>=display_l)&&(row<= display_r))
                        {fprintf(postfp,"\n%d %d b%d",
                               col-display_t+1,row-display_l+1,color);
                        }
                          col--;
                      } 
                    }                      

	       }
          
  
        }
    }
   }
  }
 }
if(lower_structure<=-10)
 {
  for(diagonal=start_diag;diagonal<=end_diag;diagonal++)
     {
      position=diag_start[diagonal];
      for(count=0;( (count<diag_count[diagonal]) 
        &&(diag[position+count].row<=display_r) ) ;count++)
       {location=position+count;
	  new_color=diag[location].color;
             if(new_color==1)
	       {new_row=diag[location].row;
                col=diag[location].column;
                if(col>=display_t)
                 {
                for(row=new_row;
                    row<(new_row+diag[location].length);row++)
		      {if((col>=display_t)&&(col<=display_b)
                          &&(row>=display_l)&&(row<= display_r))
                        { column=gray3_content[row];
                          if(column==0)
                            old_color=5;
                          else
                           {if(column<0)
                             { old_color=6;
                             }
                             else
                               {old_color=7;
                               }
                          }
                         fixcolor_plot_postscript(old_color);
                         if(!nongray_for_3_flag)
                             gray3[old_color-5]++;
                         fprintf(postfp,"\n%d %d b1",
                               col-display_t+1,row-display_l+1);
                        }
                          col--;
                      }
                    }                      

	       }
          
  
        }
    }
   }

}

void post_define_dots(float hscale,float vscale,int big_dots_flag,
       float dot_magnifier,int post_file,int zoom_labels_count,
       int stripe_flag,int last_structure)
{float new_scale;
 int i;
  if(big_dots_flag)
    {new_scale=(hscale+vscale)/2.0;
     /* 1.0 below is size in points of smallest dot */
     if((dot_magnifier*new_scale)<1.0)
         dot_magnifier=dot_magnifier/new_scale;
    }
 fprintf(postfp,"\n%s   Define dots \n","%%"); 
/* define all dot sizes s */
  fprintf(postfp,
       "\n\n%s Dot shape and Size definitions\n","%%"); 
  fprintf(postfp,
       "\n%s Adjust sl, si,zsi to change size of dots ","%%");
  fprintf(postfp,"\n%s 1.0 is maximum size without overlapping dots \n","%%");
  fprintf(postfp,"%s s0 = 2.0 doubles size of optimal dots, b1\n","%%");
  fprintf(postfp,"%s s1 = 2.0 doubles size of near optimal dots dots, b2\n",
      "%%");
  fprintf(postfp,"%s Use %.4f to make dots 1/72 inch or 1 point","%%",
       2.0/(hscale+vscale));
  fprintf(postfp,
     "\n%s s2 is 1st structure or file, s3 for 2nd structure","%%");
  if(stripe_flag) 
    {fprintf(postfp,
      "\n%s When partial overlap is shown in multiple colors on the ","%%");
     fprintf(postfp,
      "\n%s same dot, adjust zs2 for structure 1, zs3 for structure 2","%%");
     fprintf(postfp,
      "\n%s etc rather than s2,s3, etc ","%%");
 
    }
  fprintf(postfp,"\n%s Ignore Si,NSi,len  nzsi \n","%%");
  fprintf(postfp,"\n\n%s Adjust dot sizes here\n","%%");
  if((post_file!=MAIN_POST)&&(zoom_labels_count>0))
   {fprintf(postfp,"\n%s sl is for labels ,Normal is .5 \n","%%");
    fprintf(postfp,"/sl { .5 } def\n");
   } 
  for(i=0;i<=last_structure;i++)
       {
        fprintf(postfp,"/s%d { %.4f } def\n",i,dot_magnifier);
       }
  fprintf(postfp,"\n");
  if(stripe_flag)
   {for(i=2;i<=last_structure;i++)
      { fprintf(postfp,"/zs%d { %.4f } def\n",i,dot_magnifier);
      }
   }
 fprintf(postfp,"\n\n");
 if((post_file!=MAIN_POST)&&(zoom_labels_count>0)) 
  {fprintf(postfp,"/Sl {.50 sl mul } def\n");
    fprintf(postfp,"/NSl {-1.0 Sl mul } def\n");
   } /* define details of box */

 for(i=0;i<=last_structure;i++) 
       {fprintf(postfp,"/S%d {.50 s%d mul } def\n",i,i);
        fprintf(postfp,"/NS%d {-1.0 S%d mul } def\n",i,i);
       }
 if(stripe_flag)
  { 
   for(i=1;i<=last_structure;i++) 
       {
        fprintf(postfp,"/nzs%d {-1.0 zs%d mul } def\n",i,i);
       }
  }
/*
for(i=0;i<=last_structure;i++)
     {fprintf(postfp,"\n/b%d {",i); 
      fprintf(postfp,"\n  exch len exch sub 1 add ");
      fprintf(postfp,"\n moveto");
      fprintf(postfp,"\n   S%d S%d  rlineto",i,i);
      fprintf(postfp,"\n   NS%d  0 rlineto",i);
      fprintf(postfp,"\n   NS%d NS%d rlineto",i,i);
      fprintf(postfp,"\n   0   NS%d  rlineto",i);
      fprintf(postfp,"\n   S%d    0 rlineto",i);
      fprintf(postfp,"\n   S%d S%d rlineto",i,i);
      fprintf(postfp,"\n   0    S%d  rlineto",i);
      fprintf(postfp,"\n   closepath");
      fprintf(postfp,"\n   fill");
      fprintf(postfp,"\n   stroke");
      fprintf(postfp,"\n}def");
     }
*/
for(i=0;i<=last_structure;i++)
     {fprintf(postfp,"\n/b%d {",i); /*   */
      fprintf(postfp,"\n  exch len exch sub 1 add ");
      fprintf(postfp,"\n moveto");
      fprintf(postfp,"\n   S%d S%d  rmoveto",i,i);
      fprintf(postfp,"\n   NS%d  0 rlineto",i);
      fprintf(postfp,"\n   NS%d  0 rlineto",i);
      fprintf(postfp,"\n   0   NS%d  rlineto",i);
      fprintf(postfp,"\n   0   NS%d  rlineto",i);
      fprintf(postfp,"\n   S%d    0 rlineto",i);
      fprintf(postfp,"\n   S%d    0 rlineto",i);
      fprintf(postfp,"\n   0    S%d  rlineto",i);
      fprintf(postfp,"\n   0    S%d  rlineto",i);
      fprintf(postfp,"\n   closepath");
      fprintf(postfp,"\n   fill");
      fprintf(postfp,"\n   stroke");
      fprintf(postfp,"\n}def");
     }
 if(stripe_flag)
   {for(i=2;i<=last_structure;i++)
     { fprintf(postfp,"\n /qb%d {",i); /* for dots that vary in size */
        fprintf(postfp,"\n /s exch def");
         fprintf(postfp,"\n  exch len exch sub 1 add ");
        fprintf(postfp,"\n newpath");
         fprintf(postfp,"\n moveto");
         fprintf(postfp,"\n   .5 zs%d mul .5 zs%d mul rlineto",
               i,i); /* top right, OK */
          fprintf(postfp,"\n   -.5 s sub zs%d mul -.5 s add zs%d mul rlineto",
                i,i); /* top left */
         fprintf(postfp,"\n   s .5 sub zs%d mul s .5 add nzs%d mul rlineto",
                        i,i);
         fprintf(postfp,"\n   .5 s add zs%d mul  .5 s sub zs%d mul rlineto",
               i,i);
         fprintf(postfp,"\n .5 s sub zs%d mul .5 s add zs%d mul  rlineto",
                  i,i); /* top left again */

         fprintf(postfp,"\n   closepath");
         fprintf(postfp,"\n   fill");
         fprintf(postfp,"\n   stroke");
         fprintf(postfp,"\n}def");
     }
   }
 fprintf(postfp,"\n%s sbox is for labels\n","%%"); 
 if((post_file!=MAIN_POST)&&(zoom_labels_count>0))
    {fprintf(postfp,"\n/sbox {"); /* for labels */
     fprintf(postfp,"\n exch len exch sub 1 add ");
     fprintf(postfp,"\n moveto");
     fprintf(postfp,"\n  Sl Sl rlineto");
     fprintf(postfp,"\n   NSl   0  rlineto");
     fprintf(postfp,"\n   NSl    NSl rlineto");
     fprintf(postfp,"\n   0   NSl  rlineto");
     fprintf(postfp,"\n   Sl    0 rlineto");
     fprintf(postfp,"\n   Sl  Sl  rlineto");
     fprintf(postfp,"\n   0    Sl  rlineto");
     fprintf(postfp,"\n   closepath");
     fprintf(postfp,"\n   fill");
     fprintf(postfp,"\n   stroke");
     fprintf(postfp,"\n}def");
    }
}


void general_post(int big_dots_flag,int post_file,int zoom_labels_count,
                 int display_l,int display_r,
                 int display_t,int display_b,float color_table[][3],
                 int *zoom_labels_row,int*zoom_labels_col,
                     char *zoom_labels_string[],int grid_flag,
                 int last_structure,int *start,int *count,
                 struct helix* diag,int len,
                 int lower_structure,char *sequence_name,
                 char file_data_ct[][80],float dot_magnifier,
                 int forced_stripe_flag,int nongray_for_3_flag,
                 int *gray3_content)
{ float hscale,vscale;
  int i,color,stop_point,offset;
  float center;
  int stripe_flag;
  int gray3[3];
  int counter[MAXIMUM_CT_FILES+1];
  int other_stop_point;
  float new_scale;
  for(i=0;i<=last_structure;i++)
     counter[i]=0;
  if(post_file==MAIN_POST) 
     {hscale=540.0/(len+2); /* 540 is 72*7.5 for 7.5 inches of paper */
      vscale=hscale;
    }
  else
    {vscale=540.0/(display_b-display_t+2);
     hscale=540.0/(display_r-display_l+2);
    }
  fprintf(postfp,"%c!PS-Adobe-3.0 EPSF-2.0\n",'%');
  fprintf(postfp,"%sBoundingBox: 0 0 612 792\n","%%");
  fprintf(postfp,"\n%s Created by ct_boxplot Version 2.0 Darrin Stewart","%%");
  fprintf(postfp,"\n%s   and Michael Zuker","%%");
  display_time();
  fprintf(postfp,"%s\n%s\n","%%","%%");
  fprintf(postfp," %f %f scale\n",1.0,1.0);
  fprintf(postfp,"%s\n","%%");
  fprintf(postfp,"%s To change image size ....\n","%%");
  fprintf(postfp,"%s For full size image Use Bound ngBox:0 0 612 792\n","%%");
  fprintf(postfp,"%s                     and 1.0 1.0 scale\n","%%");
  fprintf(postfp,"%s For image .25 of full size, replace 612 with .25*612\n","%%");
  fprintf(postfp,"%s                replace 792 with .25*792 and replace\n","%%");
  fprintf(postfp,"%s                1.0 with .25: Other scales are similar.\n","%%");
  fprintf(postfp,"%s                Valid scales are .05 to 1.00\n","%%");
  fprintf(postfp,"%s ONLY scale and Bound ngBox lines above should be changed.\n",
       "%%");
  fprintf(postfp,"%s\n","%%");
  fixcolor_plot_postscript_define(last_structure,lower_structure,color_table,
             nongray_for_3_flag); 
  /* start of change */
  if(big_dots_flag)
    {new_scale=(hscale+vscale)/2.0;
     /* 1.0 below is size in points of smallest dot */
     if((dot_magnifier*new_scale)<1.0)
         dot_magnifier=dot_magnifier/new_scale;
    }
  if(forced_stripe_flag)
     stripe_flag=TRUE;
  else
   {
     if((post_file==MAIN_POST)||(last_structure<4))
       {stripe_flag=FALSE;
       }
     else
      {if(!nongray_for_3_flag)
        { if(  ((540.0/(display_b-display_t)/(last_structure-1)<STRIPE_PS) ) ||
          ((540.0/(display_r-display_l)/(last_structure-1.)<STRIPE_PS) )  )
            {stripe_flag=FALSE;
             }
          else
           {stripe_flag=TRUE;
           }
        }
       else
        stripe_flag=FALSE;
      }
  }
  post_define_dots(hscale,vscale,big_dots_flag,dot_magnifier,post_file,
               zoom_labels_count,stripe_flag,last_structure);
  if(forced_stripe_flag)
     stripe_flag=TRUE;
  else
   {
     if((post_file==MAIN_POST)||(last_structure<4))
       {stripe_flag=FALSE;
       }
     else
      {if(!nongray_for_3_flag)
        { if(  ((540.0/(display_b-display_t)/(last_structure-1)<STRIPE_PS) ) ||
          ((540.0/(display_r-display_l)/(last_structure-1.)<STRIPE_PS) )  )
            {stripe_flag=FALSE;
             }
          else
           {stripe_flag=TRUE;
           }
        }
       else
        stripe_flag=FALSE;
      }
  }
  fprintf(postfp,"\n/len { %d } def",display_b-display_t+1);
  fprintf(postfp,"\n\n%s Set background color\n", "%%");
  fixcolor_plot_postscript(COL_BACKGROUND);
  fprintf(postfp,"\n0 0 moveto");
  fprintf(postfp,"\n611 0 rlineto");
  fprintf(postfp,"\n0 791 rlineto");
  fprintf(postfp,"\n-611 0 rlineto");
  fprintf(postfp,"\n0 -791 rlineto");
  fprintf(postfp,"\nclosepath");
  fprintf(postfp,"\nfill");
  fprintf(postfp,"\n\n%s Display text \n", "%%");
  fixcolor_plot_postscript(COL_TEXT);
  display_copyright_message_ps(postfp,TRUE,"ct_boxplot");
  /* center the title */
  fprintf(postfp,"\n /Helvetica findfont");
  fprintf(postfp,"\n 20 scalefont");
  fprintf(postfp,"\n setfont");
  center=72*8.5/2.-(strlen(sequence_name)+23)/2.*20.*7./16.;
  if(center<5)
     center=5.;
  fprintf(postfp,"\n %f 720 moveto",center);
  fprintf(postfp,"\n (Structure dot plot for %s) show",sequence_name);
  fprintf(postfp,"\n /Helvetica findfont");
  fprintf(postfp,"\n 10 scalefont");
  fprintf(postfp,"\n setfont");
  fixcolor_plot_postscript(COL_TEXT);
  fprintf(postfp,"\n\n%s Translate for box coordinates \n", "%%");
  fprintf(postfp,"\n 25 136 translate");
 /* label top ticks */
  fprintf(postfp,"\n\n %f setlinewidth",0.7/((hscale+vscale)/2.0)); 
  fprintf(postfp,"\n\n%s Horizontal Tic marks \n", "%%");
  plot_hor_tics_zoom_post(hscale,grid_flag,display_l,display_r);
  /* label side ticks */ 
   fprintf(postfp,"\n\n%s Vertical Tic Marks \n", "%%");
   plot_ver_tics_zoom_post(display_t,display_b,vscale,grid_flag);
   fprintf(postfp,"\n\n%s Draw box \n", "%%");
   fprintf(postfp,"\n\n %f setlinewidth",0.7);
/* draw box                left and bottom */
     fprintf(postfp,"\n 0 540.5 moveto 0 0 lineto stroke");
     fprintf(postfp,"\n 0 0.0 moveto 540.5 0 lineto stroke");    
  /* draw right and top */
  fprintf(postfp,"\n 0 540.5 moveto 540.5 540.5 lineto stroke");
  fprintf(postfp,"\n 540.5 0  moveto 540.5 540.5 lineto stroke");
/* set scale */
 fprintf(postfp,"\n\n%s Scale to box size coordinates \n", "%%");
  fprintf(postfp,"\n %f %f scale",hscale,vscale);
  fprintf(postfp,"\n\n %f setlinewidth",0.7/((hscale+vscale)/2.0));
  fprintf(postfp,"\n\n%s Diagonal line \n", "%%");
  draw_diagonal_line_post(display_l,display_r,display_t,
            display_b); 
  fprintf(postfp,"\n\n%s Draw dots \n", "%%");
  draw_postscript_dots(counter,stripe_flag,
       last_structure,start,count,
        diag,display_l,display_r,display_t,display_b,lower_structure,
        nongray_for_3_flag,gray3_content,gray3);
  if((post_file!=MAIN_POST)&&(zoom_labels_count>0))
    {fprintf(postfp,"\n\n%s Draw Labels \n", "%%");
     fixcolor_plot_postscript(COL_DOTS);
     fprintf(postfp,"\n /Helvetica findfont");
     fprintf(postfp,"\n 8 scalefont");
     fprintf(postfp,"\n setfont");
      for(i=0;i<zoom_labels_count;i++)
      {fprintf(postfp,"\n%d %d sbox",
              zoom_labels_row[i]-display_t+1,
              zoom_labels_col[i]-display_l+1);
        fprintf(postfp,"\n %d %d moveto",
              zoom_labels_col[i]-display_l+2,
              display_b-zoom_labels_row[i]+1);
       fprintf(postfp,"\n %f %f scale",1/hscale,1/vscale);
       fprintf(postfp,"\n 0 -1.5 rmoveto");
       fprintf(postfp,"\n (%s) show ",zoom_labels_string[i]);
       fprintf(postfp,"\n %f %f scale",hscale,vscale);
      }
    }
  fprintf(postfp,"\n %f %f scale",1./hscale,1./vscale); 
  fprintf(postfp,"\n -25 -136 translate");
  if((nongray_for_3_flag)||(lower_structure<=-10))
       {fixcolor_plot_postscript(5);
        fprintf(postfp,"\n 450 75 moveto");
        fprintf(postfp,"\n ( 1 & 2       %d) show",gray3[0]);
        fixcolor_plot_postscript(6);
        fprintf(postfp,"\n 450 60 moveto");
        fprintf(postfp,"\n ( 1 & 3       %d) show",gray3[1]);
        fixcolor_plot_postscript(7); 
        fprintf(postfp,"\n 450 45 moveto");
        fprintf(postfp,"\n ( 2 & 3       %d) show",gray3[2]);        
       }
  fixcolor_plot_postscript(COL_COMPLETE_OVERLAP);
   fprintf(postfp,"\n 35 120  moveto ");
   fprintf(postfp,"\n (Full Overlap) show");
   fprintf(postfp,"\n 125 120  moveto ");
   fprintf(postfp,"\n (%d) show",counter[0]);
   if(!stripe_flag)
    {
     if(last_structure>3)
       {fixcolor_plot_postscript(COL_PART_OVERLAP);
        fprintf(postfp,"\n 35 105 moveto ");
        fprintf(postfp," (Partial Overlap) show"); 
        fprintf(postfp,"\n 125 105 moveto ");
        fprintf(postfp," (%d) show",counter[1]);
       }
    }
  stop_point=last_structure/2+1;
  if(last_structure >16)
     {stop_point=9;
      other_stop_point=15;
      fixcolor_plot_postscript(16);
      fprintf(postfp,"\n 450 20 moveto");
      fprintf(postfp,"\n (%d Others ) show",last_structure-15);
     }
  else
    other_stop_point=last_structure;
  for(color=2,offset=0;color<=stop_point;
             color++,offset=offset+15)
         {fixcolor_plot_postscript(color);
          fprintf(postfp,"\n 200 %d moveto",120-offset);
          fprintf(postfp,"\n (%s) show",file_data_ct[color-2]);
          fprintf(postfp,"\n 350 %d moveto",120-offset);
          fprintf(postfp,"\n (%d*) show",counter[color]);
        }
      for(color=(stop_point+1),offset=0;color<=other_stop_point;
             color++,offset=offset+15)
         {fixcolor_plot_postscript(color);
          fprintf(postfp,"\n 400 %d moveto",120-offset);
          fprintf(postfp,"\n (%s) show",file_data_ct[color-2]);
          fprintf(postfp,"\n 550 %d moveto",120-offset);
          fprintf(postfp,"\n (%d*) show",counter[color]);
        }
  fixcolor_plot_postscript(COL_TEXT);
  fprintf(postfp,"\n 30 5 moveto");
  fprintf(postfp,"\n (*Counts for each structure include overlap dots.) show");
  fprintf(postfp,"\n showpage\n");
  fprintf(postfp,"%sEOF\n","%%");
  fclose(postfp);
}  
void finish_making_post(int big_post_dots,
          int post_file, char *post_filename,
          int zoom_labels_count,int display_l,
          int display_r,int display_t,int display_b,float color_table[][3],
          int *zoom_labels_row,int *zoom_labels_col,char *zoom_labels_string[],
          int grid_flag,int last_structure, int *start, int *count,
          struct helix *diag,int len,
          int lower_structure,char *sequence_name,
          char file_data_ct[][80],float dot_magnifier,
          int forced_stripe_flag,int nongray_for_3_flag,
          int *gray3_content)
{/* post_file is 1 for main, 2 for zoom */
   int error;
   printf("\n The name for the postscript file is     %s \n",
              post_filename);
   error=open_post(post_filename);
   if(error)
     {printf("\n Error creating postscript file \n");
      return;
     }
   general_post(big_post_dots,post_file,zoom_labels_count,
       display_l,display_r,display_t,display_b,color_table,
       zoom_labels_row,zoom_labels_col,zoom_labels_string,grid_flag,
        last_structure,start,count, diag,
        len,lower_structure,sequence_name,file_data_ct,dot_magnifier,
        forced_stripe_flag,nongray_for_3_flag,gray3_content);
}


