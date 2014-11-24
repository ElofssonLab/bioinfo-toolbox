/* Darrin Stewart and Michael Zuker */
/* Feb 4, 1999  */
/* ct_boxplot_general.c */
/* find distance between tick marks */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ct_boxplot.h"



/* strip off .ct from end of title */
void fix_title(char *title)
{int i;
   i=strlen(title);
   if(i>5)
     {i=i-3;
      if((title[i]=='.')&&(title[i+1]=='c')&&(title[i+2]=='t'))
          {title[i]='\0';
           return;
          }
     } 
    return;
}


/* Find the distance in bases between tick marks */
int step_fun(int dis_l,int dis_r)
{ int step;
  step=(dis_r-dis_l+1)/8; /* make 8 steps  (9 tick marks total) */
  if(step<40)
    { if(step<3)
        step=2;
       else if(step>5)
          {step=(step+5)/10; /* make a large step a multiple of ten */
           step=step*10;
          }
         else
            if(step==4) /* make steps of 4 look like 5 */
               step=5;
    }
  else if(step>90)
          {step=(step+50)/100; /* round to nearest 100 */
           step=step*100;
          }
  else
       {step=(step+25)/50; /* round to nearest 50 */
        step=step*50;
       }
  return step;
}

/* make the tick marks start at a nice place */
/* ends are treated separately and always drawn */
int start_fun(int dis_l,int step)
{int start;
  start=dis_l+step;
  if(step>90)
    {start=(start+60)/100; /* start at a multiple of 100 */
     start=start*100;
    }
  else if(step>=40)
    {start=(start+32)/50; /* start at a multiple of 50 */
     start=start*50;
    }
  else if(step>5)
   {start=(start+5)/10; /* start at a multiple of ten */
    start=start*10;
   }
  return start;
}

/* Conversion Routines to display integers and floats */

char *num_string_int(int x,char *str) /* converts integer to string */
{    sprintf(str,"%d",x);
     return str;
}


/* given a row and column, return the color of that position */
int color_read(int row,int col,int color,int *g_diag_start, 
      struct helix *g_diag)
{int diag,diag_start,diag_end,i,temp_row;
 if(row>col)
   {temp_row=row;
    row=col;
    col=temp_row;
   }
 diag=col+row-1;
 if(diag<0)
  return FALSE;
 /* printf("\n color read entered with row=%d col=%d color=%d",row,col,color);*/
 diag_start=g_diag_start[diag];
 diag_end=g_diag_start[diag+1];
 /* printf("\n diag_start is %d diag_end is %d",diag_start,diag_end); */
 /* Return true if this sequence has a a dot at this point */
 /* Must match color, or use -1 to indicate any match is ok */
 for(i=diag_start;i<diag_end;i++)
     {temp_row=g_diag[i].row;
     /* printf("\n temp_row is %d col is %d length is %d", temp_row,
         g_diag[i].row,g_diag[i].length); */
      if(temp_row>row) 
        return FALSE;   /*  entries in diagonal are sorted */
      if((temp_row<=row)&&(g_diag[i].row+g_diag[i].length-1>=row)
          &&((g_diag[i].color==color)||(color==-1)))
        {/*printf("\n returning TRUE");*/ 
         return TRUE;
        }
     }
 return FALSE;
} 

int *store_gray_content(int *gray3_content,int helices_in_plot_file,int len,
       struct helix *diag, int *diag_start)
{int new_color,row,col,position;
 int new_row;
  gray3_content=(int *)malloc((len+1)*sizeof(int));
  if(gray3_content==NULL)
   {printf("\n Error allocating memory for storage of gray_content\n");
    exit(1);
   }
 for(row=1;row<=len;row++)
    gray3_content[row]=0; /* set to zero first */
 for(position=0;position<helices_in_plot_file;position++)
     {new_row=diag[position].row;
      col=diag[position].column;
      new_color=diag[position].color;
      if(new_color==1) /* look for gray */
            { 
              for(row=new_row;row<new_row+diag[position].length;
                                row++)
                     {/* if(seq 4 and 2  are here store 5) */
                      /* if(seq 4 and 3 are here store 6) */
                      /* else store 1 */  
                       if(color_read(row,col,4,diag_start,diag))
                           {if(color_read(row,col,2,diag_start,diag))
                              {gray3_content[row]=-1*col;
                             
                              }
                            else
                               {gray3_content[row]=col;
                             
                               }
                           }
                     col--;
		     }
	  }
		                               
     }

 return gray3_content;
}

