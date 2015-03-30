/* boxplot_ps.inc */
/* by Darrin Stewart and Michael Zuker */
/* May 17, 2000 */


FILE *postfp;  /* for postscript file */

int global_postscript_points_plotted;



void fixcolor_plot_postscript(int color)
{ /*  postscript colors */
  if((color>=0)&&(color<TOTAL_COLORS))
     fprintf(postfp,"\n  color_%d\n",color); 
}

void fixcolor_plot_postscript_optimal(void)
{ /*  postscript optimal */
     fprintf(postfp,"\n  color_OP\n"); 
}

void fixcolor_plot_postscript_define(void) 
{ /*  postscript colors */
  int color;
  fprintf(postfp,"\n%s Set colors","%%");
  fprintf(postfp,"\n%s 0 is for text, 1 is optimal dots","%%");
  fprintf(postfp,"\n%s 2 is almost optimal dots, 3-8 higher energy dots ","%%");
  fprintf(postfp,"\n%s 9 is Labels, 10 is background, 11 is grid lines","%%");
  fprintf(postfp,"\n%s OP is for lower triangle","%%");
  fprintf(postfp,"\n%s original colors from boxplot-3.0.col\n\n","%%");

  
  for(color=0;color<TOTAL_COLORS;color++)
    {if((global_number_of_colors>=7)&&(((color>=1)&&(color<=LAST_PROB_COLOR))))
        {fprintf(postfp,"/color_%d {  %.3f %.3f %.3f setrgbcolor} def \n",
             color,
             global_prob_color[color][0],
             global_prob_color[color][1],
             global_prob_color[color][2]);
        }
     else
       fprintf(postfp,"/color_%d {  %.3f %.3f %.3f setrgbcolor} def \n",color,
         global_color[color][0],global_color[color][1],global_color[color][2]);
    }
 fprintf(postfp,"/color_OP {  %.3f %.3f %.3f setrgbcolor} def \n",
         global_color[COLOR_OPTIMAL][0],
         global_color[COLOR_OPTIMAL][1],
         global_color[COLOR_OPTIMAL][2]);
 
}




void single_hor_tic_post(int j_col,float scale,int make_grid_flag)
{
/* offset=1.*(global_length+2)/540.;  move boundary 1 points over */
 fprintf(postfp,"\n %f 1. scale",scale);
 fprintf(postfp,"\n\n %f setlinewidth",0.40/scale);
 fprintf(postfp,"\n %d 541 moveto",j_col+1);
 fprintf(postfp,"\n 0 8 rlineto");
 if(make_grid_flag==TRUE)
   {fprintf(postfp,"\n stroke");
    fixcolor_plot_postscript(COLOR_GRID);
    fprintf(postfp,"\n %d 539 moveto",j_col+1);
    fprintf(postfp,"\n 0 %f rlineto",-1.*(global_length+1)*scale);
    fprintf(postfp,"\n stroke");
    fixcolor_plot_postscript(COLOR_TEXT);
   }
 fprintf(postfp,"\n %d 550 moveto",j_col+1);
 fprintf(postfp,"\n %f 1. scale",1./scale);
 fprintf(postfp,"\n -6 0 rmoveto");
 fprintf(postfp,"\n (%d) show",j_col);
}

void plot_hor_tics_post(float scale,int make_grid_flag)
{ 
  int j_col,step,start;
  step=step_fun(1,global_length);
  start=start_fun(1,step);
  single_hor_tic_post(global_length,scale,FALSE);
  single_hor_tic_post(1,scale,FALSE);
  for(j_col=start;j_col<(global_length-2*step/3);j_col=j_col+step)
      { single_hor_tic_post(j_col,scale,make_grid_flag);
      }

}

void single_ver_tic_post(int i_row,float scale,int make_grid_flag)
{
/*  offset=1.*(global_length+2)/540.;  move boundary 1 points over */
  fprintf(postfp,"\n 1.0 %f scale",scale);
 fprintf(postfp,"\n\n %f setlinewidth",0.40/scale);
  fprintf(postfp,"\n 541.5 %d moveto",global_length-i_row+2);
  fprintf(postfp,"\n 8 0 rlineto");
  if(make_grid_flag==TRUE)
   {fprintf(postfp,"\n stroke");
    fixcolor_plot_postscript(COLOR_GRID);
    fprintf(postfp,"\n 539.5 %d moveto",global_length-i_row+2);
    fprintf(postfp,"\n %f 0 rlineto",-1.0*scale*(global_length+1));
    fprintf(postfp,"\n stroke");
    fixcolor_plot_postscript(COLOR_TEXT);
   }
  fprintf(postfp,"\n 550.5 %d moveto",global_length-i_row+2);
  fprintf(postfp,"\n 1.0 %f scale",1./scale);
  fprintf(postfp,"\n (%d) show",i_row);
}

void plot_ver_tics_post(float scale,int make_grid_flag)
{ int i_row,step,start;
 step=step_fun(1,global_length);
 start=start_fun(1,step);
 single_ver_tic_post(1,scale,FALSE);
 single_ver_tic_post(global_length,scale,FALSE);
 for(i_row=start;i_row<(global_length-step/2);i_row=i_row+step)
     {single_ver_tic_post(i_row,scale,make_grid_flag);
      }
}

void single_hor_tic_zoom_post(int j_col,float scale,int make_grid_flag)
{ fprintf(postfp,"\n %f 1. scale",scale);
 fprintf(postfp,"\n\n %f setlinewidth",0.40/scale);
  fprintf(postfp,"\n %d 540 moveto",j_col-display_l+1);
  fprintf(postfp,"\n 0 8 rlineto");
  if(make_grid_flag==TRUE)
   {fprintf(postfp,"\n stroke");
    fixcolor_plot_postscript(COLOR_GRID);
    fprintf(postfp,"\n %d 539 moveto",j_col-display_l+1);
    fprintf(postfp,"\n 0 -538 rlineto");
    fprintf(postfp,"\n stroke");
    fixcolor_plot_postscript(COLOR_TEXT);
   }
  fprintf(postfp,"\n %d 550 moveto",j_col-display_l+1);
  fprintf(postfp,"\n %f 1. scale",1./scale);
  fprintf(postfp,"\n -4 0 rmoveto");
  fprintf(postfp,"\n (%d) show",j_col);
}

void plot_hor_tics_zoom_post(float scale,int make_grid_flag)
{ 
  int j_col,step;
  int start;
  step=step_fun(display_l,display_r);
  start=start_fun(display_l,step);
  single_hor_tic_zoom_post(display_r,scale,FALSE);
  single_hor_tic_zoom_post(display_l,scale,FALSE);
  for(j_col=start;j_col<(display_r-2*step/3);
                        j_col=j_col+step)
      {single_hor_tic_zoom_post(j_col,scale,make_grid_flag);
      }
}

void single_ver_tic_zoom_post(int i_row,float scale,int make_grid_flag)
{ fprintf(postfp,"\n 1. %f scale",scale);
 fprintf(postfp,"\n\n %f setlinewidth",0.40/scale);
  fprintf(postfp,"\n 540.5 %d moveto",display_b-i_row+1);
  fprintf(postfp,"\n 8 0 rlineto");
  if(make_grid_flag==TRUE)
   {fprintf(postfp,"\n stroke");
    fixcolor_plot_postscript(COLOR_GRID);
    fprintf(postfp,"\n 539.5 %d moveto",display_b-i_row+1);
    fprintf(postfp,"\n -538 0 rlineto");
    fprintf(postfp,"\n stroke");
    fixcolor_plot_postscript(COLOR_TEXT);
   }
  fprintf(postfp,"\n 550.5 %d moveto",display_b-i_row+1);
  fprintf(postfp,"\n 1. %f scale",1./scale);
  fprintf(postfp,"\n (%d ) show",i_row);
}

void plot_ver_tics_zoom_post(float scale,int make_grid_flag)
{ int i_row,step,start;
 step=step_fun(display_t,display_b);
 start=start_fun(display_t,step);
 single_ver_tic_zoom_post(display_t,scale,FALSE);
 single_ver_tic_zoom_post(display_b,scale,FALSE);
 for(i_row=start;i_row<(display_b-step/2);i_row=i_row+step)
   {single_ver_tic_zoom_post(i_row,scale,make_grid_flag);  
   }
}

void draw_diagonal_line_post(int ps_file)
{float right,left,top,bottom;
 int dis_l,dis_r,dis_t,dis_b;
 float offset;
 if(ps_file!=MAIN_POST)
   {
     if(!((display_r>=display_t)&&(display_b>=display_l)))
             /* cases all below, all above */
       return; /* There is no diagonal line */
     dis_l=display_l;
     dis_r=display_r;
     dis_t=display_t;
     dis_b=display_b;
   }
 else
   {dis_l=1;
    dis_r=global_length;
    dis_t=1;
    dis_b=global_length;
   }

 /* find right end of line */
 if(dis_r>dis_b)
   { /* diagonal crosses bottom border */
     bottom=(float)dis_b;
     right=(float)bottom;
   }
 else
   { /* diagonal crosses right border */
     right=(float)dis_r;
     bottom=(float)right;
   }
 /* find left end of line */
 if(dis_l>dis_t)
   {/* diagonal crosses left side */
    left=(float)dis_l;
    top=(float)dis_l;
   }
 else
   {/* diagonal crosses top side */
    top=(float)dis_t;
    left=(float)top;
   }
 if(ps_file!=MAIN_POST)
   {
    fprintf(postfp,"\n %f  %f moveto",left-dis_l,dis_b-top+2.0);
    fprintf(postfp,"\n %f %f  lineto stroke",right-dis_l+2,dis_b-bottom);
    /* also draw left and bottom for main window */

   }
 else
   {offset=1.*(global_length+2)/540.; /* move boundary 1 points over */
    fprintf(postfp,"\n %f  %f moveto",left-dis_l+1-offset,global_length-top+3+offset);
    fprintf(postfp,"\n %f %f  lineto stroke",right-dis_l+3+offset,global_length-bottom+1-offset
);
   /* also draw left and bottom for main window */ 
   fprintf(postfp,"\n %f %f  moveto ",right-dis_l+3+offset,global_length-bottom+1-offset);
   fprintf(postfp,"\n %f %f lineto",left-dis_l+1-offset,global_length-bottom+1-offset);
   fprintf(postfp,"\n %f %f lineto",left-dis_l+1-offset,global_length-top+3+offset);
   fprintf(postfp,"\n %f %f lineto",right-dis_l+3+offset,global_length-top+3+offset);
   fprintf(postfp,"\n %f %f lineto closepath stroke",right-dis_l+3+offset,global_length-bottom+1-offset);
  
   }
} 


void draw_postscript_dots_main_opt_prob(void)
{int row,col,new_row;
 int location;
 int color; 
 color=1;
 fixcolor_plot_postscript_optimal();
 for(location=0;location<global_helices_in_plot_file;location++)
       {  if(global_diag[location].energy<=-20001)
	       {new_row=global_diag[location].row;
                col=global_diag[location].column;              
                for(row=new_row;
                    row<(new_row+global_diag[location].length);row++)
		      {
                        fprintf(postfp,"\n%d %d b%d",col-1,row+1,color);
                         col--;
                      } 
		                           

	       }
	       

       }
}

void draw_postscript_dots_main(void)
{int row,col,new_row;
 int location,new_color;
 int color; 
 global_postscript_points_plotted=0;
 for(color=global_number_of_colors;color>=1;color--)
 {fixcolor_plot_postscript(color);
  for(location=0;location<global_helices_in_plot_file;location++)
       {  if((global_diag[location].energy<=global_energy_cutoff)&&
          ((global_diag[location].length>=global_chain_len)||
           (global_diag[location].energy==global_optimal_energy)))

  
	  {new_color=global_diag[location].color;
           if(new_color==color)
	       {new_row=global_diag[location].row;
                col=global_diag[location].column;
                if((!global_opt_prob_flag)
                      ||(global_diag[location].energy>=-20000))
		    {
                for(row=new_row;
                    row<(new_row+global_diag[location].length);row++)
		      {fprintf(postfp,"\n%d %d b%d",row-1,col+1,color);
                       global_postscript_points_plotted++;
                       if(!global_opt_prob_flag)
		        {
                         if((new_color==1)&&(!global_adjust_col_flag))
                          {fixcolor_plot_postscript_optimal();
                           fprintf(postfp,"\n%d %d b%d",col-1,row+1,color);
                           fixcolor_plot_postscript(color);
                          }
		        }
                       col--;
                      } 
		    }                       

	       }
          }
  
       }
 }
if(global_opt_prob_flag)
   draw_postscript_dots_main_opt_prob();
}

void draw_postscript_dots_zoom(void)
{int row,col,new_row,count;
 int diag,position,location,new_color;
 int start_diag,end_diag;
 int color; 
 float size;
 global_postscript_points_plotted=0;
 start_diag=display_t+display_l-1;
 end_diag=display_b+display_r-1;
 for(color=global_number_of_colors;color>=1;color--)
 {/* printf("\n color is %d",color); */
  fixcolor_plot_postscript(color);
  for(diag=start_diag;diag<=end_diag;diag++)
    {/* printf("\n diagonal is %d",diag); */
     position=global_diag_start[diag];
     for(count=0;( (count<global_diag_count[diag])
         &&(global_diag[position+count].row<=display_b) ) ;count++)
       {location=position+count;
        if((!global_opt_prob_flag)
                      ||(global_diag[location].energy>=-20000))
		    {
        if( (global_diag[location].energy==global_optimal_energy)||
          (global_diag[location].energy<=global_energy_cutoff)&&
          (global_diag[location].length>=global_chain_len))
	  {new_color=global_diag[location].color;
             if(new_color==color)
	       {new_row=global_diag[location].row;
                col=global_diag[location].column;
                for(row=new_row;
                    row<(new_row+global_diag[location].length);row++)
		      {if((row>=display_t)&&(row<=display_b)
                          &&(col>=display_l)&&(col<= display_r))
                        { global_postscript_points_plotted++;
                           if(global_contains_prob!=TRUE)
                            fprintf(postfp,"\n%d %d b%d",
                               row-display_t+1,col-display_l+1,color);
                          else
			    {size=sqrt((-1.0*(float)global_diag[location].energy)/10000.);
                             if(size>1.2)
                                size=1.2;
                             fprintf(postfp,"\n%d %d %.5f pb",
                                  row-display_t+1,col-display_l+1,size);
                            }
                        }
                          col--;
                      }                        

	       }
          }
		    }
        }
    }
 }
 /* draw optimal points below main diagonal */
 color=1;
 fixcolor_plot_postscript_optimal();
 for(diag=start_diag;diag<=end_diag;diag++)
    {position=global_diag_start[diag];
     for(count=0;( (count<global_diag_count[diag])
         &&(global_diag[position+count].row<=display_r) ) ;count++)
       {location=position+count;
	new_color=global_diag[location].color;
             if(new_color==color)
	       {new_row=global_diag[location].row;
                col=global_diag[location].column;
                if(col>=display_t)
		  {
                for(row=new_row;
                    row<(new_row+global_diag[location].length);row++)
		      {if((row>=display_l)&&(row<=display_r)
                          &&(col>=display_t)&&(col<=display_b))
                        {
                          if(global_contains_prob!=TRUE)
                            fprintf(postfp,"\n%d %d b%d",
                               col-display_t+1,row-display_l+1,color);
                          else
			    {if(global_opt_prob_flag)
                              {if(global_diag[location].energy<=-20001)
                                fprintf(postfp,"\n%d %d %.5f pb",
                                  col-display_t+1,row-display_l+1,1.0);
                              }
                            else
                              {size=sqrt((-1.0
                                *(float)global_diag[location].energy)/10000.);
                              if(size>1.2)
                                size=1.2;
                              fprintf(postfp,"\n%d %d %.5f pb",
                                  col-display_t+1,row-display_l+1,size);
                              }
                            }
                         
                      }
                    col--;                        	
	           }
                 }
  
        }}
    }


}






void display_time(void)
{time_t now;
 now=time(NULL);
 fprintf(postfp,"\n%s %s","%%",ctime(&now));
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


void extra_ps_char(void)
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
  fprintf(postfp,"\n/delta { /Symbol findfont findISO { reencodeISO /Symbol-ISO exch definefont");  
  fprintf(postfp,"\n                              sf scalefont setfont }");  
  fprintf(postfp,"\n                           { sf scalefont setfont } ifelse");  
  fprintf(postfp,"\n  (d) show } def");  
  fprintf(postfp,"\n");
}


void clean_title_ps(char *title)
{/* place backslash in front of ( or ) or " or ' for postscript file */
 char clean_title[130];
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


void display_title_message_ps(void)

{ char year[8];
  float ver_offset;
  ver_offset=700.;
  display_time_year(year);
  /* Try not to make text overwrite image below */
  fprintf(postfp,"\n color_%d\n",COLOR_COPY);
  fprintf(postfp,"/sf 7.5 def\n");
  fprintf(postfp,"/Helvetica findfont findISO { reencodeISO /Symbol-ISO exch definefont \n");
  fprintf(postfp,"               sf scalefont setfont }\n");
  fprintf(postfp,"               { sf scalefont setfont } ifelse\n");
  fprintf(postfp,"20 %f moveto\n",ver_offset+66.);
  fprintf(postfp," (%s) show\n",copyright_message1);
  fprintf(postfp,"26 %f moveto\n",ver_offset+55);
  fprintf(postfp," (\251 %s %s) show\n",year,copyright_message2);

}



void post_define_dots(float hscale,float vscale,int post_adjust,
   int make_label_flag,int zoom_labels_count)
{float dot_size,new_scale;
 int i;
  fprintf(postfp,"\n%s   Define dots \n","%%"); 
/* define all dot sizes s */
 fprintf(postfp,
       "\n\n%s Dot shape and Size definitions\n","%%");
 if(global_post_file==MAIN_POST)
    zoom_labels_count=0;
 if((global_post_file==MAIN_POST)||(!global_contains_prob))
   {
    fprintf(postfp,"%s s1 is Optimal, s2 next to optimal through s7\n","%%");
    fprintf(postfp,"%s sbox is for labels\n","%%");
    fprintf(postfp,"\n%s Change si below to adjust size of dots\n","%%");
    fprintf(postfp,"%s 1.0 is maximum size without overlapping dots \n","%%");
    fprintf(postfp,"%s sl is for labels \n","%%");
    fprintf(postfp,"%s s1 = 2.0 doubles size of optimal dots, b1\n","%%");
    fprintf(postfp,"%s s2 = 2.0 doubles size of near optimal dots dots, b2\n",
      "%%");
    fprintf(postfp,"%s Use %.4f to make dots 1/72 inch or 1 point\n","%%",
       2.0/(hscale+vscale));
    fprintf(postfp,"%s Ignore Si,NSi,len \n\n","%%");
   }
  if((global_contains_prob)&&(global_post_file!=MAIN_POST))
   {fprintf(postfp,"%s sp is size of all dots\n","%%");
    fprintf(postfp,"%s prior to scaling base on probability\n","%%");
    fprintf(postfp,"%s sp = 2.0 doubles size of dots\n","%%");
    fprintf(postfp,"%s Use %.4f to make dots 1/72 inch or 1 point\n","%%",
       2.0/(hscale+vscale));
    fprintf(postfp,"%s Ignore SP,NSP,len \n\n","%%");
   }
 if((make_label_flag)||(zoom_labels_count>0))
   {fprintf(postfp,"/sl { .5 } def\n");
   }
 
 if((global_contains_prob)&&(global_post_file!=MAIN_POST))
   {if(post_adjust)
       {new_scale=(hscale+vscale)/2.0;
        dot_size=global_ps_dot_min[1][0];
        if(new_scale*dot_size<global_ps_dot_min[1][1])
           dot_size=global_ps_dot_min[1][1]/new_scale;
       }
      else
       dot_size=1.0;
     fprintf(postfp,"/sp { %.4f } def\n\n",dot_size);
     fprintf(postfp,"/SP { .5 sp mul } def\n");
     fprintf(postfp,"/NSP {-1.0 SP mul } def\n");
   }
 if((global_post_file==MAIN_POST)||(!global_contains_prob))
   {for(i=1;i<=global_number_of_colors;i++)
       {if(post_adjust)
         {new_scale=(hscale+vscale)/2.0;
          dot_size=global_ps_dot_min[i][0];
          if(new_scale*dot_size<global_ps_dot_min[i][1])
             dot_size=global_ps_dot_min[i][1]/new_scale;
          }
        else
          dot_size=1.0;
        fprintf(postfp,"/s%d { %.4f } def\n",i,dot_size);
       }
   }
 /* define S and NS */
 fprintf(postfp,"\n\n");
 if((make_label_flag)||(zoom_labels_count>0)) 
  {fprintf(postfp,"/Sl {.50 sl mul } def\n");
    fprintf(postfp,"/NSl {-1.0 Sl mul } def\n");
   }
 if((global_post_file==MAIN_POST)||(!global_contains_prob))
   { for(i=1;i<=global_number_of_colors;i++) 
       {fprintf(postfp,"/S%d {.50 s%d mul } def\n",i,i);
        fprintf(postfp,"/NS%d {-1.0 S%d mul } def\n",i,i);
       }
   }
 /* define details of box */
 if((make_label_flag)||(zoom_labels_count>0))
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
  if((global_post_file==MAIN_POST)||(!global_contains_prob))
   {for(i=1;i<=global_number_of_colors;i++)
     {if(post_adjust)
       {new_scale=(hscale+vscale)/2.0;
        dot_size=global_ps_dot_min[i][0];
        if(new_scale*dot_size<global_ps_dot_min[i][1])
           dot_size=global_ps_dot_min[i][1]/new_scale;
        dot_size=dot_size/2.0;
       }
      else
       dot_size=.5;
      fprintf(postfp,"\n/b%d {",i); /* for non-probability or main window */
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
   }
  if((global_contains_prob)&&(global_post_file!=MAIN_POST))
    {fprintf(postfp,"\n/pb {");
          /* for boxes whose size depends on probability */
     fprintf(postfp,"\n /s exch def");
     fprintf(postfp,"\n exch len exch sub 1 add bx");
     fprintf(postfp,"\n} def");
     fprintf(postfp,"\n/bx {");
     fprintf(postfp,"\n newpath");
     fprintf(postfp,"\n moveto");
     fprintf(postfp,"\n SP s mul  0 rlineto");
     fprintf(postfp,"\n 0 SP  s mul  rlineto");
     fprintf(postfp,"\n NSP s mul   0 rlineto");
     fprintf(postfp,"\n NSP s mul   0 rlineto"); 
     fprintf(postfp,"\n 0     NSP s mul rlineto");
     fprintf(postfp,"\n 0     NSP s mul rlineto");
     fprintf(postfp,"\n SP s mul     0 rlineto");
     fprintf(postfp,"\n SP s mul     0 rlineto");
     fprintf(postfp,"\n 0   SP s mul rlineto");
     fprintf(postfp,"\n closepath");
     fprintf(postfp,"\n fill");
     fprintf(postfp,"\n stroke");
     fprintf(postfp,"\n} def");
    }
}

void draw_single_label(int label_i,int label_j,float hscale,float vscale)
{   fprintf(postfp,"\n%d %d sbox",label_i-display_t+1,
               label_j-display_l+1);
    fprintf(postfp,"\n %d %d moveto",label_j-display_l+2,
              display_b-label_i+1);
    fprintf(postfp,"\n %f %f scale",1/hscale,1/vscale);
    fprintf(postfp,"\n 0 -1.5 rmoveto");
    fprintf(postfp,"\n (\\(%d,%d\\)) show ",label_i,label_j);
    fprintf(postfp,"\n %f %f scale",hscale,vscale);
}


void general_post(char *title,int make_grid_flag,int make_label_flag,
                 int label_i,int label_j,int post_adjust,
                 int zoom_labels_count,
                 int *zoom_labels_row,
                 int *zoom_labels_col)
{/* make label refers to the single label possible from boxplot_ng.c */
 /* zoom_labels_count can be 0 or larger and comes from boxplot.c */
 float hscale,vscale,x_float;
  char str[100];
  int color,offset;
  float center;
  int label_count;
  int name_len;
  char post_energy_file[120];
  if(make_label_flag)
    {strcpy(post_energy_file,global_post_filename);
     name_len=strlen(post_energy_file);
     post_energy_file[name_len-3]='\0';
     make_energy_file(post_energy_file,label_i,label_j);
    }
  if(global_post_file==MAIN_POST) 
     {zoom_labels_count=0;
       hscale=540.0/(global_length+2);
         /* 540 is 72*7.5 for 7.5 inches of paper */
      vscale=hscale;
    }
  else
    {vscale=540.0/(display_h+1);
     hscale=540.0/(display_w+1);
    }
  fprintf(postfp,"%c!PS-Adobe-3.0 EPSF-2.0\n",'%');
  fprintf(postfp,"%sBoundingBox: 0 0 612 792\n","%%");
  fprintf(postfp,"\n%s Created by boxplot_ng version 2.0","%%");
  fprintf(postfp,"\n%s Darrin Stewart","%%");
  fprintf(postfp,"\n%s  and Michael Zuker","%%");
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
  fixcolor_plot_postscript_define(); 
  if(global_post_file==MAIN_POST)
    fprintf(postfp,"\n/len { %d } def",global_length);
  else
    fprintf(postfp,"\n/len { %d } def",display_h);
  extra_ps_char();
  post_define_dots(hscale,vscale,post_adjust,make_label_flag,zoom_labels_count);
  fprintf(postfp,"\n\n%s Set background color\n", "%%");
  fixcolor_plot_postscript(COLOR_BACKGROUND);
  fprintf(postfp,"\n0 0 moveto");
  fprintf(postfp,"\n611 0 rlineto");
  fprintf(postfp,"\n0 791 rlineto");
  fprintf(postfp,"\n-611 0 rlineto");
  fprintf(postfp,"\n0 -791 rlineto");
  fprintf(postfp,"\nclosepath");
  fprintf(postfp,"\nfill");
  fprintf(postfp,"\n\n%s   Draw Text \n","%%"); 
  display_title_message_ps();
  fixcolor_plot_postscript(COLOR_TEXT);
  fprintf(postfp,"\n/sf 20 def");
  fprintf(postfp,"\n /Helvetica findfont findISO { reencodeISO /Symbol-ISO exch definefont ");
  fprintf(postfp,"\n                              sf scalefont setfont }");
  fprintf(postfp,"\n        { sf scalefont setfont } ifelse");
  center=72*8.5/2.-strlen(title)/2.*20.*7./16.;
  if(center<1)
     center=1;
  fprintf(postfp,"\n %f 720 moveto",center);
  /* center the title */
  clean_title_ps(title);

  fprintf(postfp,"\n (%s) show",title);
  fprintf(postfp,"\n /Helvetica findfont");
  fprintf(postfp,"\n 10 scalefont");
  fprintf(postfp,"\n setfont");
  fixcolor_plot_postscript_optimal();
  if(global_contains_prob)
    {if(global_opt_prob_flag)
      {strcpy(str,"Lower Triangle:  Optimal Structure");
      }
     else
       {strcpy(str,"Lower Traingle: Highest Color Range");
       }
    }
  else
    { strcpy(str,"Lower Triangle: Optimal Energy"); 
    }

  fprintf(postfp,"\n 30 90 moveto");
  fprintf(postfp,"\n (%s) show",str);  
  fixcolor_plot_postscript(COLOR_TEXT); 
  if(global_contains_prob==TRUE)
     {fprintf(postfp,"\n /Helvetica findfont");
      fprintf(postfp,"\n 10 scalefont");
      fprintf(postfp,"\n setfont");
      fprintf(postfp,"\n 220 690 moveto");
      strcpy(str,num_string_float(-1.*((float)global_energy_cutoff)/10000.));
      if(global_mi_flag)
        strcat(str," <= Mutual Information <= ");
      else
        strcat(str,"    <= Probability <= ");
      strcat(str,num_string_float(-1.*((float)global_optimal_energy)/10000.));
     }
  else
     {fprintf(postfp,"\n/sf 10 def");
      fprintf(postfp,"\n 130 690 moveto");
      fprintf(postfp,"\ndelta");
      fprintf(postfp,"\n /Helvetica findfont");
      fprintf(postfp,"\n 10 scalefont");
      fprintf(postfp,"\n setfont");
      strcpy(str,"G in Plot File = ");
      strcat(str,num_string_fancy_int(global_energy_cutoff-global_optimal_energy));
      strcat(str," kcal/mole");
     }
  if(global_chain_len>1)
    {strcat(str,"                 Filter: ");
     strcat(str,num_string_int(global_chain_len));
   }
  fprintf(postfp,"\n (%s) show",str); 
/* display colors */
  if(!global_contains_prob)
     {fixcolor_plot_postscript(COLOR_OPTIMAL);
      fprintf(postfp,"\n 340 90  moveto ");
      strcpy(str," Optimal Energy  =  ");
      strcat(str,num_string_fancy_int(global_optimal_energy));
      strcat(str,"  kcal/mole");
      fprintf(postfp,"\n ( %s ) show",str);
      x_float=(float)global_optimal_energy;
      for(color=2,offset=10;color<=global_number_of_colors;color++,offset=offset+10)
         {fixcolor_plot_postscript(color);
          strcpy(str,num_string_fancy_float(
            x_float+(color-2)*global_color_increment));
          strcat(str," < Energy <=");
          strcat(str,
                 num_string_fancy_float(
               x_float+(color-1)*global_color_increment));
          strcat(str,"  kcal/mole");
          fprintf(postfp,"\n 340 %d  moveto",90-offset);
          fprintf(postfp,"\n ( %s ) show",str);
        }
     }
     else
       {fixcolor_plot_postscript(COLOR_OPTIMAL);
        strcpy(str,
        num_string_float(global_prob[global_number_of_colors-4][0]));
        if(global_mi_flag)
          strcat(str,"   < M.I.  <=2.0000");
        else
          strcat(str,"   < Prob. <=1.0000");
        fprintf(postfp,"\n 340 90 moveto");
        fprintf(postfp,"\n (%s) show",str);
        for(color=2,offset=12;color<=global_number_of_colors-1;
                color++,offset=offset+10)
	  {fixcolor_plot_postscript(color);
           strcpy(str,
             num_string_float(global_prob[global_number_of_colors-4][color-1]));
 
           if(global_mi_flag)
             strcat(str,"   < M.I.  <=");
           else  
	     strcat(str,"   < Prob. <=");
	   strcat(str,
             num_string_float(global_prob[global_number_of_colors-4][color-2])); 
           fprintf(postfp,"\n 340 %d moveto",90-offset);
           fprintf(postfp,"\n (%s) show",str);
          }
        fixcolor_plot_postscript(color);
        if(global_mi_flag)
            strcpy(str,"0.0000   < M.I.  <=");
        else
            strcpy(str,"0.0000   < Prob. <=");
 
        strcat(str,
       num_string_float(global_prob[global_number_of_colors-4][global_number_of_colors-2]));
        fprintf(postfp,"\n 340 %d moveto",90-offset);
        fprintf(postfp,"\n (%s) show",str);
       }
  fixcolor_plot_postscript(COLOR_TEXT);
  fprintf(postfp,"\n 25 100 translate");
 /* label top ticks */
  fprintf(postfp,"\n\n%s   Draw Tic Marks \n","%%"); 
  fprintf(postfp,"\n 0.05 setlinewidth");
  if(global_post_file==MAIN_POST)
    plot_hor_tics_post(hscale,make_grid_flag);
  else
    plot_hor_tics_zoom_post(hscale,make_grid_flag);
  /* label side ticks */
  if(global_post_file==MAIN_POST)
    plot_ver_tics_post(vscale,make_grid_flag);
  else
    plot_ver_tics_zoom_post(vscale,make_grid_flag);
/* draw box */ 
  /* left and bottom */
 fprintf(postfp,"\n\n .7 setlinewidth");
 fixcolor_plot_postscript(COLOR_TEXT);
  if(global_post_file!=MAIN_POST)
    {fprintf(postfp,"\n 0 540.5 moveto 0 0 lineto stroke");
     fprintf(postfp,"\n 0 0.5 moveto 540.5 0 lineto stroke");  
   /* draw right and top */
     fprintf(postfp,"\n 0 540.5 moveto 540.5 540.5 lineto stroke");
     fprintf(postfp,"\n 540.5 0  moveto 540.5 540.5 lineto stroke");
    }
/* set scale */
  fprintf(postfp,"\n %f %f scale",hscale,vscale);
  fprintf(postfp,"\n\n %f setlinewidth",0.7/((hscale+vscale)/2.0));
  fprintf(postfp,"\n\n%s   Draw diagonal line \n","%%");
  draw_diagonal_line_post(global_post_file); 
  fprintf(postfp,"\n\n%s   Draw DOTS \n","%%"); 
  if(global_post_file==MAIN_POST)
     draw_postscript_dots_main();
  else
     draw_postscript_dots_zoom();
  if((make_label_flag)||(zoom_labels_count>0))
    {fprintf(postfp,"\n\n%s   Draw Labels \n","%%"); 
     fixcolor_plot_postscript(COLOR_LABEL);
     fprintf(postfp,"\n /Helvetica findfont");
     fprintf(postfp,"\n 8 scalefont");
     fprintf(postfp,"\n setfont");
     if(global_post_file==MAIN_POST)
       { fprintf(postfp,"\n%d %d sbox",
              label_i-1,
              label_j+1);
        fprintf(postfp,"\n %d %d moveto",
              label_j+2,
              global_length-label_i+2);
       }
      else
	{if(make_label_flag)
          {draw_single_label(label_i,label_j,hscale,vscale);
          }
         else
          {for(label_count=0;label_count<zoom_labels_count;label_count++)
            {draw_single_label(zoom_labels_row[label_count],
                               zoom_labels_col[label_count],
                hscale,vscale);
            }
          }
        }
  }
  fixcolor_plot_postscript(COLOR_TEXT);
  fprintf(postfp,"\n %f %f scale",1./hscale,1./vscale);
  fprintf(postfp,"\n -25 -100 translate");
  fprintf(postfp,"\n 27 78 moveto");
  fprintf(postfp,"\n /Helvetica findfont");
  fprintf(postfp,"\n 10 scalefont");
  fprintf(postfp,"\n setfont");
  strcpy(str,"Upper Triangle Base Pairs Plotted: ");
  strcat(str,num_string_int(global_postscript_points_plotted));
  fprintf(postfp,"\n ( %s ) show",str); 
  fprintf(postfp,"\n showpage\n");
  fprintf(postfp,"%sEOF\n","%%");
 fclose(postfp);
}


int open_post(char *filename)
{    /* open specified output file */
       /*  post_filename[strlen(filename) - 1] = NULL; */
   if ((postfp = fopen(filename, "w")) == NULL) /* open a file */
         {printf ("\n * Could not open file:  %s", filename);
          return(1);
          }
    return (0);
}   



