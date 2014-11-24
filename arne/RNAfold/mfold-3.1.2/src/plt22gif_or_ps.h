/* plt22gif_or_ps.h */

/* This file is used by color_table.c, plt22gif.c, plt22ps.c */
#include <stdio.h>
#include <math.h>

#define MAXIMUM_BASES 20000
/* $$$$ Set this maximum sequence length carefully */

FILE *ctab; /* for color table */
char ctab_filename[90];
int ann_to_color[MAXIMUM_BASES]; /* ann_to_color[i] holds the color for */
/* for base i */

int global_total_bases_in_file;

float *global_location_base_x; /* for gives x corridinate of the dot or  */
/* letter for base i */
float *global_location_base_y; /* same for y corridinate  */
char *global_char_base; /* This is the character to placed at x,y  */

/* The 3 values are computed and stored, they are drawn last in order to */
/* to be placed on top of all other drawn characters */

const char *color_table_char[101] =
{                                          /* log probabilities: */
    "c0",     /* red */                /* > .999 */
    "c1",
    "c2",
    "c3",
    "c4",
    "c5",
    "c6",
    "c7",                              /* > .99  */
    "c8",     /* yellow */
    "c9",
    "c10",
    "c11",
    "c12",
    "c13",
    "c14",                              /* > .9   */
    "c15",
    "c16",     /* green */
    "c17",
    "c18",
    "c19",
    "c20",     /* cyan */               /* ~ .5   */
    "c21",
    "c22",
    "c23",
    "c24",     /* blue  */
    "c25",                              /* < .1   */
    "c26",
    "c27",
    "c28",
    "c29",
    "c30",
    "c31",     /* dark magenta */
    "c32",                              /* < .01  */
    "c33",
    "c34",
    "c35",
    "c36",
    "c37",
    "c38",
    "c39",     /* black */              /* < .001 */    
    "c40"      /* must repeat black for position NUM_COLORS*1.00 */
    "c50","c51","c52","c53","c54","c55","c56","c57","c58","c59",
    "c60","c61","c62","c63","c64","c65","c66","c67","c68","c69",
    "c70","c71","c72","c73","c74","c75","c76","c77","c78","c79",
    "c80","c81","c82","c83","c84","c85","c86","c87","c88","c89",
    "c90","c91","c92","c93","c94","c95","c96","c97","c98","c99",
    "c100",

};

struct color_table_struct_ps { float red;
                        float green;
                        float blue;
                      };

struct color_table_struct_ps color_table_ps[NUM_COLORS+1];
struct color_table_struct_ps main_color_table_ps[MAIN_COLORS];

struct color_table_struct_gif { int red;
                                int green;
                                int blue;
                              };

struct color_table_struct_gif main_color_table_gif[MAIN_COLORS];

FILE *ann_file; /* used for ann or ss-count file */

/* the extra colors required for the colored bases are created here */


void set_main_gif_colors()
{int red,green,blue,i;
 int remainder;
  for(i=1;i<MAIN_COLORS;i++)
    {red=main_color_table[i]/65536;
     remainder=main_color_table[i]%65536;
     green=remainder/256;
     blue=remainder%256;
     main_color_table_gif[i].red=red;
     main_color_table_gif[i].green=green;
     main_color_table_gif[i].blue=blue;
    }
}

void set_extra_gif_colors(FILE *giffp,int log_flag)
{ int red,green,blue,i;
  int remainder;
  if(log_flag)
     {
     for(i=0;i<=LOG_NUM_COLORS;i++)
      {red=log_color_table[i]/65536;
       remainder=log_color_table[i]%65536;
       green=remainder/256;
       blue=remainder%256;
       fprintf(giffp,"colorallocate im %s %d %d %d\n",color_table_char[i],
                 red,green,blue);
      }
    }
  else
    {
     for(i=0;i<=NUM_COLORS;i++)
      {red=color_table[i]/65536;
       remainder=color_table[i]%65536;
       green=remainder/256;
       blue=remainder%256;
       fprintf(giffp,"colorallocate im %s %d %d %d\n",color_table_char[i],
                 red,green,blue);
      }
    }
}

void set_main_ps_colors(void)
{int i,red,green,blue,remainder;
 for(i=1;i<MAIN_COLORS;i++)
    {red=main_color_table[i]/65536;
     remainder=main_color_table[i]%65536;
     green=remainder/256;
     blue=remainder%256;
     main_color_table_ps[i].red=(float)red/255.;
     main_color_table_ps[i].green=(float)green/255.;
     main_color_table_ps[i].blue=(float)blue/255.;
    }
}



void set_extra_ps_colors(int prob_flag)
{ int red,green,blue,i,remainder;
  if(prob_flag)
    {for(i=0;i<=LOG_NUM_COLORS;i++)
      {red=log_color_table[i]/65536;
       remainder=log_color_table[i]%65536;
       green=remainder/256;
       blue=remainder%256;
        color_table_ps[i].red=(float)red/256.;
        color_table_ps[i].green=(float)green/256.;
        color_table_ps[i].blue=(float)blue/256.;
      }
    }
  else
   {
    for(i=0;i<=NUM_COLORS;i++)
    {red=color_table[i]/65536;
     remainder=color_table[i]%65536;
     green=remainder/256;
     blue=remainder%256;
     color_table_ps[i].red=(float)red/255.;
     color_table_ps[i].green=(float)green/255.;
     color_table_ps[i].blue=(float)blue/255.;
    }
   } 
 set_main_ps_colors();
} 

void ann_fix_name(char *filename)
{ int i;  /* this function strips off the _structure from the name */
  int len;
  len=strlen(filename);
  /* printf("\n len is %d",len); */
  i=len-1;
  for(i=len-1;((i>1)&&isdigit(filename[i]));i--)
    {/* printf("\n i is %d and character is %c",i,filename[i]); */
    }
  if((i>1)&&(filename[i]=='_'))
     {filename[i]='\0';
      return;
     }
}

void open_ann(char *filename,int global_ann,int global_prob_ann,
      int no_name_change,char *specific_ann_file)
{/* if global_ann is true, strip off _structure if it is there */
 /* Add  .ann when global_ann is true, else .ss-count */
 /* if global_prob_ann is true, leave on Underline structure, */
 /* Then try with it off, add .ss-count */
  /* if no_name_change is FALSE, add .ann or .ss-count only */

  /*  */
  /*  */
  /*  */
  /* Take care of name change option first */
 if(no_name_change!=TRUE)
   {strcpy(filename,specific_ann_file);
    if(global_ann==TRUE)
      strcat(filename,".ann");
    else
      strcat(filename,".ss-count");  
    printf("\n Trying to open  %s",filename);   
    if((ann_file=fopen(filename,"r"))==NULL)
     {printf("\n Could not open file: %s \n",filename);
      exit(1);
     }
   }
 else /* No name change */
   {/* Use default names */
    if(global_prob_ann!=TRUE) /* This is not a probability annotation */
      {/* strip off structure */
         ann_fix_name(filename);
         if(global_ann==TRUE)
           strcat(filename,".ann");
         else
           strcat(filename,".ss-count");  
        printf("\n trying to open %s",filename);
        if((ann_file=fopen(filename,"r"))==NULL)
          {printf("\n Could not open file: %s \n",filename);
           exit(1);
          }
      }
    else
      { /* global_prob_ann is true , try it with structure on first */
         strcpy(specific_ann_file,filename); /* make a copy */
         strcat(filename,".ann");
         printf("\n Trying to open %s",filename);
         if((ann_file=fopen(filename,"r"))!=NULL)
           return; /* File was opened successfully */
         printf("\n Could not open file: %s \n",filename);
         printf("\n I will try something else.");
	 /* strip off structure and try it */
         ann_fix_name(specific_ann_file);
         strcat(specific_ann_file,".ann");
         printf("\n trying to open  %s",specific_ann_file);
         if((ann_file=fopen(specific_ann_file,"r"))==NULL)
           {printf("\n Could not open file: %s \n",specific_ann_file);
            printf("\n I give up\n");
            printf("\n Error !!!!!!!!!!!!!");
            exit(1);
           }
      }
   }
}


void finish_col_table_ss_int_gif(int maximum_pairings,char *ctab_filename,
       FILE *ctab,int save_tgd_flag,char *tgd_filename,int png_mode)
{char string[100],action[200];
 int width,height,error,color;
 int i; 
 int i2;
 int center1,center2;
 int stopping_point;
 float percent_single;
 if(maximum_pairings<20)
   {stopping_point=maximum_pairings;
    width=350;
    center1=350/2;
    height=120+30*(maximum_pairings+1);
   }
 else
   {stopping_point=maximum_pairings/2;
    width=700;
    center1=350/2;
    center2=350+350/2;
    height=120+(30*stopping_point+1);
   }
 fprintf(ctab,"create im %d %d\n",width,height);
 /* declare colors */
 set_main_gif_colors();
 fprintf(ctab,"colorallocate im wh %d %d %d\n",
                 main_color_table_gif[-1*COLOR_BACKGROUND].red,
	         main_color_table_gif[-1*COLOR_BACKGROUND].green,
                 main_color_table_gif[-1*COLOR_BACKGROUND].blue); /* white */
 fprintf(ctab,"colorallocate im bk %d %d %d\n",
                 main_color_table_gif[-1*COLOR_TEXT].red,
	         main_color_table_gif[-1*COLOR_TEXT].green,
                 main_color_table_gif[-1*COLOR_TEXT].blue);  /* black */
 set_extra_gif_colors(ctab,FALSE);
 fprintf(ctab,"stringcenter im gdFont12x24 %d %d \"%s\"  bk\n",
               width/2,30,"Color Table for SS-Count");
 fprintf(ctab,"stringcenter im gdFontGiant %d %d \" %s \" bk\n",
                 center1+12,65,"SS-Count   %   Hex Color");
 for(i=0;i<=stopping_point;i++)
      {if(maximum_pairings==0)
          percent_single=500.0;
       else
          percent_single=100.0*(float)(i)/(float)(maximum_pairings);
       fprintf(ctab,"stringcenter im gdFontGiant %d %d \" %d \" bk\n",center1-60,95+30*i,i); 
       fprintf(ctab,"stringcenter im gdFontGiant %d %d \" %.1f\" bk\n",center1,95+30*i,percent_single);
       color=(NUM_COLORS*
          (float)(maximum_pairings-i)/(float)maximum_pairings+.5);
       fprintf(ctab,"filledrectangle im %d %d %d %d %s\n",50,83+30*i,
                    95,103+30*i,color_table_char[color]);
       fprintf(ctab,"rectangle im %d %d %d %d %s\n",95,83+30*i,
                300,103+30*i,color_table_char[color]);
        fprintf(ctab,"stringcenter im gdFontGiant %d %d \" %.6x \" %s\n",center1+80,95+30*i,color_table[color],color_table_char[color]); 
	
      }
 if(maximum_pairings!=stopping_point)
   {fprintf(ctab,"stringcenter im gdFontGiant %d %d \" %s \" bk\n",
                 center2,65,"SS-Count  %  Hex Color");
     for(i=(stopping_point+1);i<=maximum_pairings;i++)
      {i2=i-stopping_point-1;
       if(maximum_pairings==0)
          percent_single=500.0;
       else
          percent_single=100.0*(float)(i)/(float)(maximum_pairings);
       fprintf(ctab,"stringcenter im gdFontGiant %d %d \" %d \" bk\n",center2-60,95+30*i2,i); 
       fprintf(ctab,"stringcenter im gdFontGiant %d %d \" %.1f \" bk\n",center2-10,95+30*i2,percent_single);
       color=(NUM_COLORS*
          (float)(maximum_pairings-i)/(float)maximum_pairings+.5);
       fprintf(ctab,"filledrectangle im %d %d %d %d %s\n",400,83+30*i2,
                    445,103+30*i2,color_table_char[color]);
       fprintf(ctab,"rectangle im %d %d %d %d %s\n",445,83+30*i2,
                    650,103+30*i2,color_table_char[color]);
       fprintf(ctab,"stringcenter im gdFontGiant %d %d \" %.6x \" %s\n",center2+80,95+30*i2,color_table[color],color_table_char[color]); 
	
      }
   }
 fprintf(ctab,"interlace im 1\n"); /* turn on interlace */  
 strcpy(string,"gif im ");
 strcat(string,ctab_filename); 
 strcat(string,"\n");
 fprintf(ctab,string);
 fclose(ctab);
 if(png_mode)
    strcpy(action,"tgd_png ");
 else
   strcpy(action,"tgd ");
 strcat(action,tgd_filename);
  /*  printf("\n action is %s\n",action);*/
 error=system(action);
 /* uncomment line below to save tgd file */
 /* save_tgd_flag=TRUE; */
  if(error!=0)
    {printf("\n Error converting the color table file to GIF format");
     printf("\n error is %d",error);
     printf("\n The program executes tgd to convert file to GIF format.");
     printf("\n Either it was not installed or it was not installed ");
     printf("\n correctly. ");
     printf("\n \n See plt22gif.doc for more information");
     printf("\n ");
     printf("\n %s should be deleted",tgd_filename);
    }
  else
    {if(save_tgd_flag==TRUE)
      {printf("\n The file %s was preserved.",tgd_filename);
       printf("\n Please delete it when finished!!!!\n");
      }
     else
       {     
        strcpy(action,"rm ");
        strcat(action,tgd_filename);
        error=system(action);
        if(error!=0)
          printf("\n error deleting %s",tgd_filename);  
        }
    }
 
 return;
}     

void  finish_col_table_ann_html(int maximum_pairings,char *ctab_filename,
       FILE *ctab)
{
 int i,color; 
 int start_color[NUM_COLORS+1];/* start_color[i] is the p-num at which color i starts */
 int end_color[NUM_COLORS+1];/* end_color[i] is the p-num at which color i ends */
 int i2;
 int last_color_used;
 int stopping_point;
 float percent_start,percent_end;
 /* set up table of where each color starts and ends */
 last_color_used=0;
 for(i=0;i<=NUM_COLORS;i++)
   {start_color[i]=0;
    end_color[i]=0;
   }
 last_color_used=0;
 color=0;
 start_color[0]=0;
 if(maximum_pairings>NUM_COLORS)  
  {for(i=2;i<=maximum_pairings;i++)
   {color=(int)(NUM_COLORS*((float)i)/(float)(maximum_pairings));
   /* printf("\n f i is %d, color is %d",i,color); */
    if(color>last_color_used)
      {start_color[color]=i;
       end_color[last_color_used]=i-1;
       last_color_used=color;
      }
   }
  end_color[color]=maximum_pairings;
  end_color[NUM_COLORS-1]=maximum_pairings;
  }
 fprintf(ctab,"<HTML><HEAD><TITLE> Color Table: p-num</TITLE>\n");
 fprintf(ctab,"</HTML>\n");
 fprintf(ctab,"<BODY BGCOLOR=\"60ad97\" ALINK=\"00FFFF\" TEXT=\"#000000\" LINK=\"#905000\" VLINK=\"##905000\">\n");
 fprintf(ctab,"<H1 ALIGN=\"CENTER\">Color Table: p-num</H1>\n");
 fprintf(ctab,"<TABLE ALIGN=\"CENTER\" BORDER=\"8\" WIDTH=\"100%%\" ");
 fprintf(ctab,"CELLSPACING=\"2\" BGCOLOR=\"WHITE\" >\n");
 if(maximum_pairings>NUM_COLORS)
  {stopping_point=NUM_COLORS/2-1; 
   fprintf(ctab," <TR>\n");
   fprintf(ctab,"   <TH WIDTH=\"10%%\">Color</TH>\n");
   fprintf(ctab,"   <TH WIDTH=\"10%%\">p-num</TH>\n");
   fprintf(ctab,"   <TH WIDTH=\"3%%\">%%</TH>\n"); 
   fprintf(ctab,"   <TH WIDTH=\"12%%\">Hex Color</TH>\n");
   fprintf(ctab,"   <TH WIDTH=\"12%%\"><FONT COLOR=\"FFFFFF\">. </TH>\n");
   fprintf(ctab,"   <TH WIDTH=\"10%%\">Color</TH>\n"); 
   fprintf(ctab,"   <TH WIDTH=\"10%%\">p-num</TH>\n");
   fprintf(ctab,"   <TH WIDTH=\"3%%\">%%</TH>\n"); 
   fprintf(ctab,"   <TH WIDTH=\"12%%\">Hex Color</TH>\n");
   fprintf(ctab," </TR>");
   for(i=0,i2=stopping_point+1;i<=stopping_point;i++,i2++)
      {fprintf(ctab,"\n <TR>\n");
       fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\" BGCOLOR=\"#%.6x\"><FONT COLOR=\"#%.6x\">^ </FONT> </TD>\n",color_table[i],color_table[i]);
       fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\">");
       if(start_color[i]==end_color[i])
         fprintf(ctab,"%d</TD>\n",start_color[i]);
       else
         fprintf(ctab,"%d-%d</TD>\n",start_color[i],end_color[i]);
       fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\">");
       percent_start=(100.*(float)start_color[i])/((float)maximum_pairings);
       if(start_color[i]==end_color[i])
         {fprintf(ctab,"%.1f </TD>\n",percent_start);
         } 
       else
          {percent_end=(100.*(float)end_color[i])/((float)maximum_pairings);
           fprintf(ctab,"%.1f-%.1f </TD>\n",percent_start,percent_end);
         }
       fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\">");
       fprintf(ctab,"<FONT COLOR=\"#%.6x\">%.6x </FONT></TD>",color_table[i],color_table[i]);
       fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"15%%\"><FONT COLOR=\"FFFFFF\"> . </TD>");
       if(i2<=last_color_used)
	 {
          fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\" BGCOLOR=\"#%.6x\"><FONT COLOR=\"#%.6x\">^ </FONT> </TD>\n",color_table[i2],color_table[i2]);
          fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\">");
          if(start_color[i2]==end_color[i2])
            fprintf(ctab,"%d</TD>\n",start_color[i2]);
          else
            fprintf(ctab,"%d-%d</TD>\n",start_color[i2],end_color[i2]);
          fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\">");
          percent_start=(100. * (float)start_color[i2])/((float)maximum_pairings);
          if(start_color[i2]==end_color[i2])
            fprintf(ctab,"%.1f </TD>\n",percent_start);
          else
            { percent_end=(100. * (float)end_color[i2])/((float)maximum_pairings);
              fprintf(ctab,"%.1f-%.1f </TD>\n",percent_start,percent_end);
            }
          fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\">");
          fprintf(ctab,"<FONT COLOR=\"#%.6x\">%.6x </FONT></TD>",color_table[i2],color_table[i2]);
         }
       else
         { fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\"><FONT COLOR=\"FFFFFF\"> . </TD>");
           fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\"><FONT COLOR=\"FFFFFF\"> . </TD>");
           fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\"><FONT COLOR=\"FFFFFF\"> . </TD>");
           fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\"><FONT COLOR=\"FFFFFF\"> . </TD>");
         }
      fprintf(ctab," </TR>\n");
      }
  }
 else
 { fprintf(ctab," <TR>\n");
   fprintf(ctab,"   <TH WIDTH=\"10%%\">Color</TH>\n");
   fprintf(ctab,"   <TH WIDTH=\"10%%\">p-num</TH>\n");
   fprintf(ctab,"   <TH WIDTH=\"3%%\">%%</TH>\n"); 
   fprintf(ctab,"   <TH WIDTH=\"12%%\">Hex Color</TH>\n");
   fprintf(ctab," </TR>");
   for(i=0;i<=maximum_pairings;i++)
      {if(i>1)
        {color=(int)(NUM_COLORS*((float)(i)/(float)(maximum_pairings)));
        /* printf("\n f i is %d color is %d",i,color); */
        }
       else
        color=1;
        fprintf(ctab,"\n <TR>\n");
       fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\" BGCOLOR=\"#%.6x\"><FONT COLOR=\"#%.6x\">^ </FONT> </TD>\n",color_table[color],color_table[color]);
       fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\">");
       fprintf(ctab,"%d</TD>\n",i);
       fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\">");
       percent_start=(100.*(float)i)/((float)maximum_pairings);
       fprintf(ctab,"%.1f </TD>\n",percent_start);
       fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\">");
       fprintf(ctab,"<FONT COLOR=\"#%.6x\">%.6x </FONT></TD>",color_table[color],color_table[color]);
       fprintf(ctab," </TR>\n");
      }
  }
 fprintf(ctab,"</TABLE>\n");
 fprintf(ctab,"<H3>\n");
 fprintf(ctab,"The color of the <i>i<sup>th</sup></i>base depends on the total number of dots ");
 fprintf(ctab,"in the <i>i<sup>th</sup></i> row and <i>i<sup>th</sup></i> column of the unfiltered <i> energy dot plot </i>(This is the <i>p-num</i>).</H3>");
 fprintf(ctab,"\n<H3>The pairing at the <i>i<sup>th</sup></i> base is more poorly determined if");
 fprintf(ctab," this value is large.\n </H3>");
 fprintf(ctab,"\n<H3>The Maximum <i>p-num</i> is %d </H3>",maximum_pairings);
 fprintf(ctab,"\n<H3> See <A HREF=\"%s\"> structure annotation:</A></H3>",http_color_doc);
 fprintf(ctab,"</BODY>\n");
 fclose(ctab);
 return;
}

void  finish_col_table_ss_int_html(int maximum_pairings,char *ctab_filename,
       FILE *ctab)
{
 int i,color,color2; 
 int i2;
 int stopping_point;
 float percent_single,percent_single2;
 if(maximum_pairings<20)
   {stopping_point=maximum_pairings;
   }
 else
   {stopping_point=maximum_pairings/2;
   }
 fprintf(ctab,"<HTML><HEAD><TITLE> SS-Count Color Table</TITLE>\n");
 fprintf(ctab,"</HTML>\n");
 fprintf(ctab,"<BODY BGCOLOR=\"#60ad97\" ALINK=\"00FFFF\" TEXT=\"#000000\" LINK=\"#905000\" VLINK=\"#905000\">\n");
 fprintf(ctab,"<H1 ALIGN=\"CENTER\">Color Table: ss-count</H1>\n");
 fprintf(ctab,"<TABLE ALIGN=\"CENTER\" BORDER=\"8\" WIDTH=\"100%%\" ");
 fprintf(ctab,"CELLSPACING=\"2\" BGCOLOR=\"WHITE\">\n");
 fprintf(ctab," <TR>\n");
 if(maximum_pairings<20)
   {
    fprintf(ctab,"    <TH WIDTH=\"10%%\">Color</TH>\n");
    fprintf(ctab,"    <TH WIDTH=\"10%%\">ss-count</TH>\n");
    fprintf(ctab,"    <TH WIDTH=\"5%%\" >%%</TH>\n"); 
    fprintf(ctab,"    <TH WIDTH=\"10%%\">Hex Color</TH>\n");
   }
 else
    {fprintf(ctab,"   <TH WIDTH=\"10%%\">Color</TH>\n");
     fprintf(ctab,"   <TH WIDTH=\"10%%\">ss-count</TH>\n");
     fprintf(ctab,"   <TH WIDTH=\"3%%\" >%%</TH>\n"); 
     fprintf(ctab,"   <TH WIDTH=\"12%%\">Hex Color</TH>\n");
     fprintf(ctab,"   <TH WIDTH=\"12%%\"><FONT COLOR=\"FFFFFF\">.</TH>\n");
     fprintf(ctab,"   <TH WIDTH=\"10%%\">Color</TH>\n"); 
     fprintf(ctab,"   <TH WIDTH=\"10%%\">ss-count</TH>\n");
     fprintf(ctab,"   <TH WIDTH=\"3%%\">%%</TH>\n"); 
     fprintf(ctab,"   <TH WIDTH=\"12%%\">Hex Color</TH>\n");
    }
 fprintf(ctab," </TR>");
 if(maximum_pairings<20)
   {
    for(i=0;i<=stopping_point;i++)
      {color=(NUM_COLORS*
          (float)(maximum_pairings-i)/(float)maximum_pairings+.5);
       if(maximum_pairings==0)
          percent_single=500.0;
       else
          percent_single=100.0*(float)(i)/(float)(maximum_pairings);
       fprintf(ctab,"\n <TR>\n");
       fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\" BGCOLOR=\"#%.6x\"><FONT COLOR=\"#%.6x\">^ </FONT> </TD>\n",color_table[color],color_table[color]);
       fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\">%d</TD>\n",i);
       fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"5%%\">%.1f </TD>\n",percent_single);
       fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\">");
       fprintf(ctab,"<FONT COLOR=\"#%.6x\">%.6x </FONT></TD>",color_table[color],color_table[color]);
       fprintf(ctab," </TR>\n");
      }
   }
 else
   { 
    for(i=0,i2=stopping_point+1;i<=stopping_point;i++,i2++)
      {color=(NUM_COLORS*
          (float)(maximum_pairings-i)/(float)maximum_pairings+.5);
       color2=(NUM_COLORS*
          (float)(maximum_pairings-i2)/(float)maximum_pairings+.5);
       percent_single=100.0*(float)(i)/(float)(maximum_pairings);
       percent_single2=100.0*(float)(i2)/(float)(maximum_pairings);
       fprintf(ctab,"\n <TR>\n");
       fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\" BGCOLOR=\"#%.6x\"><FONT COLOR=\"#%.6x\">^ </FONT> </TD>\n",color_table[color],color_table[color]);
       fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\">%d</TD>\n",i);
       fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\">%.1f </TD>\n",percent_single);
       fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\">");
       fprintf(ctab,"<FONT COLOR=\"#%.6x\">%.6x </FONT></TD>",color_table[color],color_table[color]);
       fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"15%%\"><FONT COLOR=\"FFFFFF\"> . </TD>");
       if(i2<=maximum_pairings)
	 {
          fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\" BGCOLOR=\"#%.6x\"><FONT COLOR=\"#%.6x\">^ </FONT> </TD>\n",color_table[color2],color_table[color2]);
          fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\"> %d</TD>\n",i2);
          fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\"> %.1f </TD>\n",percent_single2);
          fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\">");
          fprintf(ctab,"<FONT COLOR=\"#%.6x\">%.6x </FONT></TD>",color_table[color2],color_table[color2]);
         }
       else
         { fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\"> <FONT COLOR=\"FFFFFF\"> . </TD>");
           fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\"> <FONT COLOR=\"FFFFFF\"> . </TD>");
           fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\"> <FONT COLOR=\"FFFFFF\"> . </TD>");
           fprintf(ctab,"    <TD ALIGN=\"CENTER\" WIDTH=\"10%%\"> <FONT COLOR=\"FFFFFF\"> . </TD>");
         }
      fprintf(ctab," </TR>\n");
      }
   }
 
 fprintf(ctab,"</TABLE>\n");
 fprintf(ctab,"<H3>\n");
 fprintf(ctab,"The color of the <i>i<sup>th</sup></i> base depends on the number of structures in which it is single stranded.</H3>");
 fprintf(ctab,"\n<H3> Total Structures: %d</H3>",maximum_pairings);
 fprintf(ctab,"\n<H3> See <A HREF=\"%s\"> structure annotation:</A></H3>",
  http_color_doc);
 fprintf(ctab,"</BODY>\n");
 fclose(ctab);
 return;
}

/* The above link is used on html color tables generated during
   annotation */
void make_col_table_ss_int(int maximum_pairings,char *filename,
    int save_tgd_flag,char color_table_type,int png_mode)
{
 char tgd_filename[100];
 strcpy(ctab_filename,filename);
 strcat(ctab_filename,".col");
 if(color_table_type=='g')
   {strcpy(tgd_filename,ctab_filename);
    strcat(tgd_filename,".tgd");
    if(png_mode)
      strcat(ctab_filename,".png");
    else
      strcat(ctab_filename,".gif");
   }
 else
   strcat(ctab_filename,".html");
 printf("\n Creating color table %s\n",ctab_filename);
 if(color_table_type=='g')
     {
      if((ctab=fopen(tgd_filename,"w"))==NULL)
          {printf("\n Could not open file: %s \n",tgd_filename);
           return;
          }
      finish_col_table_ss_int_gif(maximum_pairings,ctab_filename,
       ctab,save_tgd_flag,tgd_filename,png_mode);
     }
 else
    {
     if((ctab=fopen(ctab_filename,"r"))!=NULL)
       {printf("\n %s already exists, not created\n",ctab_filename);
        close(ctab);
        return;
        }
     if((ctab=fopen(ctab_filename,"w"))==NULL)
          {printf("\n Could not open file: %s \n",ctab_filename);
           return;
          }
      finish_col_table_ss_int_html(maximum_pairings,ctab_filename,
       ctab);
    }
 return;
}     

void make_col_table_ann(int maximum_pairings,char *filename,
    int save_tgd_flag,char color_table_type)
{
 char tgd_filename[100];
 strcpy(ctab_filename,filename);
 strcat(ctab_filename,".col");
 if(color_table_type=='g')
   {
    strcpy(tgd_filename,ctab_filename);
    strcat(tgd_filename,".tgd");
    strcat(ctab_filename,".gif");
   }
 else
   {strcat(ctab_filename,".html");
    if((ctab=fopen(ctab_filename,"r"))!=NULL)
     {printf("\n %s Already exists, not created\n",ctab_filename);
      close(ctab);
      return;
     }
   }

 if((ctab=fopen(ctab_filename,"w"))==NULL)
          {printf("\n Could not open file: %s \n",ctab_filename);
           return;
          }
 if(color_table_type=='g')
     { printf("\n This color table is not implemented yet");
     /* if((ctab=fopen(tgd_filename,"w"))==NULL)
          {printf("\n Could not open file: %s \n",tgd_filename);
           return;
          }
      finish_col_table_ss_int_gif(maximum_pairings,ctab_filename,
       ctab,save_tgd_flag,tgd_filename); */
     }
 else
    {if((ctab=fopen(ctab_filename,"w"))==NULL)
          {printf("\n Could not open file: %s \n",ctab_filename);
           return;
          }
      printf("\n color annotation file is %s",ctab_filename);
      finish_col_table_ann_html(maximum_pairings,ctab_filename,
       ctab);
    }
 return;
}     

void define_map_ann(char *filename,int global_ann,int global_prob_ann,int 
           no_name_change,char *specific_ann_file,int table_flag,
           char color_table_type,int save_tgd_flag,int png_mode)
    /* This function is called for ann or ss-count */
{int maximum_pairings;
 int base,pairs,i;
 float color;
 char rec[90];
 int upper_limit;
 upper_limit=MAXIMUM_BASES;
 open_ann(filename,global_ann,global_prob_ann,
       no_name_change,specific_ann_file);
 /*   */
 /* Read probability Annotation First */
 /*   */
 if(global_prob_ann==TRUE)
   {maximum_pairings=1;
    global_total_bases_in_file=0;
     while(fgets(rec,90,ann_file)!=NULL)
       {sscanf(rec,"%d%f",&base,&color);
       if((color>1.)||(color<0.))
           {printf("\n probability %f not in 0<=prob<=1 in file ",color);
            printf("\n Error !!!!!!!!!!!!!!\n");
            exit(1);
           }
        if(base >= upper_limit)
         {printf("\n Error Too Many bases -----------------");
	  printf("\n Adjust the constant MAXIMUM_BASES in ");
          printf("\n plt22gif_or_ps.h and recompile");
          printf("\n current limit is %d, base was %d",upper_limit,base);
          exit(1);
         }
        ann_to_color[base]=get_color((double)(color),TRUE,TRUE);
        /* printf("\n color %d is %d",base,ann_to_color[base]); */
        if(base>global_total_bases_in_file)
          global_total_bases_in_file=base;
       }
     if(global_total_bases_in_file==0)
       {printf("\n Error ___________________No data in ann file ");
        exit(1);
       }
   }
 else 
  {if(global_ann==TRUE)
   { /* ann option */
     maximum_pairings=1; /* if all zero , this avoids zero division below */
     global_total_bases_in_file=0;
     while(fgets(rec,90,ann_file)!=NULL)
       {sscanf(rec,"%d%f",&base,&color);
        pairs=(int)(color);
        if((color<1.0)&&(color>0.0))
	   {printf("\n\n !!!!!!!!!!!Error, !!!!");
            printf("\n\n  Attempted to treat a probabilty file ");
            printf("\n     as a p-num file. \n\n");
            exit(1);
           }
        if(base >= upper_limit)
         {printf("\n Error Too Many bases -----------------");
	  printf("\n Adjust the constant MAXIMUM_BASES in ");
          printf("\n plt22gif_or_ps.h and recompile");
          printf("\n current limit is %d, base was %d",upper_limit,base);
          exit(1);
         }
        if(pairs==1)
           pairs=0; /* p-num = 0 or 1 are both equally well defined */
        ann_to_color[base]=pairs;
	/*  printf("\n for base %d pairs is %d",base,pairs); */
        if(pairs>maximum_pairings)
          maximum_pairings=pairs;
        global_total_bases_in_file++;
       }
     if(global_total_bases_in_file==0)
       {printf("\n Error ___________________No data in ann file ");
        exit(1);
       }
     /* new changes */
     for(i=1;i<=global_total_bases_in_file;i++)
       {color=NUM_COLORS*((float)ann_to_color[i]/
           (float)maximum_pairings); /* avoid division by zero */
        /* scaled color based on the possible number of times a base can */
        /* pair */ 
       /* printf("\n pairings is %d  ",ann_to_color[i]); */
        ann_to_color[i]=(int)color;
	/*  printf(" i is %d, color is %d, maximum_pa is %d",
               i,ann_to_color[i],maximum_pairings); */
       }
   if(table_flag)
      {if(color_table_type=='g')
	{ printf("\n The -t g option with a p-num file does not work.");
          printf("\n The code to create the .gif format for the color table\n");
         printf("\n will be created if necessary. ");
	 /*
         make_col_table_ss_int(maximum_pairings,filename,save_tgd_flag,'g'); */
        }
       else
	 {  
          if(color_table_type=='h')
            make_col_table_ann(maximum_pairings,filename,save_tgd_flag,'h');
	  } 
      }  
   }
 else
   {if(fgets(rec,90,ann_file)==NULL) 
     { printf("\n Error _____________________No data in ss-count file");
       exit(1);
     }
    sscanf(rec,"%d",&maximum_pairings);
    global_total_bases_in_file=0;
    while(fgets(rec,90,ann_file)!=NULL)
      {sscanf(rec,"%d%d",&base,&pairs);
       if(base >= upper_limit)
         {printf("\n Error Too Many bases -----------------");
	  printf("\n Adjust the constant MAXIMUM_BASES in ");
          printf("\n plt22gif_or_ps.h and recompile");
          printf("\n current limit is %d, base was %d\n",upper_limit,base);
          exit(1);
         }
       /* printf("\n base is %d pairs is %d ",base,pairs); */
       color=NUM_COLORS*
          (float)(maximum_pairings-pairs)/(float)maximum_pairings;
       ann_to_color[base]=(int)(color+.5);
        /* scaled color based on the possible number of times a base can */
       /* be in a basepair */
       if(base>global_total_bases_in_file)
         global_total_bases_in_file=base;
      }
    if(table_flag)
      {if(color_table_type=='g')
         make_col_table_ss_int(maximum_pairings,filename,save_tgd_flag,
            'g',png_mode);
       else
	 {if(color_table_type=='h')
            make_col_table_ss_int(maximum_pairings,filename,save_tgd_flag,'h',
              FALSE);
	 }
      }      
   }
  }
   printf("\n There are  %d bases.",global_total_bases_in_file);
   global_location_base_y=(float *)malloc(
            (global_total_bases_in_file+1)*sizeof(float));
   global_char_base=(char *)malloc(
              (global_total_bases_in_file+1)*sizeof(char)); 
   global_location_base_x=(float *)malloc(
              (global_total_bases_in_file+1)*sizeof(float));
   if((global_location_base_x==NULL)||(global_location_base_y==NULL)
      ||(global_char_base==NULL))
     {printf("\n -------Error-----------------");
      printf("\n insufficient memory to allocate space for bases ");
      printf("\n from plt22gif,plt22ps with colored bases or dots \n");
      exit(1);
     }
 fclose(ann_file);
} 


