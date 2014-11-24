/* boxplot_input.inc */
/* By Darrin Stewart and Michael Zuker */
/* May 16, 2000 */

/* strip off .plot from end of title */
void fix_title(char *title)
{int i;
   i=strlen(title);
   if(i>5)
     {i=i-5;
      if((title[i]=='.')&&(title[i+1]=='p')&&(title[i+2]=='l')&&
          (title[i+3]=='o')&&(title[i+4]=='t'))
          {title[i]='\0';
           return;
          }
     } 
    return;
}




int count_helices(void)
{ int helices;
  char rec[100];
  helices=0;
  while(fgets(rec,90,fp)!=NULL)
    {helices++;
    }
  if(helices<1)
   {printf("\n Error  reading plot file, no data in file");
     printf("\n Abnormal termination from boxplot\n");
     exit(1);
   }
 printf("\n number of helices is %d",helices-1);
 return helices;
}

void check_parameters(int,char **);
void initialize_data2(void);

void open_file(int argc, char **argv)
{int helix_total; 
 helix_total=0;
 check_parameters(argc,argv);
 printf("\n Trying to open  %s",global_plot_filename);   
 if((fp=fopen(global_plot_filename,"r"))==NULL)
  {printf("\n ---------Could not open plot file: %s.----------------------- ",
         global_plot_filename);
   printf("\n\n -End of session . . . good bye \n\n");
   exit(1);
  } 
 helix_total=count_helices();
 fclose(fp);
 global_helix_array_size=helix_total;
 global_diag=malloc(sizeof(struct helix)*(helix_total+1));
 if(global_diag==NULL)
   {printf("\n Error in boxplot, insufficient memory to ");
    printf("\n execute with %d helices in plot file: %s",helix_total,
                     global_plot_filename);
    exit(1);
   }
 /* printf("\n helix total is %d",helix_total); */
 if((fp=fopen(global_plot_filename,"r"))==NULL)
  {printf("\n ---------Could not open file.----------------------- ");
   printf("\n\n -End of session . . . good bye \n\n");
   exit(1);
  }  
}

int compare_helices(const void *h1,const void *h2)
{ int dif;
  struct helix *helix1,*helix2;
  helix1=(struct helix *)h1;
  helix2=(struct helix *)h2; 
  dif=(*helix1).diagonal-(*helix2).diagonal; /* return <0 if diag1 < diag2 */
  if(dif==0) /* on same diagonal */
     dif=(*helix1).row-(*helix2).row; /* return <0 if row1<row2 or >0 if > */
  return dif;
} 

void sort_helices(void)
{qsort(global_diag,global_helices_in_plot_file,
       sizeof(struct helix),compare_helices);
}
  
void print_helix(int number_of_helices)
{int i;
 for(i=0;i<number_of_helices;i++)
    { printf("\n diag=%d, row=%d, color=%d, column=%d, energy=%d, length=%d",
            global_diag[i].diagonal,
            global_diag[i].row,
            global_diag[i].color,
            global_diag[i].column,
            global_diag[i].energy,
            global_diag[i].length);
    }
}

void fix_counts(void)
{int i;
 for(i=0;i<global_helices_in_plot_file;i++)
   {global_diag_count[global_diag[i].diagonal]++; /* update the count for this diagonal */
   }
}

int find_color(int,int,float,int);

void fix_colors(void) /* set the colors for each line of input */
{ int i;
  for (i=0;i<global_helices_in_plot_file;i++)
    {global_diag[i].color=find_color(global_diag[i].energy,
          global_contains_prob,global_color_increment,global_optimal_energy);
    }
}

void initialize_len(int adjust_col_flag, int adjust_amount)
{/* global_worst_energy and global_optimal_energy are also set */
 int level,length,row,col,energy,diag;
 int i;
 /* diag is the diagonal that each row,col corresponds to */
 float rec_energy;
 global_opt_prob_flag=FALSE;
 global_helices_in_plot_file=0;

 if(fgets(rec,90,fp)==NULL)
     {printf("\n \n Error------   Plot file is completely empty\n");
      exit(1);
     }

 global_length=0;
 /* Read the plot file */
 global_worst_energy=INT_MIN;
 global_optimal_energy=INT_MAX;
 while(fgets(rec,90,fp)!=NULL)
   {sscanf(rec,"%d%d%d%d%f",&level,&length,&row,&col,&rec_energy);
    if(adjust_col_flag)
      {printf("\n (%d,%d) converted to ",row,col);
        col=col-adjust_amount;
       if((col-length+1)<1)
         {printf("\n error, column %d rounded up to %d",col+adjust_amount,
                       length);
          col=length;
         }
       printf(" (%d,%d) ",row,col);
       if(global_length<(row+length-1))
          global_length=row+length-1;
      }
    if(global_length<col)
       global_length=col;
    global_diag[global_helices_in_plot_file].length=length;
    global_diag[global_helices_in_plot_file].row=row;
    global_diag[global_helices_in_plot_file].column=col;
    if(global_contains_prob)
      energy=(int)(-1*(10000.*rec_energy+.50));
    else
      energy=(int)(rec_energy);
    global_diag[global_helices_in_plot_file].energy=energy;
    if((global_contains_prob)&&(energy<(-1*2*10000+1)))
     {global_opt_prob_flag=TRUE;
     }
    else
      {if(energy<global_optimal_energy)
          global_optimal_energy=energy;
      } 
    if(energy>global_worst_energy)
          global_worst_energy=energy;
    diag=row+col-1;
    global_diag[global_helices_in_plot_file].diagonal=diag;
    global_helices_in_plot_file++;
    if(global_helices_in_plot_file==global_helix_array_size)
         {printf("\n Error!!!___________________________________");
          printf("\n There was an error reading the plot file");
          printf("\n expected %d helices",global_helix_array_size);
          printf("\n The lines in file was less than the number of helices  in boxplot.c");
          printf("\n Strange error\n");
          exit(1);
         }
   }
 printf("\n There are %d helices in the plot file.",global_helices_in_plot_file);

 if(global_helices_in_plot_file==0)
   {printf("\n \n Error------   Plot file contains no helices\n");
    exit(1);
   }
 /* glob_diag_count counts the total number of entries in each diagonal */
 /* printf("\n Farthest basepair is  %d\n",global_length); */
 fclose(fp);
 /* Sort the plot file */
 /* printf("\n just before sort"); */
 /* print_helix(global_helices_in_plot_file); */
 sort_helices();
 /* printf("\n just after sort");*/
 global_diag_start=(int *)malloc((2*global_length)*sizeof(int));
 if(global_diag_start==NULL)
    printf("\n Error insufficient memory for global_diag_start");
 global_diag_count=(int *)malloc((2*global_length)*sizeof(int));
 if(global_diag_count==NULL)
    printf("\n Error insufficient memory for global_diag_count"); 
 /* print_helix(global_helices_in_plot_file); */
 /*   printf("\n just before set count to zero"); */
    for(i=0;i<=(2*global_length-1);i++)
       global_diag_count[i]=0;
    /* fix the count of diagonal for each entry in plot file */
    /*   printf("\n worst_energy=%d, optimal_energy=%d, ",
            global_worst_energy,global_optimal_energy); */
  initialize_data2();

  /*  printf("\n color increment is %f",global_color_increment); */ 
  if(global_opt_prob_flag==TRUE)
    {printf("\n   Optimal structure found in data");
    }
  fix_colors();
  /* printf("\n just after set colors"); */
  /*  print_helix(global_helices_in_plot_file); */
  global_diag_start[1]=0;
  fix_counts();
  for(i=2;i<=(2*global_length-1);i++)
   {global_diag_start[i]=global_diag_start[i-1]+global_diag_count[i-1];
   /*   printf("\n for diag=%d count is %d start is %d", 
         i,global_diag_count[i],global_diag_start[i]);*/
   }
 /* global_diag_start counts forward into each array for the start of */
 /* each diagonal */
}

/*_______________________________________________________________________*/
/* Conversion Routines to display integers and floats */
char global_str[30];

char *num_string_int(int x) /* converts integer to string */
{    sprintf(global_str,"%d",x);
     return global_str;
}

char *num_string_float(float x) /* convert float to string */
{     sprintf(global_str,"%.4f",x);
     return global_str;
}

char *num_prob_string(int x) /* divides by 10000 and multiplies by -1 */
{    float x_float;
     x_float=(float)x;
     x_float=x_float/10000.;
     x_float=-1*x_float;
     /* round to thousandth position */
     sprintf(global_str,"%4.4f",x_float);
     return global_str;
} 

char *num_energy_string(int x) /* divides number by 10 and converts to */
			     /* string */
{     float x_float;
      x_float=(float)x;
      x_float=x_float/10.;
      /* round to tenth position */
        sprintf(global_str,"%7.1f",x_float);
      return global_str;
}
char *num_prob_string_float(float x) /* divides by 1000 and multiplies by -1 */
{    x=x/-10000.;
     /* round to thousandth position */
     sprintf(global_str,"%4.4f",x);
     return global_str;
} 

char *num_energy_string_float(float x) /* divides number by 10 and converts to */
			     /* string */
{     x=x/10.;
      /* round to tenth position */
      sprintf(global_str,"%7.1f",x);
      return global_str;
}

char *num_string_fancy_int(int x) /* converts integer to string , handles prob */
			/* or energies*/
{    float float_x;
     float_x=(float)x;
     if(global_contains_prob==TRUE)
         return strcpy(global_str,num_prob_string_float(float_x));
      else
         return strcpy(global_str,num_energy_string_float(float_x));
}
char *num_string_fancy_float(float x) /* converts float to string */
{     if(global_contains_prob==TRUE)
        return strcpy(global_str,num_prob_string_float(x));
      else
        return strcpy(global_str,num_energy_string_float(x));
}

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


/* other functions _______________________________*/

void finish_setting_energy_cutoff(float temp_filter)
    {  global_energy_cutoff=(int)(temp_filter*10);
       if(global_energy_cutoff>=0)
           global_energy_cutoff=global_optimal_energy+global_energy_cutoff;
       if(global_energy_cutoff>global_worst_energy)
          global_energy_cutoff=global_worst_energy;
       if(global_energy_cutoff<global_optimal_energy)
          global_energy_cutoff=global_optimal_energy;
       global_color_increment=(float)(global_energy_cutoff-
                  global_optimal_energy)/(float)(global_number_of_colors-1);

}

/* given a row and column, return the energy of that position */
int energy_read(int row,int col)
{int diag,diag_start,diag_end,i,temp_row;
 diag=col+row-1;
 diag_start=global_diag_start[diag];
 diag_end=global_diag_start[diag+1];
 for(i=diag_start;i<diag_end;i++)
     {temp_row=global_diag[i].row;
      if(temp_row>row) 
        return 0;   /*  entries in diagonal are sorted */
      if((!global_opt_prob_flag)||(global_diag[i].energy>=-20000))
	{
      if((temp_row<=row)&&(global_diag[i].row+global_diag[i].length-1>=row))  
        return global_diag[i].energy;
        }
     }
 return 0;
} 
int energy_read_opt_prob(int row,int col)
{int diag,diag_start,diag_end,i,temp_row;
 diag=col+row-1;
 diag_start=global_diag_start[diag];
 diag_end=global_diag_start[diag+1];
 for(i=diag_start;i<diag_end;i++)
     {temp_row=global_diag[i].row;
      if(temp_row>row) 
        return 0;   /*  entries in diagonal are sorted */
      if(global_diag[i].energy<-20000)
	{
      if((temp_row<=row)&&(global_diag[i].row+global_diag[i].length-1>=row))  
        return global_diag[i].energy;
        }
     }
 return 0;
} 


FILE *engfp;

int open_energy(char *filename)
{    /* open specified output file */
       /*  post_filename[strlen(filename) - 1] = NULL; */
   if ((engfp = fopen(filename, "w")) == NULL) /* open a file */
         {printf ("\n * Could not open file:  %s", filename);
          return(1);
          }
    return (0);
}


void make_energy_file(char *filename,int label_i,int label_j)
{/* make energy file */
 int error;
 float energy;
 strcat(filename,".gifeng");
 printf("\n gifeng is %s at location (%d,%d) \n",filename,label_i,label_j);
 error=open_energy(filename);
 if(error)
   {printf("Error creating .gifeng file name was %s",filename);
    return;
   }
 if((label_i<1)||(label_i>global_length)||(label_j<1)||(label_j>global_length))
     energy=9999.9;
 else
     energy=((float)energy_read(label_i,label_j)+.0001)/10.;
  fprintf(engfp,"%.1f\n",energy);
  fclose(engfp);
}

