/* ct_boxplot_read_ct.h */
/* by Darrin Stewart and Michael Zuker */
/* Feb 4, 1999 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

char sequence[MAXIMUM_SIZE];

void print_helix(int number_of_helices,struct helix *diag)
{int i;
 for(i=0;i<number_of_helices;i++)
    { printf("\n diag=%d, row=%d, color=%d, column=%d, length=%d",
            diag[i].diagonal,
            diag[i].row,
            diag[i].color,
            diag[i].column,
            diag[i].length);
    }
}

void fix_name(char *sequence_name)
{ /* remove new line and other annoying characters */
  int i;
  int len;
  int non_space;
  for(i=0;i<strlen(sequence_name);i++)
	{if(sequence_name[i]=='\n')
            sequence_name[i]=' ';
        }
  len=strlen(sequence_name);
  i=len-1;
  /* remove extra spaces at end of name */
  while( (i>4)&&(sequence_name[i]==' '))
      {
       sequence_name[i]='\0';
       i--;
      }
  /* remove leading spaces at start of name */
  /* find first non-space */
 len=strlen(sequence_name);
 len--;
 i=0;
 non_space=0;
 while((sequence_name[i]==' ')&&(i<len))
  {non_space++;
   i++;
  }
 if(non_space>0)
   {for(i=non_space;i<=strlen(sequence_name);i++)
     sequence_name[i-non_space]=sequence_name[i];
   }
  
  
}

int compare_helices(const void *h1,const void *h2)
{ int dif;
  struct helix *helix1,*helix2;
  helix1=(struct helix *)h1;
  helix2=(struct helix *)h2; 
  dif=(*helix1).diagonal-(*helix2).diagonal; /* return <0 if diag1 < diag2 */
  if(dif==0)                            /* on same diagonal */
     dif=(*helix1).row-(*helix2).row; /* return <0 if row1<row2 or >0 if > */
  return dif;
}   

void sort_helices(struct helix *diag,int helices_in_plot_file)
{qsort(diag,helices_in_plot_file,
       sizeof(struct helix),compare_helices);
}
void append_structure_name(char *data_name,char *file_name)
{/* place name from file at start of description,data_name */
 char temp_name[120];
 int length;
 char number[7];
 strcpy(temp_name,file_name);
 fix_title(temp_name); /* remove .ct from name */
 length=strlen(temp_name);
 /* look for _#, _## or _### at end of name */
 if(isdigit(temp_name[length-1]))
    {number[0]=' ';
     number[1]='[';
     if(temp_name[length-2]=='_')
      {/* one digit name */
        number[2]=temp_name[length-1];
        number[3]=']';
        number[4]='\0';
        strcat(data_name,number);
      }
     else
      { if(isdigit(temp_name[length-2]))
	 {if(temp_name[length-3]=='_')
	   {/* two digit name */
             number[2]=temp_name[length-2];
             number[3]=temp_name[length-1];
             number[4]=']';
             number[5]='\0';
             strcat(data_name,number);
           }
         else
           {if((isdigit(temp_name[length-3]))&&(temp_name[length-4]=='_'))
	     {/* three digit name */
              number[2]=temp_name[length-3];
              number[3]=temp_name[length-2];
              number[4]=temp_name[length-1];
              number[5]=']';
              number[6]='\0';
              strcat(data_name,number);
             }
           }
	 }  
      }

    }
}
void read_ct_file(FILE *fp,int file,int *len,struct helix *diag,
             char file_data_ct[MAXIMUM_CT_FILES][80],
           int *helices_in_plot_file,int *first_sequence_length,
           char *global_sequence_name,int sequence_name_set,char *filename)
{char record[160];
 int file_length;
 int row,previous,next,column,junk;
 int helix_len;
 int helices_in_this_file; 
 char base;
 char energy_string[140];
 int first_char;
 int helix_started;
 int j,i;
 int sequence_length;
 int previous_column;
 int name_start;
 int energy_start;
 int record_length;
 char c1,c2;
 char junk_length[150];
 int last_char;
 char sequence_name[150];
 sequence_length=0;
 helices_in_this_file=*helices_in_plot_file;
 /* diag is the diagonal that each row,col corresponds to */
 if(fgets(record,90,fp)==NULL)
     {printf("\n \n Error------   CT file is completely empty\n");
      exit(1);
     }
 else
    {/* record should look like  741 ENERGY =  -238.3    mok-lacZ + gtgtaa */
      /* 741 is lenght of sequence, mok-lacZ + gtgtaa is the name */
      /* and should be the same in all files on a give run */
      /* The energy varies and is stored as file_data_ct */
     sscanf(record,"%d",&file_length); /* assign length */
     sscanf(record,"%s",junk_length); /* assign length string */
     first_char=0;
     name_start=0;
     energy_start=0;
     record_length=strlen(record);
     while((isspace(record[first_char]))&&(first_char<record_length))
	   first_char++;/* Advance through spaces of record */
     while(isdigit(record[first_char])) /*advance through numbers of record */
           first_char++;
     while((isspace(record[first_char]))&&(first_char<record_length))
           first_char++; /* advance through spaces */
     energy_start=first_char;
     while((!isspace(record[first_char]))&&(first_char<record_length))
           first_char++; /* advance through energy */
     while((isspace(record[first_char]))&&(first_char<record_length))
           first_char++; /* advance through spaces */
     while((!isspace(record[first_char]))&&(first_char<record_length))
           first_char++; /* advance through = */
     while((isspace(record[first_char]))&&(first_char<record_length))
           first_char++; /* advance through spaces */
     while((!isspace(record[first_char]))&&(first_char<record_length))
           first_char++; /* advance through = number*/
     name_start=first_char;
     /* this should copy energy into file_data_ct[file]*/
     sprintf(file_data_ct[file],"%d",file+1);
     /* if the filename ends with _# append that number to */
     /* file_data_ct */
     append_structure_name(file_data_ct[file],filename);
     if(name_start>30)
       {/* failed to find energy */
        name_start=6;
        strcpy(energy_string,"Energy is Undef.");
       }
     else
       {
         for(i=energy_start,j=0;i<name_start;i++,j++)
           energy_string[j]=record[i];
         energy_string[j]='\0';
       }
     last_char=strlen(energy_string)-1;
     if(energy_string[last_char]=='\n')
         energy_string[last_char]='\0';
     strcat(file_data_ct[file]," ");
     strcat(file_data_ct[file],energy_string);
     if((file==0)&&(!sequence_name_set)) /* copy the name in */
       {for(i=name_start,j=0;(i<record_length)&&(j<79);i++,j++)
          sequence_name[j]=record[i];
        sequence_name[j]='\0';
        fix_name(sequence_name);
        strcpy(global_sequence_name,sequence_name);
       }
     }
 /* Read the ct file */
 helix_started=FALSE;
 helix_len=0;
 previous_column=-10;
 sequence_length=0; 
 while(fgets(record,90,fp)!=NULL)
   {sscanf(record,"%d %c %d %d %d %d",
            &row,&base,&previous,&next,&column,&junk);
   /* printf("\n reading in %d %s %d %d %d %d previous col is %d",row,base,
        previous,next,column,previous_column); */
    if(row!=(sequence_length+1))
      {printf("\n\nError in ct file, Format is wrong");
       printf("\n Record %d of the file contains \n %s",
             sequence_length+1,record);
       printf("\n but it should start with %d\n",sequence_length+1);
       exit(1);
      }
    if(file==0)
      sequence[sequence_length]=base;
    else
      {if(sequence_length>*first_sequence_length)
       {printf("\n The sequences are not the same length!!!");
        printf("\n Error, program abnormally terminated");
        printf("\n lengths are %d and at least %d",*first_sequence_length,
             sequence_length);
        printf("\n Abnormal termination\n");
        exit(1);
       }
       if(base!=sequence[sequence_length])
         {c1=toupper(base);
          c2=toupper(sequence[sequence_length]);
          if( !(((c1=='U')&&(c2=='T'))||((c1=='T')&&(c2=='U'))))
            {
             printf("\n\n Warning, At position %d, sequence %s has a %c",
             sequence_length+1,file_data_ct[0],sequence[sequence_length]);
             printf("\n but sequence %s has a %c \n",
                file_data_ct[file],base);
             printf("\n Program continues anyway \n");
            }

         }
     }
    sequence_length++;
    if(helix_started==TRUE)
      {if(column==0)
         {/* finish off helix */  
          diag[*helices_in_plot_file].length=helix_len;
          helix_started=FALSE;
          (*helices_in_plot_file)++;
          /* printf("\n finished Helix %d of length %d",
              helices_in_plot_file,helix_len); */
         }
       else
          if(column<row)  /* col<row has already been added */
             {/* finish off helix */   /* all helices appear twice */
               diag[*helices_in_plot_file].length=helix_len;
               helix_started=FALSE; /* just finish current helix */
               (*helices_in_plot_file)++;
	       /*   printf("\n finished helix of length %d",helix_len); */
              }
           else
            if(column!=(previous_column-1))
	        {/* finish off helix */
		  /*  printf("\n Finished helix of length %d",helix_len); */
                 diag[*helices_in_plot_file].length=helix_len;
                 (*helices_in_plot_file)++;
                 helix_started=FALSE;
                 /* start new helix */
                 previous_column=column;
                 helix_len=1;
                 if(*len<column)
                    *len=column;
		 /* printf("\n new helix at row %d column %d",row,column); */
                 diag[*helices_in_plot_file].row=row;
                 diag[*helices_in_plot_file].column=column;
                 diag[*helices_in_plot_file].diagonal=row+column-1;
                 diag[*helices_in_plot_file].color=file+2;
                 helix_started=TRUE;
                }

            else
                {helix_len++; /* add one to length */
                 previous_column=column;
                }
 
      } 
   else
     { /* helix has not been started */
      if(column>row) /* start a helix */
        {previous_column=column;
         helix_len=1;
         /* printf("\n new helix at row %d column %d",row,column); */
         diag[*helices_in_plot_file].row=row;
         diag[*helices_in_plot_file].column=column;
         diag[*helices_in_plot_file].diagonal=row+column-1;
         diag[*helices_in_plot_file].color=file+2;
         helix_started=TRUE;
         if(*len<column)
            *len=column;
        }
     }
   }
 /* printf("\n len is %d",*len); */
 /* end of while loop */
 /* if a helix was started, record its length */
    if(helix_started==TRUE)
    {diag[*helices_in_plot_file].length=helix_len;
     diag[*helices_in_plot_file].diagonal=row+column-1;
     (*helices_in_plot_file)++;
     /* printf("\n finished heliX %d of length %d",
         helices_in_plot_file,helix_len); */
    }
    if(*helices_in_plot_file==MAXIMUM_HELICES)
         {printf("\n Error!!!___________________________________");
          printf("\n There are more than %d helices in the plot file",MAXIMUM_HELICES);
          printf("\n Increase the constant MAXIMUM_HELICES in ct_boxplot.c");
          printf("\n and recompile....\n");
          exit(1);
         }
   helices_in_this_file=*helices_in_plot_file-helices_in_this_file;
   printf("\n There are %d helices in this ct file with length %d.",
    helices_in_this_file,sequence_length);
   if(file==0)
      *first_sequence_length=sequence_length;
   else
     {if(sequence_length!=*first_sequence_length)
       {printf("\n\n The sequences are not the same length!!!");
        printf("\n First sequence length is %d, other is %d",
                  *first_sequence_length,sequence_length);
        printf("\n Error, program abnormally terminated\n");
        exit(1);
       }
     }
   if(file_length!=sequence_length)
     {printf("\n Error, There are %d records, but the first line",
       sequence_length);
      printf("\n indicates that there should be %d records\n",
             file_length);
      printf("\n Abnormal Termination\n");
      exit(1);
      }
}

void fix_counts(int helices_in_plot_file,int *diag_count,
    struct helix *diag)
{int i;
 for(i=0;i<helices_in_plot_file;i++)
   {diag_count[diag[i].diagonal]++;
   /*    printf("\n diag_count[%d] is now %d",
          diag[i].diagonal,diag_count[diag[i].diagonal]); */
    /* update the count for this diagonal */
   }
} 



void sort_all_helices_set_start(int *diag_count,int *diag_start,
           int helices_in_plot_file,
           struct helix *diag,int len)
{int i;
 printf("\n There are  %d  helices in the ct files.",helices_in_plot_file);
 if(helices_in_plot_file==0)
   {printf("\n \n Error------   No helices_have_been_read\n");
    exit(1);
   }
 /* diag_count counts the total number of entries in each diagonal */
 /* Sort the plot file */
 /* printf("\n just before sort");*/ 
 /* print_helix(helices_in_plot_file);  */
 sort_helices(diag,helices_in_plot_file);
 /* printf("\n just after sort");*/
 /* print_helix(helices_in_plot_file); */
 /*   printf("\n just before set count to zero"); */
    for(i=0;i<=(2*len-1);i++)
       diag_count[i]=0;
    /* fix the count of diagonal for each entry in plot file */

  /*  print_helix(helices_in_plot_file); */
  diag_start[1]=0;
  fix_counts(helices_in_plot_file,diag_count,diag);
  for(i=2;i<=(2*len-1);i++)
   {diag_start[i]=diag_start[i-1]+diag_count[i-1];
   /* printf("\n for diag=%d count is %d start is %d", 
         i,diag_count[i],diag_start[i]); */
   } 
  /* diag_start counts forward into each array for the start of */
 /* each diagonal */
}

/* add overlaps for an individual diagonal named current_diagonal */
void add_overlap_diag(int *diagonal_fill,int current_diagonal,
                     int *diag_start, int *diag_count,
                     struct helix *diag,
                     int number_of_ct_files,int *new_helices_in_plot_file)
{/* diagonal_fill is a diagonal line that is filled to indicate */
  /* how many dots at each position */
  /* diagonal_count will be updated */
  int count,i,starting_position,farthest_position;
  int row,length,current_row,helix_started;
  int start_row_of_helix;
  int closest_position;
  int length_of_current_helix;
  count=diag_count[current_diagonal];
  if(number_of_ct_files>1)
    {if(count<2)
      return;
    }
  else
    {if(count<1)
      return;
    }
  /* add each entry along the diagonal */
  starting_position=diag_start[current_diagonal];
  /* process each helix on the diagonal */
  farthest_position=0;
  closest_position=MAXIMUM_HELICES;
  for(i=0;i<count;i++)
    {/* printf("\n start row is %d, start column is %d",
        diag[starting_position+i].row,diag[starting_position+i].column);*/
     row=diag[starting_position+i].row;
     length=diag[starting_position+i].length;
     /* process a single helix */
     if(farthest_position<(row+length-1))
        farthest_position=row+length-1;
     if(closest_position>row)
       closest_position=row;
     for(current_row=row;current_row<=(row+length-1);current_row++)
       { diagonal_fill[current_row]++;
       /* printf("\n value is %d at row %d of diagonal %d",
           diagonal_fill[current_row],current_row,current_diagonal);*/
       }
    }
  /* The number at each entry along the diagonal indicates some */
  /* overlap ,  0 or 1 is no overlap */
  /* Value equal to number_of_ct_files indicates full overlap */
  /* Other indicates partial overlap */
  /* Make overlaping helices for partial overlap*/
   if(number_of_ct_files>2)
    {helix_started=FALSE;
     length_of_current_helix=0;
     for(i=closest_position;i<=(farthest_position+1);i++)
       {if((diagonal_fill[i]>1)&&(diagonal_fill[i]<number_of_ct_files))
         {/* continue helix or start new helix */
           if(helix_started==TRUE)
              length_of_current_helix++;
           else
	    {length_of_current_helix=1;
             start_row_of_helix=i;
             helix_started=TRUE;
	     /*  printf("\n starting new partial helix for diagonal %d",
                current_diagonal);*/
            }
         }
      else
         {/* end helix if one was started or do nothing*/
           if(helix_started==TRUE)
            {diag_count[current_diagonal]++;
             diag[*new_helices_in_plot_file].row=start_row_of_helix;
             diag[*new_helices_in_plot_file].diagonal=current_diagonal;
             diag[*new_helices_in_plot_file].color=1;
             diag[*new_helices_in_plot_file].column=current_diagonal-
                     start_row_of_helix+1;
             diag[*new_helices_in_plot_file].length=length_of_current_helix;
             (*new_helices_in_plot_file)++;
	     /*  printf("\n added partial helix, there are now %d",*new_helices_in_plot_file);
             printf("\n row %d was col %d  was length was %d",start_row_of_helix,
                 current_diagonal-row+1,length_of_current_helix);*/
             helix_started=FALSE;
             length_of_current_helix=0;
	  }
        
       }
     }
    }
  /* do full overlap last */
  helix_started=FALSE;
  length_of_current_helix=0;
  /* printf("\n farthest_position is %d",farthest_position); */
  for(i=closest_position;i<=(farthest_position+1);i++)
    {if(diagonal_fill[i]==number_of_ct_files)
      {/* continue helix or start new helix */
        if(helix_started==TRUE)
            {length_of_current_helix++;
            /* printf("\n adding 1 to length of full helix to be %d",length_of_current_helix); */
            }
        else
	  {length_of_current_helix=1;
           start_row_of_helix=i;
           helix_started=TRUE;
	   /*  printf("\n starting new full helix for diagonal %d",
              current_diagonal);*/
          }
      }
     else
       {/* printf("\n diagonal fill is %d",diagonal_fill[i]);*/
        if(helix_started==TRUE)
          {diag_count[current_diagonal]++;
           diag[*new_helices_in_plot_file].row=start_row_of_helix;
           diag[*new_helices_in_plot_file].diagonal=current_diagonal;
           diag[*new_helices_in_plot_file].color=0;
           diag[*new_helices_in_plot_file].column=current_diagonal-start_row_of_helix+1;
           diag[*new_helices_in_plot_file].length=length_of_current_helix;
           (*new_helices_in_plot_file)++;
	   /*  printf("\n added full helix, there are now %d",*new_helices_in_plot_file);
           printf("\n row %d was col %d  was length was %d",start_row_of_helix,
                 current_diagonal-row+1,length_of_current_helix);*/
           helix_started=FALSE;
           length_of_current_helix=0;
	  }
        else
	  {/* printf("\n helix_started=FALSE");*/
          }
        
       }
      diagonal_fill[i]=0; /* clear the array for the next pass */
    }

}

		     /* add overlaps for all diagonals */
void  add_overlaps(int *diag_count,
                  int *diag_start,int *helices_in_plot_file,
                   struct helix *diag,int len,int number_of_ct_files)
{int new_helices_in_plot_file;
 int i;
 int *diagonal_fill;
 diagonal_fill=(int *)malloc(len*sizeof(int));
 /* The length of longest diagonal, this was too big */
 /* Should probably use square_root of 2 * len/2+1 */
 /* Then clear it */
 if(diagonal_fill==NULL)
   {printf("\n Error, Insufficient memory to allocate diagonal_fill");
   exit(1);
   }
 for(i=0;i<len;i++)
     diagonal_fill[i]=0;
 /* Save old information */
 new_helices_in_plot_file= *helices_in_plot_file;
 /* add overlap helices */ 
 for(i=1;i<=(2*len-1);i++)
  {/* printf("\n processing overlap for diagonal %d",i);*/
   add_overlap_diag(diagonal_fill,i,diag_start,diag_count,
    diag,number_of_ct_files,&new_helices_in_plot_file);
  }
 /* update for new starts */
 diag_start[1]=0;
 for(i=2;i<=(2*len-1);i++)
   diag_start[i]=diag_start[i-1]+diag_count[i-1];
 *helices_in_plot_file=new_helices_in_plot_file;
 free(diagonal_fill);
}

void read_of_all_ct_files(FILE *fp,int *diag_count,int *diag_start,
           int number_of_ct_files,
           int *helices_in_plot_file,int *len,
           char filename_ct[MAXIMUM_CT_FILES][80],
            struct helix *diag,
            char file_data_ct[MAXIMUM_CT_FILES][80],
           char *sequence_name,int sequence_name_set)
{ char filename[80];
  int first_sequence_len;
  int i,successful_opens;
  successful_opens=0;  
  for(i=0;i<number_of_ct_files;i++)
   {strcpy(filename,filename_ct[i]); 
    printf("\n Trying to open  %s",filename);   
    if((fp=fopen(filename,"r"))==NULL)
      {printf("\n ---------Could not open file %s------------------- ",
                filename);
       strcpy(file_data_ct[i],"File not found: ");
       strcat(file_data_ct[i],filename);
      }
    else
      {successful_opens++;
       read_ct_file(fp,i,len,diag,file_data_ct,
            helices_in_plot_file,&first_sequence_len,sequence_name,
            sequence_name_set,filename);
       fclose(fp);
      }
  }
 if((successful_opens==0)||(*len==0))
   {printf("\n Unable to open any input files successfully\n");
    exit(1);
   }  
 printf("\n Farthest basepair is  %d\n",*len);
 sort_all_helices_set_start(diag_count,diag_start,
      *helices_in_plot_file,diag,*len); 
 add_overlaps(diag_count,diag_start,helices_in_plot_file,
              diag,*len,number_of_ct_files);
 sort_helices(diag,*helices_in_plot_file);
 printf("\n After overlap there are %d helices",
 *helices_in_plot_file);
}
  


