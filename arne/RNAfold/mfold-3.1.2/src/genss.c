

/* genss.c - accesses coordinates in a plt2 file output by Naview and pastes
                them into a ct file to create an ss file for STRED (XRNA)      */

/* by cchang@abj.bio.sunysb.edu on 1/25/96:				*/
/*									*/
/*	Add the function to indicate the begin/end point of a		*/
/*	helix at the fifth column in the output *.ss file:		*/
/*									*/
/*	Fifth col. =							*/
/*		1:	begin helix, 5' side				*/
/*		2:	end helix, 5' side				*/
/*		3:	begin helix, 3' side				*/
/*		4:	end helix, 3' side				*/
/*		8:	begin line					*/
/*		9:	end line					*/
/*		12:	single base pair				*/
/*									*/
/*	Note:	8 & 9 are only put at very beginning and ending of	*/
/*		the whole structure.					*/

#include <stdio.h>
#define  SCF  10   /* scale factor for plt2 coordinates ==> STRED coordinates */

int split_struct();
void genss_errors ();

main (argc, argv)
int argc;
char *argv[];
{
    char  p2_fil[40], ct_fil[40], ss_fil[40];		/* file names */
    char  rec[120], acgu, ci;
    float x, y;
    int   n, nct, bi, bj, hi;
    FILE  *p2f, *ctf, *ssf;				/* file pointers */

    int   flag, begin_by_8, end_by_9, first9, last8;
	int	c1,c5,c6;
	char	c2;
	float	c3,c4;

	c1=-999;	c6=-999;	begin_by_8 = 0;
	first9=-999;	last8=-999;	flag=0;

/* check for correct command line syntax */
    if (argc != 4)
        genss_errors (1, argv[0]);
    strcpy (ct_fil, argv[1]);                          /* retrieve file names */
    strcpy (p2_fil, argv[2]);
    strcpy (ss_fil, argv[3]);                          

/* Open files	*/
    if ((ctf = fopen(ct_fil, "r")) == NULL)                   /* open ct file */
        genss_errors (2, ct_fil);
    if ((p2f = fopen(p2_fil, "r")) == NULL)                 /* open plt2 file */
        genss_errors (2, p2_fil);
    if ((ssf = fopen(ss_fil, "w")) == NULL)               /* open output file */
        genss_errors (2, ss_fil);

/* pass1 begins	*/
    if (fgets(rec, 120, ctf) == NULL)                        /* ct file empty? */
        genss_errors (3, ct_fil);
    sscanf(rec, "%d", &nct);         /* pick up sequence length from 1st line */

    for (n = 1; n <= nct; n++)  {
        if (fscanf(ctf, "%d%*c%c%*d%*d%d%d", &bi, &ci, &bj, &hi) == EOF)
            genss_errors (5, ct_fil);                 /* unexpected EOF (ct) */

	if(c1 == -999 && flag == 0){
		begin_by_8 = 1;
	}

	if(c6 != 0 && bj == 0){
			last8 = bi;
	}

	if(c6 == 0 && bj != 0 && begin_by_8){
			first9 = bi-1;
			begin_by_8 = 0;
	}

	c1=bi;	c6=bj;

    }

	if(flag == 0)
		end_by_9 = 1;
	else
		end_by_9 = 0;

/* pass1 ends */

	rewind(ctf);
	c1=-999;

/* pass2 begins	*/

    printf ("\n Generating STRED file from  %s  and  %s\n\n", ct_fil, p2_fil);

    if (fgets(rec, 120, ctf) == NULL)                        /* ct file empty? */
        genss_errors (3, ct_fil);
    sscanf(rec, "%d", &nct);         /* pick up sequence length from 1st line */

    while (fgets(rec, 120, p2f) != NULL && strstr(rec, "CTX") == 0);
    if (strstr(rec, "CTX") == 0)           /* confirm base labelling record */
        genss_errors (4, p2_fil);

    for (n = 1; n <= nct; n++)  {
        if (fscanf(ctf, "%d%*c%c%*d%*d%d%d", &bi, &ci, &bj, &hi) == EOF)
            genss_errors (5, ct_fil);                  /* unexpected EOF (ct) */
        sscanf (rec, "%*s%f%f%*c%*c%c", &x, &y, &acgu);
        if (acgu != ci)                      /* input files do not correspond */
            genss_errors (6, p2_fil);


	switch((flag=split_struct(c6,bj,c1,bi,c5))){
		case 20:
		case 21:
		case 23:
		case 40:
		case 41:
		case 43:
		case 120:
		case 121:
		case 123:
			c5 = flag/10;	/* last one is an end.	*/
			flag %= 10;
			break;
		default:
			break;
	}

	if(c1 == -999){
		if(flag == 0){
			flag = 8;
			begin_by_8 = 1;
		}
	}
	else{
		if(c1 == first9 && begin_by_8)
			c5 = 9;
		else if(c1 == last8 && end_by_9)
			c5 = 8;

		fprintf(ssf,"%4d %c %8.2f %8.2f %d %4d\n",c1,c2,c3,c4,c5,c6);
	}

	c1=bi;		c2=ci;		c3=x*SCF;
	c4=y*SCF;	c5=flag;	c6=bj;


        if (fgets(rec, 120, p2f) == NULL || strstr(rec, "CTX") == 0)
            if (n != nct)            
                genss_errors (5, p2_fil);   /* too few base labelling records */
    }

	switch(flag){
		case 1:
		case 3:
			flag = 12;
			break;
		case 0:
			flag = 9;
			break;
		default:
			break;
	}

        fprintf (ssf, "%4d %c %8.2f %8.2f %d %4d\n",
						bi, ci, x*SCF, y*SCF, flag, bj);


    close (p2f);                                               /* close files */
    close (ctf);
    close (ssf);

/* pass2 ends	*/

}


int split_struct(lastbj, bj, lastbi, bi, lastflag)
int lastbj;
int bj;
int lastbi;
int bi;
int lastflag;
{

	if(lastbj == 0){
	  if(bj==0)					/* connecting	*/
		return 0;
	  else{
		if(bi < bj)			/* begin helix, 5' side	*/
			return 1;
		else				/* begin helix, 3' side	*/
			return 3;
	  }
	}
	else if(bj==0){
		if(lastflag == 1 || lastflag == 3)
			return 120;	/* the last is a single pair	*/
		else if(lastbi < lastbj)/* the last is end helix, 5'	*/
			return 20;
		else			/* the last is end helix, 3'	*/
			return 40;
	}
	else if(lastbj-bj == 1 || bj-lastbj == 1)	/* base pair	*/
		return 0;
	else if(bi < bj){
		if(lastflag == 1 || lastflag == 3)
			return 121;	/* last is single, this is 1	*/
		else if(lastbi < lastbj)/* the last is 2 then this is 1	*/
			return 21;
		else			/* the last is 4 then this is 1	*/
			return 41;
	}
	else {
		if(lastflag == 1 || lastflag == 3)
			return 123;	/* last is single, this is 3	*/
		else if(lastbi < lastbj)/* the last is 2 then this is 3	*/
			return 23;
		else			/* the last is 4 then this is 3	*/
			return 43;
	}
}



void genss_errors (errno, name)
int errno;
char *name;             /* report error and exit */
{
    switch (errno)  {
    case 1:                                  /* incorrect command line syntax */
        printf ("\n *** Usage: %s ct_fil plt2_fil ss_fil\n\n", name);
        exit(1);
    case 2:                                                    /* open failed */
        printf ("\n *** Open failed on file %s\n\n", name);
        exit(1);
    case 3:                                                     /* file empty */
        printf ("\n *** Input file %s is empty\n\n", name); 
        exit(1);
    case 4:                                      /* no base labelling records */
        printf ("\n *** No base labelling records in file %s\n\n", name); 
        exit(1);
    case 5:                                                 /* unexpected EOF */
        printf ("\n *** Sequence incomplete in %s\n\n", name);
        exit(1);
    case 6:                                  /* input files do not correspond */
        printf ("\n *** Files contain different sequences\n\n");
        exit(1);
    }
}
