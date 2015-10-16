#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#ifndef WHITE_SPACE
#define WHITE_SPACE " \t\r\n"
#endif
#define SEEK_START 0

/*C code is 20x faster than the python version with O3 compiling option*/
int min_int(int *ar, int n){
    int i;
    int min = ar[0];
    for (i=1;i<n;i++){
        if (ar[i] < min) {
            min = ar[i];
        }
    }
    return min;
}
int max_int(int *ar, int n){
    int i;
    int max = ar[0];
    for (i=1;i<n;i++){
        if (ar[i] > max) {
            max = ar[i];
        }
    }
    return max;
}
void PrintFasta(char **idList, char** annoList, char **seqList, int numseq, FILE *fpout){/*{{{*/
    int i;
    for (i = 0; i < numseq; i ++){
        int seqlen = strlen(seqList[i]);
        fprintf(fpout, ">%s\n", annoList[i]);
        fprintf(fpout, "%s\n", seqList[i]);
    }
}/*}}}*/
int IsInCharSet(const char ch, const char *charSet){/*{{{*/
/*****************************************************************************
 * check if the character "ch" is in charSet,
 ****************************************************************************/
    int n = strlen(charSet);
    int i;
    for(i = 0 ;i < n ; i ++) {
        if(ch == charSet[i]) {
            return 1;
        }
    }
    return 0;
}/*}}}*/
char* get_input_into_string(char* string,char* infile){/*{{{*/ 
    int i = 0;
    int string_length = 2;
    char c = 0;
    FILE *file = 0;
    int freadReturnValue = 0;
    if(infile){
        if (!(file = fopen( infile, "r" ))){
            return 0;
            fprintf(stderr,"Cannot open file '%s'\n", infile);
            exit(1);
        }
        if (fseek(file,0,SEEK_END) != 0){
            (void)fprintf(stderr, "ERROR: fseek failed\n");
            (void)exit(EXIT_FAILURE);
        }
        i= ftell (file); /*using fseek to get the number of chars of the file*/
        if (fseek(file,0,SEEK_START) != 0){
            (void)fprintf(stderr, "ERROR: fseek failed\n");
            (void)exit(EXIT_FAILURE);
        }
        string = (char*) malloc ((i+1)* sizeof(char));
        if ((freadReturnValue=fread(string,sizeof(char), i, file))!= i) {
            fprintf(stderr,"fread error at line %d in file %s\n", __LINE__, __FILE__);/*added 2010-09-27*/
        }
        string[i] = 0;
        fclose(file);
    }else{  /*stdin input*/
        if (!isatty(0)){
            string = (char*) malloc(sizeof(char*)*string_length);
            while (!feof (stdin)){
                c = getc(stdin);
                if (i == string_length){
                    string_length <<=1;
                    string = realloc(string,sizeof(char)*string_length); /*dynamic allocatio of the string*/
                }
                string[i] = c;
                i++;
            }
            string[i-1] = 0;
        }else{
            return 0;
        }
    }
    return string;
}/*}}}*/
int count_sequences_fasta(char* string){/*{{{*/ 
    int nbytes;
    int i; 
    int n = 0;
    int stop = 0;
    nbytes = strlen(string);
    for (i =0;i < nbytes;i++){
        if (string[i] == '>'&& stop == 0){
            stop = 1;
            n++;
        }
        if (string[i] == '\n'){
            stop = 0;
        }
    }
    if(!n){
        return 0;
    }
    return n;
}/*}}}*/
int FreeString(char *string){/*{{{*/
    if (string != NULL){
        free(string);
    }
    return 0;
}/*}}}*/
char *SpanExcluding(const char* strForSpan,char* strAfterSpan, const char charSet[]/*=WHITE_SPACE*/)/*{{{*/
//**********************************************************************
// SpanExcluding()
// compress the string by removing the characters in the charSet
// and store the compressed string into strAfterSpan
// the search is case-sensitive
// default charSet is WHITE_SPACE
//**********************************************************************
{
    int l = strlen(strForSpan);
    int n = strlen(charSet);
    int i,j ; 
    j = 0 ;
    for(i = 0 ; i < l ;i ++)
    {
        if(!IsInCharSet(strForSpan[i],charSet))
        {
            strAfterSpan[j] = strForSpan[i] ;
            j ++ ;
        }
    }
    strAfterSpan[j] = '\0' ;
    return strAfterSpan;
}/*}}}*/
char *GetIDFromAnnotationLine(char *line, char *seqid){/*{{{*/
    char *tmpstr = (char*) malloc(sizeof(char)*(strlen(line)+1));
    sscanf(line, " %100s ", tmpstr);
    char * pchp;
    char * pch;
    pchp = &tmpstr[0]-1;
    pch = strchr(tmpstr,'|');
    while (pch != NULL) {
        pchp = pch;
        pch = strchr(pch+1,'|');
    }
    strcpy(seqid, pchp + 1);
    free(tmpstr);
    return seqid;
}/*}}}*/
int ReadFastaFromBuffer(char* buffer, char **idList, char **annoList, char /*{{{*/ 
        **seqList){
    char *str = buffer;
    int buffersize = strlen(buffer);
    char *pch_beg = NULL;
    char *pch_end = NULL;
    char *pch_seqbeg = NULL;
#ifdef DEBUG_READFASTA
    int cntSeq = 0;
    fprintf(stdout, "Reading sequences ...\n");
#endif
    int cntseq = 0;
    while (1){
        pch_beg = strchr(str, '>'); 
        if (pch_beg == NULL){
            break;
        }
        pch_end = strstr(pch_beg+1, "\n>");
        if (pch_end == NULL){
            pch_end = buffer + buffersize - 1;
        }

        pch_seqbeg = strchr(pch_beg+1,'\n');
        if (pch_seqbeg == NULL){
            fprintf(stderr, "Fatal! fasta seqfile format error. Exit.\n");
            exit(1);
        }
        char *seqid = NULL;
        char *annoline = NULL;
        char *seq = NULL;
        char *rawseq = NULL;

        int sizeAnnoLine = pch_seqbeg - pch_beg - 1;
        int sizerawseq = pch_end - pch_seqbeg - 1;
        seqid = (char*) malloc(sizeof(char)*(sizeAnnoLine+1));
        annoline = (char*) malloc(sizeof(char)*(sizeAnnoLine+1));
        rawseq = (char*) malloc(sizeof(char)*(sizerawseq+1));
        seq = (char*) malloc(sizeof(char)*(sizerawseq+1));


        strncpy(annoline, pch_beg+1, sizeAnnoLine);
        annoline[sizeAnnoLine] = 0;
        strncpy(rawseq, pch_seqbeg+1, sizerawseq);
        rawseq[sizerawseq] = 0;
        GetIDFromAnnotationLine(annoline, seqid);

        SpanExcluding(rawseq, seq, WHITE_SPACE);

        idList[cntseq] = (char*) malloc(sizeof(char)*(strlen(seqid)+1));
        annoList[cntseq] = (char*) malloc(sizeof(char)*(strlen(annoline)+1));
        seqList[cntseq] = (char*) malloc(sizeof(char)*(strlen(seq)+1));
        strcpy(idList[cntseq], seqid);
        strcpy(annoList[cntseq], annoline);
        strcpy(seqList[cntseq], seq);

        cntseq += 1;

        free(seqid);
        free(annoline);
        free(rawseq);
        free(seq);

#ifdef DEBUG_READFASTA
        //cout << ">" <<seqid<<endl;
        //cout << seq << endl;
        if (cntseq % 10000 == 1){
            fprintf(stderr,"%d...\n", cntseq);
            fflush();
        }
#endif
        str = pch_end;
        if (str - buffer >= buffersize){
            break;
        }
    }
#ifdef DEBUG_READFASTA
        cout << "Finished." <<endl;
        cout << "Total sequences read in: " << cntSeq << endl;
#endif
    return 0;
}/*}}}*/
char ** RemoveUnnecessaryGap(char ** seqList, int numseq){
    int i,j;
    int *lengthList = (int *)malloc (sizeof(int)*numseq);
    for (i=0; i< numseq;i++){
        lengthList[i] = strlen(seqList[i]);
    }
    int min = min_int(lengthList, numseq);
    int max = max_int(lengthList, numseq);
    if (min != max){
        fprintf(stderr, "length uneqal, min = %d, max = %d\n", min, max);
        exit(1);
    }
    int lengthAln = strlen(seqList[0]);
    int *cntgapList = (int*) malloc(sizeof(int)*lengthAln);
    for (j=0; j< lengthAln;j++){
        int cnt = 0;
        for (i=0;i<numseq;i++){
            if (seqList[i][j] == '-'){
                cnt += 1; 
            }
        }
        cntgapList[j] = cnt;
    }
    char ** newSeqList = (char **) malloc(sizeof(char*)*numseq);
    char *seq = (char*) malloc(sizeof(char)*(lengthAln+1));
    for (i=0;i<numseq;i++){
        int cnt = 0;
        for (j=0;j<lengthAln;j++){
            if (cntgapList[j] < numseq){
                seq[cnt] = seqList[i][j];
                cnt += 1;
            }
        }
        seq[cnt] = '\0';
        newSeqList[i] = (char*)malloc(sizeof(char)*(cnt+1));
        strcpy(newSeqList[i], seq);
    }
    FreeString(seq);
    free(cntgapList);
    free(lengthList);
    return newSeqList;
}
int main(int argc, char** argv) {/*{{{*/
    int isNonOptionArg = 0;
    if(argc < 2) {
        return 1;
    }
    int i = 0;
    int j = 0;
    char *seqfile = NULL;
    char *outfile = NULL;
    const char control_option[] = ""; 
    i = 1;
    while(i < argc) {/*{{{*/
        if(argv[i][0] == '-' && !isNonOptionArg){ /*options*/
            isNonOptionArg = 0;
            if(IsInCharSet(argv[i][1], control_option)){
                /*if argv[i][1] is in control_option, it might be used as
                 * -aqs*/
                for(j = 1 ; j < strlen(argv[i]); j++) {
                    switch (argv[i][j]) {
                        default : fprintf(stderr,"Invalid option, non-control option '%c' can be used together with contorl-option\n", argv[i][j]);
                                  return -1;
                    }
                }
                i ++;
            } else if( (strcmp(argv[i],"-o") == 0) || (strcmp(argv[i],
                            "--out") == 0))  {
                outfile = (char*) malloc(sizeof(char) * (strlen(argv[i+1])+1));
                strcpy(outfile, argv[i+1]);
                i += 2;
            } else if (strcmp(argv[i], "--") == 0) {
                /*next item is non option argument*/
                isNonOptionArg = 1;
                i ++;
                continue;
            } else {
                fprintf(stderr,"Error! Invalid argument '%s'\n", argv[i]);
                return -1;
            }
        } else {/*non-option argument*/
            seqfile = (char*) malloc(sizeof(char) * (strlen(argv[i])+1));
            strcpy(seqfile, argv[i]);
            i ++;
        }
    }/*}}}*/

/*Step 1, Read in Fasta Seq*/
    char *buffer = NULL;
    int numseq = 0;
    buffer = get_input_into_string(buffer, seqfile);
    if (buffer == NULL){
        fprintf(stderr, "Failed to read file %s\n", seqfile);
        exit(1);
    }
    numseq = count_sequences_fasta(buffer);
    char ** idList = NULL;
    char ** annoList = NULL;
    char ** seqList = NULL;
    annoList = (char**) malloc(sizeof(char*) * numseq);
    idList = (char**) malloc(sizeof(char*) * numseq);
    seqList = (char**) malloc(sizeof(char*) * numseq);
    ReadFastaFromBuffer(buffer, idList, annoList, seqList);
    free(buffer);
/*Step 2: Remove Unnecessary gap*/
    char **newSeqList = NULL;
    newSeqList = RemoveUnnecessaryGap(seqList, numseq);

/* Step 3: print the removed uncessary gap*/
    FILE *fpout = stdout;
    if (outfile != NULL){
        fpout = fopen(outfile,"w");
        if( fpout == NULL){
            fprintf(stderr, "Failed to write to file %s, reset to STDOUT\n", outfile);
            fpout = stdout;
        }
    }
    PrintFasta(idList, annoList, newSeqList, numseq, fpout);

    for (i = 0; i < numseq; i ++){
        FreeString(idList[i]);
        FreeString(annoList[i]);
        FreeString(seqList[i]);
        FreeString(newSeqList[i]);
    }
    free(idList);
    free(seqList);
    free(annoList);
    free(newSeqList);
    FreeString(seqfile);
    FreeString(outfile);
    if(fpout != NULL && fpout != stdout) fclose(fpout);

} /*}}}*/
