#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <algorithm>

using namespace std;
#ifndef WHITE_SPACE
#define WHITE_SPACE " \t\r\n"
#endif
/*ChangeLog
 * ChangLog 2011-11-04
 *      For seqMSAFile, sequences are read in by block
 * ChangeLog 2012-03-28
 *      Add option -ignore-badseq y|n, this is useful when 
 *      the input fasta file is a pairwise sequence alignment
 *      then no sequence should be removed, otherwise, pair-matching will be
 *      broken*/

/*Valgrind checked at 
 * 2011-10-29: no errors and no memory leaking*/

int BLOCK_SIZE = 100000;
int outFormat=0;
char usage[]="\n\
usage: matchMSAtopo -msa seqMSAFile -topo topoFile\n\
\n\
Description: Output the topology alignment given sequence alignment.\n\
             All files should be in Fasta format.  Sequences are matched by\n\
             seqID identified from the annotation line\n\
\n\
  -o| FILE  Output the result to FILE, (default: stdout)\n\
  -of INT   Set output format, (default: 0)\n\
            0: Fasta format \n\
            1: Aligned format\n\
               each line with seqID: seq\n\
  -ignore-badseq y|n\n\
            whether ignore bad seq matching (default: yes)\n\
  -h|--help Print this help message and exit\n\
\n\
Created on 2010-08-01, updated 2012-03-28, Nanjiang Shu\n\
";
void PrintHelp() {
    cout << usage << endl;
}
void PrintVerboseHelp() { 
}

int GetGaplessLength(const char *seq, int n){/*{{{*/
    int cnt = 0;
    for (int i = 0; i < n ; i ++){
        if (seq[i] != '-'){
            cnt ++;
        }
    }
    return cnt;
}/*}}}*/
bool IsFileExist(const char *name)/*{{{*/
{
    struct stat st;
    return (stat(name,&st) == 0); 
}/*}}}*/ 
bool IsInCharSet(const char ch, const char *charSet, int n /*= 0 */)/*{{{*/
/*****************************************************************************
 * check if the character "ch" is in charSet,
 ****************************************************************************/
{
	if(n == 0)
        n = strlen(charSet);
    int i;
	for(i = 0 ;i < n ; i ++)
	{
		if(ch == charSet[i])
			return true;
	}
	return false;
}/*}}}*/
int GetFileSize(const char *name)/*{{{*/
{
    struct stat st;
    if (stat(name,&st) == 0){
        return st.st_size; 
    }else{
        cerr << "File " 
            << name
            << " does not exist." << endl;
        return -1;
    }
}/*}}}*/ 
string int2string(int number)/*{{{*/
{
    stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
}/*}}}*/
int my_strcpy(char* to, const char* from, int max)/*{{{*/
/******************************************************************************
 * my modification of strcpy
 * copy max characters from "from" to "to", add NULL terminator automatically
 * updated 2008-04-23, memcpy win in speed when the copying string is long
 * updated 2011-10-27:
 *   since calling strlen will be time consuming for very large from string,
 *   e.g. when "from" is the buffer of a whole trunk of file. Therefore, strlen
 *   part is removed.
 *****************************************************************************/
{
    if(max < 200) {
        strncpy(to,from,max);
    } else {
        memcpy(to,from,max);
    }
    to[max] = '\0';
    return max;
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
		if(!IsInCharSet(strForSpan[i],charSet, n))
		{
			strAfterSpan[j] = strForSpan[i] ;
			j ++ ;
		}
	}
	strAfterSpan[j] = '\0' ;
	return strAfterSpan;
}/*}}}*/
bool IsNumeric( const char* pszInput, int nNumberBase = 10)/*{{{*/
{
    istringstream iss( pszInput );
 
    if ( nNumberBase == 10 ) {
        double dTestSink;
        iss >> dTestSink;
    } else if ( nNumberBase == 8 || nNumberBase == 16 ) {
        int nTestSink;
        iss >> ( ( nNumberBase == 8 ) ? oct : hex ) >> nTestSink;
    } else{
        return false;
    }
    // was any input successfully consumed/converted?
    if ( ! iss ) {
        return false;
    }
    // was all the input successfully consumed/converted?
    return ( iss.rdbuf()->in_avail() == 0 );
}/*}}}*/
int option_parser_filename(int argc, char **argv, int beg, string &filename)/*{{{*/
/*********************************************************************
 * beg is the current index of argument list, e.g., argv[i] = "--out"
 **********************************************************************/
{
    int i ; 
    bool isNonOptionArg = false;
    bool isFileNameSet = false;

    for(i = beg +1 ; i < argc ; i++) {
        if(argv[i][0] == '-' && string(argv[i])!= "--" && !isNonOptionArg) {
            cerr << "option \'" 
                << argv[beg] 
                << "\' must be followed by a filename, not option."
                << endl;
            return -1;
        } else if(string(argv[i])== "--" && !isNonOptionArg) {
            isNonOptionArg = true;
        } else {
            filename = argv[i];
            isFileNameSet = true;
            break;
        }
    }
    if(!isFileNameSet) {
        cerr << "option \'" 
            << argv[beg] 
            << "\' must be followed by a filename."
            << endl;
        return -1;
    } else {
        return i+1;
    }
}
/*}}}*/
template <class T> int option_parser_numeric(int argc, char **argv, int beg, T &x, bool isRangeSet /*= false*/, T min /*= MIN_INT*/, T max /*= MAX_INT*/)/*{{{*/
/************************************************************************
 * beg is the current index of argument list, e.g., argv[i] = "--value"
 ***********************************************************************/
{
    int i ; 
    bool isValueSet = false;
    double tmp;
    i = beg +1;

    if (i < argc) {
        if(IsNumeric(argv[i])) {
            tmp = atof(argv[i]);
            if(isRangeSet) {
                if(tmp < min || tmp > max) {
                    cerr << "Invalid value! Value after option \'"
                         <<  argv[beg]
                         <<  "\' must be in the range of ["
                         << min << ", " 
                         << max << "]"
                         << endl;
                    return -1;
                }
            }
            x = T(tmp);
            isValueSet = true;
        }
    }

    if(!isValueSet) {
        cerr << "option \'" 
            << argv[beg] 
            << "\' must be followed by a numerical value."
            << endl;
        return -1;
    } else {
        return i+1;
    }
}
/*}}}*/
string GetIDFromAnnotationLine_old(const char *line, string &id)/*{{{*/
{
    char tmpstr[100+1] ="";
    sscanf(line, ">%100s ", tmpstr );
    char * pchp;
    char * pch;
    pchp = &tmpstr[0]-1;
    pch = strchr(tmpstr,'|');
    while (pch != NULL) {
        pchp = pch;
        pch=strchr(pch+1,'|');
    }
    id = pchp+1;
    return id;
}/*}}}*/
string GetIDFromAnnotationLine(const char *line, string &id)/*{{{*/
{ /*return the first word delimited by whitespace or |*/
    int sizeLine = strlen(line);
    char *newstr  = new char [sizeLine+1];
    strcpy(newstr,line);
    char *pstr = newstr;
    if (pstr[0] == '>'){
        pstr ++;
    }
    char delim[]= " \r\t\n";
    char *pch = NULL;
    pch = strtok(pstr, delim);
    while (pch){
        if (strcmp(pch, "") != 0){
            break;
        }
        pch = strtok(NULL, delim);
    }
    char *firstword = new char [sizeLine+1]; 
    strcpy(firstword, pch);
    pch = strchr(firstword, '|');
    if (pch != NULL){
        char *tmpstr0 = new char[sizeLine + 1];
        char *tmpstr1 = new char[sizeLine + 1];
        char *pch1 = NULL;
        pch1 = strtok(pstr, "|");
        int cnt = 0 ;
        while (pch1){
            if (cnt == 0){
                strcpy(tmpstr0, pch1);
            } else if (cnt == 1){
                strcpy(tmpstr1, pch1);
            }else{
                break;
            }
            cnt += 1;
            pch1 = strtok(NULL, "|");
        }
        if (strcmp(tmpstr0, "sp") == 0 
                || strcmp(tmpstr0, "lcl") == 0
                || strcmp(tmpstr0, "tr") == 0
                || strcmp(tmpstr0, "gi") == 0
                || strcmp(tmpstr0, "r") == 0
                || strcmp(tmpstr0, "p") == 0
                ) {
            id = tmpstr1;
        }else{
            id = tmpstr0;
        }

        delete [] tmpstr0;
        delete [] tmpstr1;
    }else {
        id = firstword;
    }
    delete [] newstr;
    delete [] firstword;
    return id;
}/*}}}*/
int GetFastaSeqFromBuffer(char* buffer, vector <string> &idList, vector <string> & annoList, vector <string> &seqList, int &seq_type)/*{{{*/
{
    char *str = buffer; 
    int buffersize = strlen(buffer);
    char *pch_beg = NULL;
    char *pch_end = NULL;
    char *pch_seqbeg = NULL;
#ifdef DEBUG_READFASTA
    int cntSeq = 0;
    cout << "Reading sequences ..." << endl;
#endif
    while (1){
        pch_beg = strchr(str, '>'); 
        if (pch_beg == NULL){
            break;
        }
        pch_end = strstr(pch_beg+1, "\n>");
        if (pch_end == NULL){
            pch_end = buffer + strlen(buffer) - 1;
        }
        int sizerecord = pch_end - pch_beg + 1 ;  
        char *seqWithAnno = new char [sizerecord+1];
        my_strcpy(seqWithAnno, pch_beg, sizerecord);
#ifdef DEBUG_GETSEQFROMBUFFER
        cout << "seqWithAnno[" << cntSeq << "]=" << seqWithAnno << endl;
#endif
        string anno = "";
        char *pch_anno_end = strchr(seqWithAnno, '\n');
        if (pch_anno_end != NULL){
            *pch_anno_end = '\0'; /*assign a value should used *pch = something*/
            anno = seqWithAnno+1;
        }
#ifdef DEBUG_GETSEQFROMBUFFER
        cout << "Anno[" << cntSeq << "]=" << anno << endl;
#endif

        string seqid = "";
        GetIDFromAnnotationLine(seqWithAnno, seqid);
        //cout << "<" << seqid << ">" << endl;

        pch_seqbeg = strchr(pch_beg+1,'\n');
        if (pch_seqbeg == NULL){
            cerr << "Fatal! fasta seqfile format error. Exit." << endl;
            exit(1);
        }

        int sizerawseq = pch_end - pch_seqbeg;
        my_strcpy(seqWithAnno, pch_seqbeg+1, sizerawseq);
#ifdef DEBUG_GETSEQFROMBUFFER
        cout << "seq[" << cntSeq << "]=" << seqWithAnno << endl;
#endif

        char *seq = new char [sizerawseq+1];
        SpanExcluding(seqWithAnno,seq , WHITE_SPACE);
        idList.push_back(seqid);
        annoList.push_back(anno);
        seqList.push_back(string(seq));
#ifdef DEBUG_READFASTA
        cntSeq += 1;
        //cout << ">" <<seqid<<endl;
        //cout << seq << endl;
        if (cntSeq % 10000 == 1){
            cout << cntSeq << "..." << flush;
            //cout << "position = " << str - buffer << endl;
        }
#endif
        delete [] seqWithAnno;
        delete [] seq;
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
int ReadFasta(string &file,  vector <string> &idList, vector <string> &annoList, vector <string> &seqList, int &seq_type)/*{{{*/
{
    int filesize = GetFileSize(file.c_str());
#ifdef DEBUG_READFASTA
    cout << "file " << file << ". size = " << filesize << endl;
#endif 
    int MAX_FASTA_FILE_SIZE = 1*1024*1024*1024;
    if (filesize <= 0) {
        cerr << "Size of the file " << file 
            << " <= 0." 
            << endl;
        return -1;
    } else if (filesize > MAX_FASTA_FILE_SIZE){
        cerr << "Size of the file " << file 
            << " is over the limit (" 
            << MAX_FASTA_FILE_SIZE << ")."  << endl;
        return -1;
    }else{
        FILE *fpin = fopen(file.c_str(), "rb");
        if (fpin != NULL){
            char *buffer = new char [filesize +1];
            int nread = fread(buffer, sizeof (char), filesize, fpin);
#ifdef DEBUG_READFASTA
            cout << "buffer of size " << filesize+1 << " allocated." << endl;
            cout << "Read in " << nread << " bytes from file " << file << endl;
#endif 
            fclose(fpin);
            if (nread != filesize) {
                cerr << "Read file " << file <<
                    " failed. nread (" << nread << ") != filesize (" <<
                    filesize << ")." << endl;
                return -1;
            }
            buffer[filesize] = '\0';
            GetFastaSeqFromBuffer( buffer, idList, annoList, seqList, seq_type);
            delete [] buffer ; 
            return idList.size();
        }else{
            cerr << "Failed to read fasta file " << file << endl;
            return -1;
        }
    }
}/*}}}*/
int ExtractFromSeqWithAnno(const char* seqWithAnno, int sizeSeqWithAnno, vector <string> &idList, vector <string>&annoList, vector <string>&seqList, int &seq_type){

    /*working even with enhanced fasta with auxiliary info enclosed by {}*/
    const char *pstr = seqWithAnno;
    const char *posAnnoEnd = strchr(pstr, '\n');
    if (posAnnoEnd == NULL){
        return -1;
    }
    int sizeAnno = posAnnoEnd - pstr -1;
    char *anno = new char[sizeAnno +1];
    my_strcpy(anno, pstr+1, sizeAnno);
    annoList.push_back(string(anno));

    string seqid;
    GetIDFromAnnotationLine(anno, seqid);
    idList.push_back(seqid);
    delete [] anno;

    int sizeSeq = sizeSeqWithAnno - sizeAnno;
    if (sizeSeq <= 0 ){
        cerr << "Sequence error: " << seqWithAnno << endl;
    } else {
        char *seq = new char[sizeSeq + 1];
        const char *pch = posAnnoEnd+1;
        bool isAuxiliaryRegion = false;
        char *pSeq = seq;
        while (*pch){
            if (!isAuxiliaryRegion){
                if(! (*pch == '\n' || *pch == ' ' || *pch == '\r' || *pch == '\t') ){
                    *pSeq = *pch;
                    pSeq ++;
                }else if (*pch == '{'){
                    isAuxiliaryRegion= true;
                }else if (*pch == '}'){
                    isAuxiliaryRegion = false;
                }
            }
            pch ++;
        }
        *pSeq = '\0';
        seqList.push_back(string(seq));
        delete [] seq;
    }

    return 0;
}

int ReadFastaFromBuffer(const char *buff,  vector <string> &idList, vector <string> &annoList, vector <string> &seqList, int &seq_type, string& unprocessedBuffer, bool isEOFreached)/*{{{*/
{
    /*Read Fasta sequence from buffer, since buffer might be obtained by BLOCK
     * read, there might be  unprocessedBuffer, unprocessedBuffer will be
     * renewed*/
    string brokenSeqWithAnnoLine="";
    const char *pstr = buff;
    const char *pbeg = 0;
    const char *pend = 0;
    while (1){
        pbeg = strchr(pstr, '>');
        pend = strstr(pstr+1, "\n>");
        if (pbeg != NULL){
            if (pend  != NULL){
                int sizeSeqWithAnno = pend - pbeg;
                char *seqWithAnno = new char [sizeSeqWithAnno +1];
                strncpy(seqWithAnno, pbeg, sizeSeqWithAnno);
                seqWithAnno[sizeSeqWithAnno] = '\0';
                if(ExtractFromSeqWithAnno(seqWithAnno, sizeSeqWithAnno, idList, annoList, seqList, seq_type) != 0){
                    fprintf(stderr,"Sequence error: %s\n", seqWithAnno);
                }
                delete [] seqWithAnno;
                pstr = pend+1;
            } else{
                brokenSeqWithAnnoLine = pbeg;
                break;
            }
        }else{
            break;
        }
    }
    if (isEOFreached && brokenSeqWithAnnoLine != "" ){
        if(ExtractFromSeqWithAnno(brokenSeqWithAnnoLine.c_str(), brokenSeqWithAnnoLine.length(), idList, annoList, seqList, seq_type) != 0){
            fprintf(stderr,"Sequence error: %s\n", brokenSeqWithAnnoLine.c_str());
        }
        brokenSeqWithAnnoLine="";
    }
    unprocessedBuffer = brokenSeqWithAnnoLine;
    return idList.size();
}/*}}}*/

int main(int argc, char** argv)/*{{{*/
{
    bool isNonOptionArg = false;

    if(argc < 2) {
        PrintHelp();
        return 1;
    }
    int i,j;
    string msaFile = "";
    string topoFile = "";
    string seqFile = "";
    string outfile = "";
    bool isIgnoreBadseq = true;
    const char control_option[] = ""; //options which control the program, and does not take parameters
    i = 1;
    while(i < argc)/*{{{*/ {
        if(argv[i][0] == '-' && !isNonOptionArg) {
            isNonOptionArg = false;
            if(IsInCharSet(argv[i][1], control_option, 0))//if argv[i][1] is in control_option, it might be used as -aqs
            {
                for(j = 1 ; j < int(strlen(argv[i])); j++) {
                    switch (argv[i][j]) {
                        default : fprintf(stderr,"Invalid option, non-control option '%c' can be used together with contorl-option\n", argv[i][j]); return -1;
                    }
                }
                i ++;
            } else if(strcmp(argv[i],"-h") == 0 ||strcmp(argv[i],"--help")==0 ) {
                PrintHelp(); 
                return 0;
            } else if(strcmp(argv[i],"-H") == 0 ) {
                PrintVerboseHelp();
                return 0;
            } else if( (strcmp(argv[i],"-o") == 0) 
                    || (strcmp(argv[i], "--o") == 0)
                    || (strcmp(argv[i], "-out") == 0)
                    || (strcmp(argv[i], "--out") == 0)
                    || (strcmp(argv[i], "-outfile") == 0)
                    )  
            {
                if( ( i = option_parser_filename(argc, argv, i, outfile)) == -1) {
                    return -1;
                }
            } else if( (strcmp(argv[i],"-fa") == 0) || (strcmp(argv[i], "--fa") == 0))  {
                if( ( i = option_parser_filename(argc, argv, i, seqFile)) == -1) {
                    return -1;
                }
            } else if( (strcmp(argv[i],"-topo") == 0) || (strcmp(argv[i], "--topo") == 0))  {
                if( ( i = option_parser_filename(argc, argv, i, topoFile)) == -1) {
                    return -1;
                }
            } else if( (strcmp(argv[i],"-msa") == 0) || (strcmp(argv[i], "--msa") == 0))  {
                if( ( i = option_parser_filename(argc, argv, i, msaFile)) == -1) {
                    return -1;
                }
            } else if( (strcmp(argv[i],"-ignore-badseq") == 0) ||
                    (strcmp(argv[i], "--ignore-badseq") == 0))  {
                string yesorno;
                if( ( i = option_parser_filename(argc, argv, i, yesorno)) == -1) {
                    return -1;
                }
                const char *ptr = yesorno.c_str();
                char c = toupper(ptr[0]);
                if ( c == 'Y'){
                    isIgnoreBadseq = true;
                }else{
                    isIgnoreBadseq = false;
                }
            } else if( (strcmp(argv[i], "--format") == 0) 
                    || (strcmp(argv[i],"-f")==0)
                    || (strcmp(argv[i],"-of")==0)
                    || (strcmp(argv[i],"--of")==0)
                    )  {
                if( ( i = option_parser_numeric(argc, argv, i, outFormat, true, 0, 1)) == -1) {
                    return -1;
                }
            } else if (strcmp(argv[i], "--") == 0)/*next item is non option argument*/ {
                isNonOptionArg = true;
                i ++;
                continue;
            } else {
                fprintf(stderr,"Error! Invalid argument '%s'\n", argv[i]);
                return -1;
            }
        } else /*non-option argument*/ {
            i ++;
        }
    }/*}}}*/

    FILE *fpout = NULL;
    if(outfile == "" || strcasecmp(outfile.c_str(), "stdout") == 0) {
        fpout = stdout;
    } else {
        fpout = fopen(outfile.c_str(), "w");
        if (fpout == NULL){
            cout << "Failed to write to file " << outfile << endl;
            fpout = stdout;
        }
    }
    /*add code here*/

    vector <string> idList_Topo;
    vector <string> annoList_Topo;
    vector <string> seqList_Topo;
    int numSeq_Topo = 0;
    int seq_type_Topo = 0;
    //numSeq_MSA = ReadFasta(msaFile, idList_MSA, annoList_MSA, seqList_MSA, seq_type_MSA) ;
    numSeq_Topo = ReadFasta(topoFile, idList_Topo, annoList_Topo, seqList_Topo, seq_type_Topo) ;
#ifdef DEBUG
    cout << "numSeq_Topo = " << numSeq_Topo << endl;
#endif
    map <string, string > topoMap;
    for (i = 0; i < numSeq_Topo; i ++){
        topoMap.insert(pair<string,string>( idList_Topo[i], seqList_Topo[i]));
    }
#ifdef DEBUG
    cout << "Start matching " << numSeq_MSA << " alinged sequences..." << endl;
#endif
    int maxIDsize = 10;
    if (outFormat == 1){
        vector <unsigned int> idsize;
        vector <string>:: iterator  iss;
        for (iss = idList_Topo.begin(); iss != idList_Topo.end();iss ++){
            idsize.push_back((*iss).size());
        }
        maxIDsize = *max_element(idsize.begin(),idsize.end());
        cout << "maxIDsize=" << maxIDsize << endl;
    }

    FILE *fpAln = fopen (msaFile.c_str(), "rb");
    if (fpAln == NULL){
        cerr << "Failed to open msaFile " << msaFile << ". Exit." << endl;
        return -1;
    }
    string unprocessedBuffer = "";
    bool isEOFreached = false;
    int cnt = 0;
    
    char *buff=new char [BLOCK_SIZE +1];
    while (1){
        vector <string> idList_MSA;
        vector <string> annoList_MSA;
        vector <string> seqList_MSA;
        int numSeq_MSA = 0;
        int seq_type_MSA = 0;

        int nread = fread(buff, sizeof (char), BLOCK_SIZE, fpAln);
        buff[nread] = '\0';
        if (nread < BLOCK_SIZE) {
            isEOFreached = true;
        }
        string newbuff = unprocessedBuffer + string(buff);
        if (newbuff == ""){
            break;
        }
        numSeq_MSA = ReadFastaFromBuffer(newbuff.c_str(),  idList_MSA, annoList_MSA,
                seqList_MSA, seq_type_MSA, unprocessedBuffer, isEOFreached);

        for (i = 0; i < numSeq_MSA; i ++){
            string alignSeq= seqList_MSA[i];
            string topoSeq="";
            string id = idList_MSA[i];
            bool isHavingToposeq = true;
            bool isBadMatch = false;

            unsigned int lengthAlignedSeq = alignSeq.size();

            if (topoMap.find(id) != topoMap.end()){
                topoSeq = topoMap[id];
                isHavingToposeq = true;
            } else {
                isHavingToposeq = false;
                cerr << "topoSeq not found for id " << id << ". Ignore" << endl;
                if (isIgnoreBadseq == true){
                    continue;
                }
            }
            const char* pAlignSeq = alignSeq.c_str(); 
            const char* pTopoSeq = topoSeq.c_str();
            int sizeGaplessSeq = GetGaplessLength(pAlignSeq, lengthAlignedSeq); 
            if (sizeGaplessSeq != int (topoSeq.size()) ) {
                isBadMatch = true;
                cerr << "Length of gaplessAAseq (" << sizeGaplessSeq << 
                    ") and topoSeq (" << topoSeq.size() <<
                    ") for id "<< id << " do not match. Ignore." 
                    << endl;
                if (isIgnoreBadseq == true){
                    continue;
                }
            }else{
                isBadMatch = false;
            }

            if (isHavingToposeq == false or isBadMatch == true){
                if (outFormat == 0) {
                    fprintf(fpout,">%s\n", annoList_MSA[i].c_str());
                    fprintf(fpout, "%s\n", "BADSEQ");
                    //printf("anno: %s\n",annoList_MSA[i].c_str());
                } else {
                    fprintf(fpout,"%-*s :%s\n", maxIDsize, id.c_str(), "BADSEQ");
                }
            }else {
                char *topoAlignSeq = new char [lengthAlignedSeq + 1];

                int jSeq = 0;;
                for (unsigned int ii = 0; ii < lengthAlignedSeq ; ii ++) {
                    if (pAlignSeq[ii] == '-') {
                        topoAlignSeq[ii] = '-';
                    } else {
                        topoAlignSeq[ii] = pTopoSeq[jSeq];
                        jSeq ++;
                    }
                }
                topoAlignSeq[lengthAlignedSeq] = '\0';
                if (outFormat == 0) {
                    fprintf(fpout,">%s\n", annoList_MSA[i].c_str());
                    fprintf(fpout, "%s\n", topoAlignSeq);
                    //printf("anno: %s\n",annoList_MSA[i].c_str());
                } else {
                    fprintf(fpout,"%-*s :%s\n", maxIDsize, id.c_str(), topoAlignSeq);
                }
                delete [] topoAlignSeq;
            }
            cnt ++;
        }
    }
    delete [] buff;
    fclose(fpAln);

    if(fpout != NULL && fpout != stdout) {
        fclose(fpout);
    }

    return 0;
}
/*}}}*/
