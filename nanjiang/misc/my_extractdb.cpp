#include <sys/stat.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>
#include <set>
#include <map>
#include <vector>
#include <sstream>
using namespace std;

typedef unsigned char BYTE;

struct dbindex{/*{{{*/
    int dbfileindex;
    long offset;
    unsigned long size;
};/*}}}*/
string usage="\n\
usage: my_extract [-l LISTFILE] [-o OUTFILE]\n\
                  -dbname STR ID [ID ...]\n\
 -o OUTFILE   Output the result to OUTFILE, (default: stdout)\n\
\n\
Created 2011-10-19, updated 2011-11-02, Nanjiang Shu\n\
";
void PrintHelp() {
    cout << usage << endl;
}
void PrintVerboseHelp() { }

int ReadInIDList(string file, set <string> &idSet)/*{{{*/
{
    ifstream ifp (file.c_str(), ios::binary);
    if (ifp.is_open()){
        ifp.seekg (0, ios::end);
        int length = ifp.tellg();
        ifp.seekg (0, ios::beg);
        char *buffer = new char [ length+1];
        ifp.read(buffer,length);
        buffer[length]='\0'; /*it is important to assign the last byte of a string to 0, otherwise it will be dangerous to use string functions*/
        ifp.close();
        char *pch;
        pch = strtok(buffer, "\n");
        while (pch != NULL){
            if ( strlen(pch)> 0 ){
                idSet.insert(pch);
            }
            pch = strtok(NULL, "\n");
        }
        delete [] buffer;
    }else{
        cerr << "Failed to open idlistfile "
            << file 
            << endl;
        return -1;
    }
    return 0;
}/*}}}*/

bool IsFileExist(const char *name)/*{{{*/
{
    struct stat st;
    return (stat(name,&st) == 0); 
}/*}}}*/ 
string int2string(int number)/*{{{*/
{
    stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
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
int ReadDatabaseIndex_text_old(string &dbname, map<string, dbindex> &dbindexmap,int &maxdbfileindex )/*{{{*/
{
    string indexfile=dbname+string(".index");
    ifstream ifp (indexfile.c_str(), ios::binary);
    string line;
    dbindex tmpdbindex;
    if (ifp.is_open()){
        while (ifp.good()) {
            getline(ifp, line);
            if (line.size() > 0 && line.substr(0,3) != string("DEF")){
                char *id = new char [ line.size()] ; 
                if (sscanf( line.c_str(), "%s %d %ld %ld", id, &(tmpdbindex.dbfileindex)
                        ,&(tmpdbindex.offset),&(tmpdbindex.size)) != 4){
                    cerr << "Read dbindex failed for the database "
                        << dbname << endl;
                    return -1;
                } else {
                    dbindexmap.insert(pair<string, dbindex>(string(id), tmpdbindex));
                    if (tmpdbindex.dbfileindex > maxdbfileindex){
                        maxdbfileindex = tmpdbindex.dbfileindex;
                    }
                }
                delete [] id;
            }
        }
        ifp.close();
        return dbindexmap.size();
    } else {
        cerr << "Failed to open database index file "
            << indexfile 
            << endl;
        return -1;
    }
}/*}}}*/

int ReadDatabaseIndex_text(string &indexfile, map<string, dbindex> &dbindexmap,int &maxdbfileindex )/*{{{*/
{
    /*read in the whole file first*/
    ifstream ifp (indexfile.c_str(), ios::binary);
    dbindex tmpdbindex;
    if (ifp.is_open()){
        // get length of file:
        ifp.seekg (0, ios::end);
        int length = ifp.tellg();
        ifp.seekg (0, ios::beg);
        char *buffer = new char [ length+1];
        ifp.read(buffer,length);
        ifp.close();
        buffer[length] = '\0';
        char *pch;
        pch = strtok(buffer, "\n");
        while (pch != NULL){
            if ( strlen(pch)> 0 && strncmp(pch, "DEF", 3) != 0){
                char *id = new char [ strlen(pch)] ; 
                if (sscanf( pch, "%s %d %ld %ld", id, &(tmpdbindex.dbfileindex)
                        ,&(tmpdbindex.offset),&(tmpdbindex.size)) != 4){
                    cerr << "Failed to read the indexfile "
                        << indexfile << endl;
                    return -1;
                } else {
                    dbindexmap.insert(pair<string, dbindex>(string(id), tmpdbindex));
                    if (tmpdbindex.dbfileindex > maxdbfileindex){
                        maxdbfileindex = tmpdbindex.dbfileindex;
                    }
                }
                delete [] id;
            }
            pch = strtok(NULL, "\n");
        }
        delete [] buffer;
        return dbindexmap.size();
    } else {
        cerr << "Failed to open database index file "
            << indexfile 
            << endl;
        return -1;
    }
}/*}}}*/
int ReadDatabaseIndex_binary(string &indexfile, map<string, dbindex> &dbindexmap,int &maxdbfileindex )/*{{{*/
{
    /*read in the whole file first . binary index*/
    ifstream ifp (indexfile.c_str(), ios::binary);
    if (ifp.is_open()){
        // get length of file:
        ifp.seekg (0, ios::end);
        int length = ifp.tellg();
        ifp.seekg (0, ios::beg);
        BYTE *buffer = new BYTE [ length+1];
        ifp.read(reinterpret_cast<char*>(buffer),length);
        ifp.close();

        BYTE *pBuffer = buffer;
        unsigned int sizedumpedtext=0;
        memcpy(&sizedumpedtext,pBuffer,sizeof(unsigned int));
        //memcpy(dumpedtext,pBuffer,sizeof(unsigned int));
        pBuffer+=sizeof(unsigned int) + sizedumpedtext*sizeof(char);

        unsigned int dumpedidlistsize=0;
        memcpy(&dumpedidlistsize,pBuffer,sizeof(unsigned int));
        pBuffer+=sizeof(unsigned int);
        char *dumpedidlist  = new char[dumpedidlistsize+1];
        memcpy(dumpedidlist,pBuffer,dumpedidlistsize);
        pBuffer+=sizeof(char)*dumpedidlistsize;

        vector <string> idList;
        char *pch = NULL;
        pch = strtok(dumpedidlist, "\n");
        while (pch != NULL){
            if (strlen(pch) > 0){
                idList.push_back(pch);
            }
            pch = strtok(NULL, "\n");
        }
        delete [] dumpedidlist;

        unsigned int numRecord=0;
        memcpy(&numRecord,pBuffer,sizeof(unsigned int));
        pBuffer+=sizeof(unsigned int);
        if (numRecord != idList.size()){
            cerr << "numRecord("
                << numRecord
                << ") != idList.size ("
                << idList.size()
                << ")."
                << endl;
        }

        unsigned char *arrayFileIndex=new unsigned char[numRecord];
        unsigned int  *arrayoffset =new unsigned int[numRecord];
        unsigned int  *arraysize =new unsigned int[numRecord];
        memcpy(arrayFileIndex,pBuffer,sizeof(unsigned char)*numRecord);
        pBuffer+=sizeof(unsigned char)*numRecord;
        memcpy(arrayoffset,pBuffer,sizeof(unsigned int)*numRecord);
        pBuffer+=sizeof(unsigned int)*numRecord;
        memcpy(arraysize,pBuffer,sizeof(unsigned int)*numRecord);
        pBuffer+=sizeof(unsigned int)*numRecord;

        if ((pBuffer-buffer) != length ){
            cerr << "Failed to read binary index file "
                << " filesize = " << length
                << " readsize = " << pBuffer-buffer
                << endl;
        }
        dbindex tmpdbindex;
        for (unsigned int i = 0; i < numRecord; i++) {
            tmpdbindex.dbfileindex=arrayFileIndex[i];
            tmpdbindex.offset=arrayoffset[i];
            tmpdbindex.size=arraysize[i];
            dbindexmap.insert(pair<string, dbindex>(idList[i], tmpdbindex));
            if (tmpdbindex.dbfileindex > maxdbfileindex){
                maxdbfileindex = tmpdbindex.dbfileindex;
            }
        }
        delete [] arrayFileIndex;
        delete [] arraysize;
        delete [] arrayoffset;
        delete [] buffer;
        return dbindexmap.size();
    } else {
        cerr << "Failed to open database index file "
            << indexfile 
            << endl;
        return -1;
    }
}/*}}}*/
int ReadDatabaseIndex(string &dbname, map<string,dbindex>&dbindexmap, int &maxDBIndexNumber)/*{{{*/
{
    string indexfile_binary=dbname+string(".indexbin");
    string indexfile_text=dbname+string(".index");
    if ( IsFileExist( indexfile_binary.c_str()) ) {
        return  ReadDatabaseIndex_binary(indexfile_binary, dbindexmap, maxDBIndexNumber); 
    }else if ( IsFileExist( indexfile_text.c_str()) ) {
        return  ReadDatabaseIndex_text(indexfile_text, dbindexmap, maxDBIndexNumber); 
    }else{
        cerr << "Index file for seqDatabase "
            << dbname
            << "does not exist." << endl;
        return -1;
    }

}/*}}}*/


int GetDBFPList( vector <FILE*> &fpList, string dbname, int maxdbfileindex)/*{{{*/
{
    for (int i=0;i<=maxdbfileindex;i++){
        string filename=dbname + int2string(i) + string(".db");
        FILE *fp = fopen(filename.c_str(),"rb");
        if (fp != NULL){
            fpList.push_back(fp);
        }else{
            cerr << "Failed to open dbfile "
                << filename
                << ". Exit."
                << endl;
            exit(1);
        }
    }
    return 0;
}/*}}}*/
int ExtractData(set <string> &idSet, map <string, dbindex>& dbindexmap, vector <FILE*> &fpList, string &outfile)/*{{{*/
{
    set <string> :: iterator iss;
    int status_fseek;
    ofstream outfs ;
    if (outfile != ""){
        outfs.open(outfile.c_str(), ios::out | ios::binary);
        if (!outfs){
            cerr << "Failed to write to outfile " << outfile << endl;
            cerr << "Reset output to stdout." << endl;
        }
    }
    
    for (iss = idSet.begin(); iss != idSet.end();iss++){
        if(dbindexmap.find(*iss) != dbindexmap.end()){
            dbindex *pDBindex = &dbindexmap[*iss];
            if ((status_fseek = fseek(fpList[pDBindex->dbfileindex],
                            pDBindex->offset, SEEK_SET)) == 0) {
                size_t nread;
                char *buffer = new char [ pDBindex->size+1];
                buffer[pDBindex->size]='\0'; /*it is important to assign the last byte of a string to 0, otherwise it will be dangerous to use string functions*/
                nread = fread(buffer, sizeof(char), pDBindex->size,
                        fpList[pDBindex->dbfileindex]);
                if (nread != pDBindex->size){
                    cerr <<"fread failed" << endl;
                }
                if (outfs){
                    outfs.write(buffer, nread);
                }else{
                    cout.write(buffer, nread);
                }
                delete [] buffer;
            } else {
                cerr <<"Fatal! fseek failed for id \'"
                    << *iss 
                    << "\' in the database. Exit."
                    << endl;
                exit(1);
        }

        } else {
            cerr << "ID \'"
                << *iss 
                << "\' not found in the database. Ignore."
                << endl;
            continue;
        }
    }
    if (outfs){
        outfs.close();
    }
    return 0;
}/*}}}*/

int main(int argc, char** argv)/*{{{*/
{
    bool isNonOptionArg = false;
    if(argc < 2) {
        PrintHelp();
        return 1;
    }
    int i,j;
    double value = 0.0;
    string control_option = "q"; //options which control the program, and does not take parameters
    bool isQuiet = false;
    string outfile = "";
    string dbname="";
    string idlistfile="";
    set <string> idSet;

    i = 1;
    while(i < argc)/*{{{*/ {
        if(argv[i][0] == '-' && !isNonOptionArg){
            isNonOptionArg = false;
            if( control_option.find(argv[i][1]) != string::npos)/*if argv[i][1] is in control_option, it might be used as -aqs*/ {
                for(j = 1 ; j < int(strlen(argv[i])); j++) {
                    switch (argv[i][j]) {
                        case 'q': isQuiet = true; break;
                        default : cerr << 
                                  "Invalid option, non-control option \'" 
                                  << argv[i][j] 
                                  << "\' can not be used together with contorl-option" 
                                  << endl; 
                                  return -1;
                    }
                }
                i ++;
            } else if( string(argv[i])=="-h" || string(argv[i])=="--help") {
                PrintHelp(); 
                return 1;
            } else if(string(argv[i])=="-H") {
                PrintVerboseHelp();
                return 0;
            } else if(string(argv[i])== "-o" ||
                    string(argv[i])=="--outfile")  {
                if( ( i = option_parser_filename(argc, argv, i, outfile)) == -1)
                    return -1;
            } else if(string(argv[i])== "-dbname" ||
                    string(argv[i])=="--dbname")  {
                if( ( i = option_parser_filename(argc, argv, i, dbname)) == -1)
                    return -1;
            } else if(string(argv[i])== "-l" ||
                    string(argv[i])=="--l")  {
                if( ( i = option_parser_filename(argc, argv, i, idlistfile)) == -1)
                    return -1;
            } else if( string(argv[i])== "--value")  {
                if( ( i = option_parser_numeric(argc, argv, i, value, true,
                                0.0, 15.0)) == -1)
                    return -1;
            } else if (string(argv[i]) == "--")/*next item is non option argument*/ {
                isNonOptionArg = true;
                i ++;
                continue;
            } else {
                cerr << "Error! Invalid argument \'"<< argv[i] << "\'" <<
                    endl;
                return -1;
            }
        } else /*positional argument*/ {
            idSet.insert(argv[i]);
            i ++;
        }
    }/*}}}*/
    if (idlistfile != ""){
        ReadInIDList(idlistfile, idSet);
    }

    if (idSet.size() <= 0 ){
        cerr << "No ID set. Exit." << endl;
        return 1;
    }

    map <string, dbindex> dbindexmap;
    int maxDBIndexNumber = 0;

    if (dbname != ""){
#ifdef DEBUG_TIME_READINDEX
        clock_t start,finish;
        double duration;
        start=clock();
#endif
        if (ReadDatabaseIndex(dbname, dbindexmap, maxDBIndexNumber) == -1 ){
            cerr << "Read Database Index failed for "<<dbname
                << ". Exit."
                << endl;
            return -1;
        }
#ifdef DEBUG_TIME_READINDEX
        finish=clock();
        duration = double(finish-start)  /double(CLOCKS_PER_SEC);
        cout << "ReadDatabaseIndex2 costs " << duration << " seconds." << endl;
#endif
    }else{
        cerr << "dbname not set. Exit." << endl;
        return -1;
    }
    vector <FILE*> fpList;
    GetDBFPList( fpList, dbname, maxDBIndexNumber);
    //cout << "maxDBIndexNumber = " << maxDBIndexNumber << endl;

    ExtractData(idSet, dbindexmap, fpList, outfile);

    for(i=0;i <= maxDBIndexNumber; i++){
        fclose(fpList[i]);
    }
    return 0;
} /*}}}*/

