#!/usr/bin/env python
# Description:
import os
import sys
import tempfile
import myfunc
import subprocess

#ChangeLog 2014-09-11
#   Added the option -ext
#ChangeLog 2014-11-07
#   changed reading of sequence database, to solve the RAM overflow problem
#ChangeLog 2014-11-24 
#   added the option -extra-description-file
#   added the option -bigmem, in case of big memory, load the whole sequence
#   database into memory
#ChangeLog 2014-12-15 
#   add option -nr 90, using CD-hit to cut the sequence identity level
#   when dbname is not set, only output the splitted file

progname =  os.path.basename(sys.argv[0])
usage = """
usage: %s -mapfile fam2seqidmapfile -dbname dbname -seqdb seqdb

OPTIONS:
  -q              Quiet mode
  -extra-description-file FILE
                  add extra description file for each seqid, to be added to the
                  new sequence database format, one line per record, seqid
                  description
  -splitdir [-tmpdir] DIR 
                  Dir to store splitted fasta files
  -ext    STR     Extension of the out file, (default: .fa)
  -gzip           Gzip the individual fasta files
  -nr     INT     Cut sequence identity to <INT%%
  -h, --help      Print this help message and exit

Created 2012-06-08, updated 2014-12-15, Nanjiang Shu 

Examples:
    %s -tmpdir splittedfasta -mapfile test.pfamid2seqid -dbname famseqdb.nr90 -nr 90 -seqdb uniprot.fasta
"""%(progname, progname)

rundir = os.path.dirname(sys.argv[0])
binpath = rundir

def PrintHelp():
    print usage

def ReadSeqDBDict(infile):#{{{
    seqdbDict = {}
    (idList, annotationList, seqList) = myfunc.ReadFasta(infile)
    for i in xrange(len(idList)):
        seqdbDict[idList[i]] = (annotationList[i], seqList[i])
    return seqdbDict
#}}}
def OutputPfamFastaFile(seqidList, pfamid, seqdbDict, hdl_seqdb, extra_desp_dict, outpath):#{{{
    outfile = "%s%s%s%s"%(outpath, os.sep, pfamid, g_params['out_ext'])
    fpout = myfunc.myopen(outfile, None, "w", True)
    isAddExtraDescription = False
    if len(extra_desp_dict) > 0:
        isAddExtraDescription = True

    for seqid in seqidList:
        if seqid.find("UniRef") != -1:
            try: 
                ss = seqid.split("_")
                seqid = ss[1]
            except IndexError:
                pass

        if g_params['isBigmem']:
            try:
                record = seqdbDict[seqid]
                (tmpanno, tmpseq) = record
                if isAddExtraDescription:
                    try:
                        extraanno = extra_desp_dict[seqid]
                    except KeyError:
                        extraanno = ""
                    if extraanno != "":
                        tmpanno = "%s %s"%(extraanno, tmpanno)
                fpout.write(">%s\n%s\n"%(tmpanno, tmpseq))
            except KeyError:
                print >> sys.stderr, "seqid %s not found in seqdb"%(seqid)
        else:
            record = hdl_seqdb.GetRecord(seqid)
            if record:
                if isAddExtraDescription:
                    try:
                        extraanno = extra_desp_dict[seqid]
                    except KeyError:
                        extraanno = ""
                    if extraanno == "":
                        fpout.write("%s"%(record))
                    else:
                        (tmpseqid, tmpanno, tmpseq) = myfunc.ExtractFromSeqWithAnno(record)
                        tmpanno = "%s %s"%(extraanno, tmpanno)
                        fpout.write(">%s\n%s\n"%(tmpanno, tmpseq))
                else:
                    fpout.write("%s"%(record))
            else:
                print >> sys.stderr, "seqid %s not found in seqdb"%(seqid)
    myfunc.myclose(fpout)
    if g_params['isGzip']:
        cmd = ["gzip", "-N", "-f", outfile]
        print " ".join(cmd)
        subprocess.check_call(cmd, stdout=open(os.devnull,"w"))


#}}}
def GetCDHitWordSize(nrlevel):#{{{
    if nrlevel>100:
        print >> sys.stderr, "%s: Error. nrlevel (%d) > 100. exit"%(sys.argv[0],nrlevel)
        sys.exit(1)
    elif nrlevel >= 70:
        return 5
    elif nrlevel >= 60:
        return 4
    elif nrlevel >= 50:
        return 3
    elif nrlevel >= 40:
        return 2
    else:
        print >> sys.stderr, "%s: Error. nrlevel (%d) < 40. exit"%(sys.argv[0],nrlevel)
        sys.exit(1)
#}}}
def CutSeqIDT(pfamid, outpath, nrlevel, wordsize): #{{{
    fafile = "%s%s%s%s"%(outpath, os.sep, pfamid, g_params['out_ext'])
    nr_fafile = "%s%s%s%s%s"%(outpath, os.sep, pfamid, ".nr%d"%(nrlevel), g_params['out_ext'])
    tmpoutfile = "%s%s%s%s%s.tmp"%(outpath, os.sep, pfamid, ".nr%d"%(nrlevel), g_params['out_ext'])
    cmd =  ["cd-hit", "-T", "8", "-i", fafile, "-o" ,tmpoutfile, "-c",
            "%.f"%(nrlevel/100.0),"-n", "%d"%(wordsize), "-M", "8000"]
    try:
        print " ".join(cmd)
        subprocess.check_call(cmd, stdout=open(os.devnull,"w"))
        if os.path.exists(tmpoutfile):
            cmd = ["/bin/mv","-f",tmpoutfile, nr_fafile]
            subprocess.check_call(cmd)
            cmd = ["/bin/rm", "-f", fafile, "%s.clstr"%(tmpoutfile), "%s.bak.clstr"%(tmpoutfile)]
            subprocess.check_call(cmd)

            if g_params['isGzip']:
                cmd = ["gzip", "-N", "-f", nr_fafile]
                print " ".join(cmd)
                subprocess.check_call(cmd, stdout=open(os.devnull,"w"))

    except subprocess.CalledProcessError, e:
        print e
#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = ""
    mapfile = ""
    dbname =  ""
    seqdb = ""
    tmpdir = ""
    out_ext = ".fa"
    extra_description_file = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-dbname", "--dbname"]:
                dbname = argv[i+1]
                i += 2
            elif argv[i] in ["-seqdb", "--seqdb"]:
                seqdb = argv[i+1]
                i += 2
            elif argv[i] in ["-mapfile", "--mapfile"]:
                mapfile = argv[i+1]
                i += 2
            elif argv[i] in ["-extra-description-file", "--extra-description-file"]:
                extra_description_file = argv[i+1]
                i += 2
            elif argv[i] in ["-ext", "--ext"]:
                g_params['out_ext'] = argv[i+1]
                i += 2
            elif argv[i] in ["-splitdir","--splitdir", "-tmpdir", "--tmpdir"]:
                tmpdir = argv[i+1]
                i += 2
            elif argv[i] in ["-l", "--l"] :
                idListFile = argv[i+1]
                i += 2
            elif argv[i] in ["-nr", "--nr"] :
                g_params['nrlevel'] = int(argv[i+1])
                i += 2
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            elif argv[i] in ["-gzip", "--gzip"]:
                g_params['isGzip'] = True
                i += 1
            elif argv[i] in ["-bigmem", "--bigmem"]:
                g_params['isBigmem'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            print >> sys.stderr, "Error! Wrong argument:", argv[i]
            return 1
    if dbname == "" and tmpdir == "":
        print >> sys.stderr, "Error! neither dbname nor splitdir are set."
        return 1

    if seqdb == "":
        print >> sys.stderr, "Error! seqdb not set. Exit."
        return 1
    if mapfile == "":
        print >> sys.stderr, "Error! mapfile not set. Exit."
        return 1

    if tmpdir == "":
        tmpdir = tempfile.mkdtemp()
    print "tmpdir=\"%s\""%(tmpdir)

    if not os.path.exists(tmpdir):
        try:
            subprocess.check_call(["mkdir", "-p", tmpdir])
        except subprocess.CalledProcessError, e:
            return 1

#     famid2seqidDict = myfunc.ReadFam2SeqidMap(mapfile)

    g_params['cdhit_wordsize'] = GetCDHitWordSize(g_params['nrlevel'])

    seqdbDict = None
    hdl_seqdb = None

    if g_params['isBigmem']:
        seqdbDict = ReadSeqDBDict(seqdb)
    else:
        hdl_seqdb = myfunc.MyDB(seqdb)
        if hdl_seqdb.failure:
            print >> sys.stderr, "Failed to load seqdb %s. exit"%(seqdb)
            return 1

    pfamidList  = []

    extra_desp_dict = {}
    if extra_description_file != "":
        hdl_extra = myfunc.ReadLineByBlock(extra_description_file)
        if hdl_extra.failure:
            print >> sys.stderr, "Failed to read extra_description_file %s."%(extra_description_file)
            return 1
        lines = hdl_extra.readlines()
        while lines != None:
            for line in lines:
                line = line.strip()
                if not line or line[0] == "#":
                    continue
                seqid = myfunc.GetSeqIDFromAnnotation(line)
                if seqid != "":
                    extra_desp_dict[seqid] = line
            lines = hdl_extra.readlines()


    hdl = myfunc.ReadLineByBlock(mapfile)
    if hdl.failure:
        print >> sys.stderr, "Failed to read mapfile %s. exit"%(mapfile)
        return 1

    cntfam = 0
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            line = line.strip()
            if not line or line[0] == "#":
                continue
            strs = line.split()
            if len(strs) > 2:
                cntfam += 1
                pfamid = strs[0]
                pfamidList.append(pfamid)
                seqidList = strs[2:]
                OutputPfamFastaFile(seqidList, pfamid, seqdbDict, hdl_seqdb, extra_desp_dict, tmpdir)

                if g_params['nrlevel'] < 100:
                    CutSeqIDT(pfamid, tmpdir, g_params['nrlevel'],
                            g_params['cdhit_wordsize'])
                if not g_params['isQuiet']:
                    print "%d\t%s output"%(cntfam, pfamid)
            else:
                msg="broken item in file %s: line \"%s\""
                print >> sys.stderr, msg%(mapfile, line)
        lines = hdl.readlines()
    hdl.close()

    tmpfamidlistfile = tmpdir + os.sep + "tmpfamidlist.famidlist"
    ftmp = open(tmpfamidlistfile, "w")
    for famid in pfamidList:
        print >> ftmp, famid
    ftmp.close()

    if dbname != "":
        newext = g_params['out_ext']
        if g_params['nrlevel'] != 100:
            newext = ".nr%d%s"%(g_params['nrlevel'], g_params['out_ext'])
        cmd = "python %s/my_formatdb.py -l %s -dbname %s -datapath %s -dataext %s" % (
                binpath, tmpfamidlistfile, dbname, tmpdir, newext)
        os.system(cmd)

    if hdl_seqdb:
        hdl_seqdb.close()

    #os.system("rm -rf %s"%(tmpdir))
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = False
    g_params['out_ext'] = ".fa"
    g_params['isBigmem'] = False
    g_params['nrlevel'] = 100
    g_params['cdhit_wordsize'] = 5
    g_params['isGzip'] = False
    return g_params
#}}}
if __name__ == '__main__' :

    g_params = InitGlobalParameter()
    sys.exit(main(g_params))
