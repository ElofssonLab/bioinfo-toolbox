#!/usr/bin/env python
# Description:
# convert topcons_single allinfo to fasta file
import os
import sys
import myfunc
import libtopologycmp as lcmp

# ChangeLog 2013-02-01
#   extraction method 1 is added, in this case, only those entries predicted
#   with the same topology for all sub-predictors are exported
# ChangeLog 2013-02-21
#   use MatchList, count the number of agreed agreed sub-predictions to
#   topcons_single.
#
min_TM_overlap = 5

progname = os.path.basename(sys.argv[0])
usage = """
Usage:  %s allinfo_file [-outpath DIR]

Description: Extract topology sequences from Topcons_single result file

OPTIONS:
  -outpath DIR    Set output path, (default: the same as input file)
  -m INT          Set extracting method, (default: 0)
                  method 0: all predicted topologies are output
                  method 1: only those predicted the same (same numTM and same
                            NtermStatus) for all sub-predictors are extracted.
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2012-03-27, updated 2013-02-21, Nanjiang Shu
"""%(progname)

BLOCK_SIZE=100000

def PrintHelp():
    print usage
def ExtractFromTopconsSingleAllinfo(recordContent):#{{{
    record = {}
    tmplines = recordContent.split("\n")
    lines = []
    for line in tmplines:
        if line:
            lines.append(line)
    numLine = len(lines)
    if numLine < 1:
        return {}
    i = 0
    anno = ""
    predtopo_topcons_single = ""
    predtopo_scampi_single = ""
    predtopo_hmmtop = ""
    predtopo_stmhmm = ""
    predtopo_memsat = ""
    rlty = -100.0
    isTM = True
#   print "recordContent=\n", recordContent
    while i < numLine:
#         print i
        if lines[i][0:4] == 'NAME':
            anno += lines[i][9:]
            i += 1
        elif lines[i][0:4] == 'ISTM':
            strs = lines[i].split()
            if strs[1] != 'YES':
                isTM = False
                break
            else:
                isTM = True
            i += 1
        elif lines[i][0:4] == 'MTHD':
            method = lines[i].split()[1]
            j = 1
            tmp_topo = ""
            while i+j < numLine and lines[i+j][0:4] == 'TOPO':
                ss = lines[i+j].split()
                if len(ss) > 1:
                    tt = ss[1]
                else:
                    tt = ''
                tmp_topo += tt
                j+=1
            i += j
            if method == "topcons-single":
                predtopo_topcons_single = tmp_topo
            elif method == "scampi-single":
                predtopo_scampi_single = tmp_topo
            elif method == "hmmtop":
                predtopo_hmmtop = tmp_topo
            elif method == "stmhmm":
                predtopo_stmhmm = tmp_topo
            elif method == "memsat":
                predtopo_memsat = tmp_topo
        elif lines[i][0:4] == 'RLTY':
            rlty = float(lines[i].split()[1])
            i += 1
        else:
            i += 1
    #print "anno", anno
    if isTM and predtopo_topcons_single == '':
#         print "anno=",anno
#         print "recordContent\"%s\""%recordContent
        print >> sys.stderr, "%s isTM=yes but no predtopo" %(
                anno.split()[0])
    if isTM:
        record['anno'] = anno
        record['predtopo_topcons_single'] = predtopo_topcons_single
        record['predtopo_scampi_single'] = predtopo_scampi_single
        record['predtopo_hmmtop'] = predtopo_hmmtop
        record['predtopo_stmhmm'] = predtopo_stmhmm
        record['predtopo_memsat'] = predtopo_memsat
        record['rlty'] = rlty
        record['seqlength'] = len(predtopo_topcons_single)
#         print record
        return record
    else:
        return {}
#}}}
def Read_topconssingle_result_from_buffer(buff, recordList, isEOFreached):#{{{
    if not buff:
        return ""
    unprocessedBuffer=""
    beg=0
    end=0
    while 1:
        beg=buff.find("NAME",beg)
        if beg >= 0:
            end=buff.find("\n\n",beg+1)
            if end >=0:
                recordContent = buff[beg:end]
                record = ExtractFromTopconsSingleAllinfo(recordContent)
                if record != {}:
                    recordList.append(record)
                beg = end
            else:
                unprocessedBuffer = buff[beg:]
                break
        else:
            unprocessedBuffer = buff[end:]
            break
    if isEOFreached and unprocessedBuffer:
        recordContent = unprocessedBuffer
        record = ExtractFromTopconsSingleAllinfo(recordContent)
        if record != {}:
            recordList.append(record)
        unprocessedBuffer = ""
    return unprocessedBuffer
    #}}}
def TopconsSingle2Fasta(infile, outpath):#{{{
    try:
        rootname=os.path.basename(os.path.splitext(infile)[0])
        outfile_topcons_single = outpath + os.sep + rootname + "_topcons_single.topo"
        outfile_scampi_single = outpath + os.sep + rootname + "_scampi_single.topo"
        outfile_hmmtop = outpath + os.sep + rootname + "_hmmtop.topo"
        outfile_stmhmm = outpath + os.sep + rootname + "_stmhmm.topo"
        outfile_memsat = outpath + os.sep + rootname + "_memsat.topo"
        outRLTYFile = outpath + os.sep + rootname + "_topcons_single.rlty"
        fpout_topcons_single = open(outfile_topcons_single, "w")
        fpout_scampi_single = open(outfile_scampi_single, "w")
        fpout_hmmtop = open(outfile_hmmtop, "w")
        fpout_stmhmm = open(outfile_stmhmm, "w")
        fpout_memsat = open(outfile_memsat, "w")
        fpout_rlty = open(outRLTYFile, "w")
        fpin = open(infile, "r")
        unprocessedBuffer=""
        isEOFreached = False
        processedTopoIDSet = set([])
        while 1:
            buff = fpin.read(BLOCK_SIZE)
            if len(buff) < BLOCK_SIZE:
                isEOFreached=True
            buff = unprocessedBuffer + buff
            recordList = []
            unprocessedBuffer = Read_topconssingle_result_from_buffer(
                    buff, recordList, isEOFreached)
            if len(recordList) > 0: 
                for record in recordList:
                    seqid = myfunc.GetSeqIDFromAnnotation(record['anno'])
                    if record['predtopo_topcons_single'] != "":
                        fpout_topcons_single.write(">%s topcons_single rlty=%.2f\n"%(
                            record['anno'], record['rlty']))
                        fpout_topcons_single.write("%s\n"%record['predtopo_topcons_single'])
                    if record['predtopo_scampi_single'] != "":
                        fpout_scampi_single.write(">%s scampi_single\n"%(record['anno'] ))
                        fpout_scampi_single.write("%s\n"%record['predtopo_scampi_single'])
                    if record['predtopo_hmmtop'] != "":
                        fpout_hmmtop.write(">%s hmmtop\n"%(record['anno'] ))
                        fpout_hmmtop.write("%s\n"%record['predtopo_hmmtop'])
                    if record['predtopo_stmhmm'] != "":
                        fpout_stmhmm.write(">%s stmhmm\n"%(record['anno'] ))
                        fpout_stmhmm.write("%s\n"%record['predtopo_stmhmm'])
                    if record['predtopo_memsat'] != "":
                        fpout_memsat.write(">%s memsat\n"%(record['anno'] ))
                        fpout_memsat.write("%s\n"%record['predtopo_memsat'])
                    if record['rlty'] != -100.0:
                        fpout_rlty.write("%s %.2f\n"%(seqid, record['rlty']))
            if isEOFreached == True:
                break
        fpin.close()
        fpout_topcons_single.close()
        fpout_scampi_single.close()
        fpout_memsat.close()
        fpout_stmhmm.close()
        fpout_hmmtop.close()
        fpout_rlty.close()

        print "Result have been output to"
        print "\t%s"%outfile_topcons_single
        print "\t%s"%outfile_scampi_single
        print "\t%s"%outfile_hmmtop
        print "\t%s"%outfile_stmhmm
        print "\t%s"%outfile_memsat
        print "\t%s"%outRLTYFile

    except IOError:
        print >> sys.stderr, "Failed to open file %s for read"%(infile)
        raise
#}}}
def IsAllIdenticalTopology(topoList):#{{{
    numSeq = len(topoList)
    if numSeq <= 1:
        return True
    else:
        posTMList = [myfunc.GetTMPosition(topo) for topo in topoList]
        NtermStateList = [ lcmp.GetNtermState(topo) for topo in topoList]
        numTMList = [len(posTM) for posTM in posTMList]
        min_TM_overlap = 5

        for i in xrange(numSeq-1):
            for j in xrange(i+1, numSeq):
                if not lcmp.IsIdenticalTopology(NtermStateList[i],
                        NtermStateList[j], numTMList[i], numTMList[j],
                        posTMList[i], posTMList[j], topoList[i], topoList[j],
                        min_TM_overlap):
                    return False
        return True
#}}}

def TopconsSingle2Fasta_method1(infile, outpath):#{{{
    try:
        rootname = os.path.basename(os.path.splitext(infile)[0])
        outfile_topcons_single = outpath + os.sep + rootname + "_topcons_single.m1.topo"
        outfile_agreement = outpath + os.sep + rootname + ".agreement.stat.txt"
        logfile1 = outpath + os.sep + rootname + ".m1.idt.log"
        logfile2 = outpath + os.sep + rootname + ".m1.nonidt.log"
        fpout_topcons_single = open(outfile_topcons_single, "w")
        fplog1 = open(logfile1, "w")
        fplog2 = open(logfile2, "w")
        fpout_agree = open(outfile_agreement, "w")

        fpin = open(infile, "r")
        unprocessedBuffer=""
        isEOFreached = False
        processedTopoIDSet = set([])
        while 1:
            buff = fpin.read(BLOCK_SIZE)
            if len(buff) < BLOCK_SIZE:
                isEOFreached=True
            buff = unprocessedBuffer + buff
            recordList = []
            unprocessedBuffer = Read_topconssingle_result_from_buffer(
                    buff, recordList, isEOFreached)
            if len(recordList) > 0: 
                for record in recordList:
# ignore cases where topcons_single does not give a prediction
                    if record['predtopo_topcons_single'] == "":
                        continue
                    seqid = myfunc.GetSeqIDFromAnnotation(record['anno'])
                    topoList = []
                    topoList.append(record['predtopo_scampi_single'])
                    topoList.append(record['predtopo_hmmtop'])
                    topoList.append(record['predtopo_stmhmm'])
                    topoList.append(record['predtopo_memsat'])
                    (matchList, numIDTtopo, numPredictor) = lcmp.MatchTopology(
                            record['predtopo_topcons_single'], topoList,
                            min_TM_overlap, seqid)
                    fpout_agree.write("%s\t%d\t%d"%(seqid, numIDTtopo,
                        numPredictor))
                    for tt in matchList:
                        fpout_agree.write("\t%d"%(tt))
                    fpout_agree.write("\t%6.2f"%(record['rlty']))
                    fpout_agree.write("\n")

                    if numIDTtopo == numPredictor and numIDTtopo >= 2:
                        msg = ">%s topcons_single RLTY=%.2f"\
                                "numIDTtopo=%d numPredictor=%d\n"
                        fpout_topcons_single.write(msg%(record['anno'] ,
                            record['rlty'], numIDTtopo, numPredictor))
                        fpout_topcons_single.write("%s\n"%record['predtopo_topcons_single'])
                        msg = "%s RLTY= %.2f numIDTtopo= %d numPredictor= %d"
                        print >> fplog1, msg%(seqid, record['rlty'],
                            numIDTtopo, numPredictor)
                    else:
                        msg = "%s RLTY= %.2f numIDTtopo= %d numPredictor= %d"
                        print >> fplog2, msg%(seqid, record['rlty'],
                            numIDTtopo, numPredictor)

            if isEOFreached == True:
                break
        fpin.close()
        fpout_topcons_single.close()
        fplog1.close()
        fplog2.close()
        fpout_agree.close()

        print "Result have been output to"
        print "\t%s"%outfile_topcons_single
        print "\t%s"%outfile_agreement
        print "\t%s"%logfile1
        print "\t%s"%logfile2

    except IOError:
        msg = "Failed to read file {} in function {}"
        print >> sys.stderr, msg.format(infile, sys._getframe().f_code.co_name)
        return 1
#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = ""
    method  = 0
    allinfoFile = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            allinfoFile = argv[i]
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-outpath", "--outpath"]:
                outpath = argv[i+1]
                i += 2
            elif argv[i] in ["-m", "--m", "-method", "--method"]:
                try:
                    method = int(argv[i+1])
                    i += 2
                except (TypeError, ValueError, IndexError):
                    print >> sys.stderr, "argument error: -m INT"
                    return 1
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            allinfoFile = argv[i]
            i += 1

    if allinfoFile == "":
        print >> sys.stderr, "Infile not set. Exit."
        return 1
    elif not os.path.exists(allinfoFile):
        print >> sys.stderr, "allinfoFile %s is empty.  Exit"%(
                allinfoFile)
        return 1

    if outpath == "":
        outpath = os.path.dirname(allinfoFile)
        if outpath == "":
            outpath = "."
    elif not os.path.exists(outpath):
        os.system("mkdir -p %s"%(outpath))

    if method == 0:
        TopconsSingle2Fasta(allinfoFile, outpath)
    elif method == 1:
        TopconsSingle2Fasta_method1(allinfoFile, outpath)
    else:
        print >> sys.stderr, "Wrong method, must be 0 or 1, but %d was set." \
                % (method)
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))
