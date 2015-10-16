#!/usr/bin/python
# Description:
import os
import sys
import myfunc
import libtopologycmp as lcmp
#import bioinformatics
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s hhalignresultfile [-seqdb STR] [-seqfile STR]
"""%(progname)

usage_ext="""
Description:

OPTIONS:
  -o  OUTFILE      output the pairaln file
  -ot OUTFILE      output table info file
  -os OUTFILE      output alignment statistics file
  -evalue   FLOAT  evalue threshold, (default: 1e-3)
  -coverage FLOAT  coverage of the aligned region to the shorter sequence, (default: 0.0)
  -seqdb STR       Set sequence database
  -l LISTFILE      Set the listfile
  -q               Quiet mode
  -h, --help       Print this help message and exit

Created 2013-05-17, updated 2013-05-22, Nanjiang Shu 
"""
usage_exp="""
Examples:
"""

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}


def GetUnalignedHeadTail(seq, pos_aln_begin, pos_aln_end):
    return (seq[0:pos_aln_begin], seq[pos_aln_end:])
def FillUnalignedGapForHead(head1, head2):
    len1 = len(head1)
    len2 = len(head2)
    if len1 < len2:
        head1 = "-"*(len2-len1) + head1
    elif len1 > len2:
        head2 = "-"*(len1-len2) + head2
    return (head1, head2)

def FillUnalignedGapForTail(tail1, tail2):
    len1 = len(tail1)
    len2 = len(tail2)
    if len1 < len2:
        tail1 = tail1 + "-"*(len2-len1) 
    elif len1 > len2:
        tail2 = tail2 +  "-"*(len1-len2) 
    return (tail1, tail2)

def ScanHHHitLine(line):#{{{
# No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
#  1 tr|A6GYR2|A6GYR2_FLAPJ Rhomboi  99.9 2.6E-31 2.6E-31  217.3   0.0  187   39-280     1-246 (246)
    try: 
# changeLog 2013-05-22 
# the parser is now working with more than 10000 (5 digits) columns
        datastr = line[35:]
        strs = datastr.split()

        prob = float(strs[0])
        evalue = float(strs[1])
        pvalue = float(strs[2])
        score = float(strs[3])
        num_align_col = int(strs[5])
        strs1 = strs[6].split("-")
        pos_query_begin = int(strs1[0]) - 1 #the position in hhalign is starting from zero
        pos_query_end = int(strs1[1]) # the end residue is encluded in hhalign, but in this script, residue of the end position is not included, e.g. (0,2) is not residue 0 and 1
        if len(strs) == 8:
            strs1 = strs[7].replace("-", " ").replace("(", " ").split()
            pos_template_begin = int(strs1[0]) - 1 
            pos_template_end = int(strs1[1])
            template_length = int(strs1[2].rstrip(')'))
        elif len(strs) == 9:
            strs1 = strs[7].split("-")
            pos_template_begin = int(strs1[0]) - 1 
            pos_template_end = int(strs1[1])
            template_length = int(strs[8].lstrip('(').rstrip(')'))
#         prob = float(line[35:40])
#         evalue = float(line[41:48])
#         pvalue = float(line[49:56])
#         score = float(line[57:63])
#         num_align_col = int(line[70:74])
#         strs = line[75:84].split("-")
#         pos_query_begin = int(strs[0]) - 1 #the position in hhalign is starting from zero
#         pos_query_end = int(strs[1]) # the end residue is encluded in hhalign, but in this script, residue of the end position is not included, e.g. (0,2) is not residue 0 and 1
#         strs = line[85:94].split("-")
#         pos_template_begin = int(strs[0]) - 1 
#         pos_template_end = int(strs[1])
#         template_length = int(line[95:].rstrip(')'))

        return (prob, evalue, pvalue, score, num_align_col, pos_query_begin,
                pos_query_end, pos_template_begin, pos_template_end,
                template_length)
    except (IndexError, ValueError, TypeError):
        print "Bad HHHitLine:\"%s\""%line
        raise
    
#}}}
def ScanHHAlignStatLine(line):#{{{
    try:
        strs = line.split()
        if len(strs) == 7:
            prob    = float(strs[0].split("=")[1])
            evalue  = float(strs[1].split("=")[1])
            score   = float(strs[2].split("=")[1])
            alignedcol = int(strs[3].split("=")[1])
            identity   = float(strs[4].split("=")[1].rstrip('%'))
            similarity = float(strs[5].split("=")[1])
            sum_prob   = float(strs[6].split("=")[1])
        return (prob, evalue, score, alignedcol, identity, similarity, sum_prob)
    except (IndexError, ValueError, TypeError):
        print >> sys.stderr, "Bad line: \"%s\""%(line)
        raise
#}}}
def ReadHHAlignResult(infile):
    try:
        fpin = open(infile,"r")
        lines = fpin.read().split("\n")
        fpin.close()
        hitList = []
        numLine = len(lines)
        cnt = 0
        while cnt < numLine:
            if cnt < numLine and lines[cnt][0:5] == "Query":
                query_description = lines[cnt][14:]
                query_seqid = myfunc.GetSeqIDFromAnnotation(query_description)
                cnt += 1
            elif cnt < numLine and lines[cnt].find("Match_columns") == 0:
                try:
                    query_legnth = int(lines[cnt].split()[1])
                except (IndexError, ValueError, TypeError):
                    print >> sys.stderr, "Bad line:\"%s\""%(lines[cnt])
                    raise
                cnt += 1
            elif lines[cnt][0:7] == " No Hit":
                cnt += 1
                while cnt < numLine and lines[cnt] != "":
                    hit = {}
                    #print lines[cnt]
                    (hit['prob'], hit['evalue'], hit['pvalue'], hit['score'],
                            hit['num_align_col'], hit['pos_query_begin'],
                            hit['pos_query_end'], hit['pos_template_begin'],
                            hit['pos_template_end'], hit['template_length']
                            ) = ScanHHHitLine(lines[cnt])
                    hit['query_length'] = query_legnth
                    hit['query_seqid'] = query_seqid
                    hit['query_description'] = query_description
                    hitList.append(hit)
                    cnt += 1
                cnt += 1
            elif cnt < numLine and lines[cnt][0:3] == "No ":
                j = 0
                jstart = cnt
                try:
                    hitIndex = int(lines[cnt].split()[1])-1
                except (IndexError, ValueError, TypeError):
                    print >> sys.stderr, "Bad hit line: %s"%lines[cnt]
                    raise
                try:
                    hit = hitList[hitIndex]
                except IndexError:
                    print >> sys.stderr, "hitIndex=%d, numHit=%d"%(hitIndex, len(hitList))
                    raise
                j += 1
                hit['hit_description'] = lines[jstart+j].lstrip(">")
                hit['hit_seqid'] = myfunc.GetSeqIDFromAnnotation(hit['hit_description'])
                j += 1
                hit['statline'] = lines[jstart+j]
                (hit['prob'], hit['evalue'], hit['score'], hit['alignedcol'],
                        hit['identity'], hit['similarity'], hit['sum_prob']
                        ) = ScanHHAlignStatLine(hit['statline'])
                j += 2
                alnseqList1 = []
                alnseqList2 = []
                while lines[jstart+j][0:2] == "Q ":
                    try:
                        l1 = lines[jstart+j][17:]
                        l2 = lines[jstart+j+4][17:]
                        s1 = l1.split()[1]
                        s2 = l2.split()[1]
                        alnseqList1.append(s1)
                        alnseqList2.append(s2)
                        if lines[jstart+j+5].strip() != "":
                            j += 8
                        else:
                            j += 7
                    except IndexError:
                        msg = "Bad hhr result file %s. l1=%s, l2=%s"
                        print >> sys.stderr, msg%(infile, l1, l2)
                        sys.exit(1)
                hit['query_alignseq'] = "".join(alnseqList1)
                hit['template_alignseq'] = "".join(alnseqList2)
                cnt += j
            elif lines[cnt][0:4] == "Done":
                break
            else:
                cnt += 1
        return hitList
    except IOError:
        print >> sys.stderr, "Failed to read infile %s"%(infile)
        return []

def HHAlign2Pairaln(infile, evalue_threshold, coverage_threshold, hdl_seq, #{{{
        fpout, fpout_tableinfo, fpout_stat) :
    if not os.path.exists(infile):
        print >> sys.stderr, "infile %s does not exist, Ignore"%(infile)
        return 1

    hhalignHitList = ReadHHAlignResult(infile)
    numHit = len(hhalignHitList)
    if numHit < 1:
        print >> sys.stderr, "No hit found for file %s. Ignore"%infile
        return 1
    elif numHit > 1:
        print >> sys.stderr, "More than 1 (%d) hit found for file %s."%(numHit, infile)
        return 1

#     for item in hhalignHitList[0]:
#         print item, hhalignHitList[0][item]
    hit = hhalignHitList[0]

    if coverage_threshold >= 0.0:
        try:
            if hit['query_length'] >= hit['template_length']:
                coverage_of_shorter_seq = myfunc.FloatDivision(
                        len(hit['template_alignseq'].replace("-", "")),
                        hit['template_length'])
            else:
                coverage_of_shorter_seq = myfunc.FloatDivision(
                        len(hit['query_alignseq'].replace("-", "")),
                        hit['query_length'])
        except KeyError:
            print >> sys.stderr, "bad hit for file %s"%(infile)
            return 1
        if coverage_of_shorter_seq < coverage_threshold:
            print >> sys.stderr, "coverage (%.3f) < %g for %s. Ignore" %(coverage_of_shorter_seq,
                    coverage_threshold, infile)
            return 1

    if hit['evalue'] > evalue_threshold:
        print >> sys.stderr, "evalue (%g) > %g for %s. Ignore" %(hit['evalue'],
                evalue_threshold, infile)
        return 1


    query_rawseq = hdl_seq.GetRecord(hit['query_seqid'])
    if query_rawseq == None:
        return 1
    hit_rawseq = hdl_seq.GetRecord(hit['hit_seqid'])
    if hit_rawseq == None:
        return 1
    (hit_seqid, hit_annotation, hit_seq) = myfunc.ExtractFromSeqWithAnno(hit_rawseq)
    (query_seqid, query_annotation, query_seq) = myfunc.ExtractFromSeqWithAnno(query_rawseq)
    (hit_unaligned_head, hit_unaligned_tail) = GetUnalignedHeadTail(hit_seq, hit['pos_template_begin'], hit['pos_template_end'])
    (query_unaligned_head, query_unaligned_tail) = GetUnalignedHeadTail(query_seq, hit['pos_query_begin'], hit['pos_query_end'])
    (hit_unaligned_head, query_unaligned_head) = FillUnalignedGapForHead(hit_unaligned_head, query_unaligned_head)
    (hit_unaligned_tail, query_unaligned_tail) = FillUnalignedGapForTail(hit_unaligned_tail, query_unaligned_tail)

# output pairaln
    softmargin = 5
    if hit['pos_query_begin'] <= softmargin or hit['pos_template_begin'] <= softmargin:
        isHeadUnaligned = False
        query_unaligned_head = query_unaligned_head.upper()
        hit_unaligned_head = hit_unaligned_head.upper()
    else:
        isHeadUnaligned = True
        query_unaligned_head = query_unaligned_head.lower()
        hit_unaligned_head = hit_unaligned_head.lower()

    if (hit['pos_query_end'] >= hit['query_length']-softmargin or hit['pos_template_end']
            >= hit['template_length'] - softmargin):
        isTailUnaligned = False
        query_unaligned_tail = query_unaligned_tail.upper()
        hit_unaligned_tail = hit_unaligned_tail.upper()
    else:
        isTailUnaligned =  True
        query_unaligned_tail = query_unaligned_tail.lower()
        hit_unaligned_tail = hit_unaligned_tail.lower()

    complete_query_alignseq = "%s%s%s"%(query_unaligned_head, hit['query_alignseq'].upper(), query_unaligned_tail)
    complete_tempalte_alignseq = "%s%s%s"%(hit_unaligned_head, hit['template_alignseq'].upper(), hit_unaligned_tail)
    #print hit['query_alignseq']
    #print hit['template_alignseq']
    if fpout != None:
        fpout.write(">%s\n"%(hit['query_description']))
        fpout.write("%s\n"%complete_query_alignseq)
        fpout.write(">%s\n"%(hit['hit_description']))
        fpout.write("%s\n"%complete_tempalte_alignseq)

# output stat
    if fpout_stat != None:
        pos_query = "%d-%d"%(hit['pos_query_begin'], hit['pos_query_end'])
        pos_template = "%d-%d"%(hit['pos_template_begin'], hit['pos_template_end'])

        fpout_stat.write("%-8s %-8s %7g %8.3f %6.1f %6.1f %6d %9s %4d %9s %4d\n"%(
            hit['query_seqid'],
            hit['hit_seqid'],
            hit['evalue'],
            coverage_of_shorter_seq,
            hit['identity'],
            hit['prob'],
            hit['num_align_col'],
            pos_query,
            hit['query_length'],
            pos_template,
            hit['template_length'],
            ))

# output tableinfo
    if fpout_tableinfo != None:
        isLocalAlignment = True
        rd = lcmp.GetAlignmentFactorFromPairAlignment(hit['query_alignseq'],
                hit['template_alignseq'], isLocalAlignment)
        # rd = lcmp.GetAlignmentFactorFromPairAlignment(complete_query_alignseq, complete_tempalte_alignseq, isLocalAlignment)
        fpout_tableinfo.write("%-16s %-15s %6.1f %6.1f %9d %6d %6d %9.1f %6d %6d %6d %6.1f %6.1f\n"% (
            hit['query_seqid'],
            hit['hit_seqid'],
            rd['seqidt0'], hit['similarity']*100,
            rd['alnLength'],
            rd['seqLength1'], rd['seqLength2'],
            hit['score'],
            rd['numIDT'], -1, rd['numGap'],
            rd['seqidt1'], rd['seqidt2']))

#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    outfile_tableinfo = ""
    outfile_stat = ""
    fileListFile = ""
    fileList = []
    seqdb = ""
    evalue_threshold = 1e-3
    coverage_threshold = 0.0


    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            fileList.append(argv[i])
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-o", "--o", "-outfile"]:
                (outfile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-ot", "--ot"]:
                (outfile_tableinfo, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-os", "--os"]:
                (outfile_stat, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-evalue", "--evalue"]:
                (evalue_threshold, i) = myfunc.my_getopt_float(argv, i)
            elif argv[i] in ["-coverage", "--coverage"]:
                (coverage_threshold, i) = myfunc.my_getopt_float(argv, i)
            elif argv[i] in ["-seqdb", "--seqdb"]:
                (seqdb, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-l", "--l"] :
                (fileListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            fileList.append(argv[i])
            i += 1

    if fileListFile != "":
        fileList += myfunc.ReadIDList(fileListFile)
    if len(fileList) < 1:
        print >> sys.stderr, "No input set. exit"
        return 1
    if seqdb == "":
        print >> sys.stderr, "Seqdb not set. exit"
        return 1

    hdl_seq = myfunc.MyDB(seqdb)
    if hdl_seq.failure:
        return 1

    fpout = sys.stdout
    fpout_tableinfo = None
    fpout_stat = None
    if outfile != "":
        fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    if outfile_tableinfo != "":
        fpout_tableinfo = myfunc.myopen(outfile_tableinfo, None, "w", False)  
    if outfile_stat != "":
        fpout_stat = myfunc.myopen(outfile_stat, None, "w", False)

    if fpout_stat != None:
        fpout_stat.write("%-8s %-8s %7s %8s %6s %6s %6s %9s %4s %9s %4s\n"%(
           "#ID1",
           "ID2",
           "Evalue",
           "Coverage",
           "IDT",
           "Prob",
           "AlnCol",
           "PosQuery",
           "LenQ",
           "PosTemp",
           "LenT",
            ))
    if fpout_tableinfo != None:
        fpout_tableinfo.write("#%-15s %-15s %6s %6s %9s %6s %6s %9s %6s %6s %6s %6s %6s\n" % (
        "Seq1","Seq2", "IDT0", "SIM0", "AlnLength", "Len1","Len2",
        "Score","N_IDT", "N_SIM", "N_GAP", "IDT1", "IDT2"))

    for infile in fileList:
        HHAlign2Pairaln(infile, evalue_threshold, coverage_threshold, hdl_seq,
                fpout,fpout_tableinfo, fpout_stat) 


    if outfile != "":
        myfunc.myclose(fpout)
    if outfile_tableinfo != "":
        myfunc.myclose(fpout_tableinfo)
    if outfile_stat != "":
        myfunc.myclose(fpout_stat)
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))
