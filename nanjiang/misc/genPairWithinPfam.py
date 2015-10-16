#!/usr/bin/env python
# Description:

# ChangeLog 2012-06-10 
#    TMpro list can be added, so that pairs are generated only for TM proteins
import os
import sys
import myfunc
import random

DATADIR3 = os.environ['DATADIR3']

progname =  os.path.basename(sys.argv[0])
usage="""
Usage: %s  [ID [ID...]] [-l idListFile] [-outpath DIR]

DESCRIPTION: 
  Generate pairs (by random selection) within pfam

OPTIONS:
  -l   LISTFILE  Set IDListFile
  -o    OUTFILE  Output the result to FILE
  -m 0|1         Method for generating pairs, (default: 0)
                 0: random selection of pairs, one protein might be used >=2 times
                 1: random selection of pairs, one protein can be used at most once
                 2. all-to-all pairs
  -mapfile FILE  Set pfamid2seqid mapfile

  -maxseq   INT  Maximum number of sequence to use within each fam, 
                 (default: 200) this is for method 0
  -maxpair  INT  Maximum number of pairs for each fam, 
                 (default: 300) this is for method 1
  -seed     INT  Set random seed, default is set by time
  -pdbtosp FILE  PDB to swissprot id (uniprot ac) maplist 
  -onlypdb       Output only pairs from the PDB

  -outwithfamid FILE
        Output the pairlist with famid

  -outwithpdb FILE
        Output the pairlist with pdb id

  -outfam2seqmap FILE
        output the new fam2seq map from which the pairlist was generated

  -tmprolist, -restrictlist FILE
        Set TMpro idlist file, restrict pair slection to only
        those in the list

  -h, --help     Print this help message and exit

Created 2012-05-30, updated 2014-06-02, Nanjiang Shu 

Examples:
    %s  -l idlist.txt -mapfile t1.pfamid2seqid -tmprolist t1.tmproidlist -o pairlist.txt -outwithfamid pairlist.withfamid.txt
"""%(progname, progname)


def PrintHelp():
    print usage

def GeneratePairWithinFam_m_0(idList, idMapDict, restrictIDSet,#{{{
        maxseq_for_fam, rand_seed, fpout, fpout_withfamid):
# select maxseq_for_fam sequences to generate pairs to balance large families
# and small famlies
    for famid in idList:
        try:
            seqidlist = idMapDict[famid]
        except KeyError:
            print >> sys.stderr, "famid %s not found in idMapDict"%(famid)
            continue

        if restrictIDSet != set([]):
            seqidlist = list(set(seqidlist) & restrictIDSet)
        numseqid = len(seqidlist)
        if numseqid > maxseq_for_fam:
            seqidlist = random.sample(seqidlist, maxseq_for_fam)
        if len(seqidlist) < 2:
            msg =  "numseq for famid %s (%d) < 2, ignore"
            print >> sys.stderr, msg%(famid, len(seqidlist))
            continue
        max_numpair = int(numseqid*numseqid/2.0/5)
        pairlist = myfunc.GenerateRandomPair(len(seqidlist), max_numpair, rand_seed)
        for pair in pairlist:
            if fpout != None:
                print >> fpout, "%s %s" % (seqidlist[pair[0]],
                        seqidlist[pair[1]])
            if fpout_withfamid != None:
                print >> fpout_withfamid, "%s %s %s" % (seqidlist[pair[0]],
                        seqidlist[pair[1]], famid)
    return 0
#}}}

def GeneratePairWithinFam_m_1(idList, idMapDict, restrictIDSet,#{{{
        maxpair_for_fam, rand_seed, fpout, fpout_withfamid,
        fpout_fam2seqmap):
    for famid in idList:
        try:
            seqidlist = idMapDict[famid]
        except KeyError:
            print >> sys.stderr, "famid %s not found in idMapDict"%(famid)
            continue

        if restrictIDSet != set([]):
            seqidlist = list(set(seqidlist) & restrictIDSet)

        if fpout_fam2seqmap != None:
            fpout_fam2seqmap.write("%s %d"%(famid, len(seqidlist)))
            for seqid in seqidlist:
                fpout_fam2seqmap.write(" %s"%(seqid))
            fpout_fam2seqmap.write("\n")


        pairlist = myfunc.GenerateRandomPair_no_repeat_use(len(seqidlist),
                maxpair_for_fam, rand_seed)
        for pair in pairlist:
            if fpout != None:
                print >> fpout, "%s %s" % (seqidlist[pair[0]],
                        seqidlist[pair[1]])
            if fpout_withfamid != None:
                print >> fpout_withfamid, "%s %s %s" % (seqidlist[pair[0]],
                        seqidlist[pair[1]], famid)
    return 0#}}}
def GeneratePairWithinFam_m_2(idList, idMapDict, restrictIDSet,#{{{
        fpout, fpout_withfamid, fpout_withpdb):
    """
    Input:
        idList     idlist for families
    """
    for famid in idList:
        try:
            seqidlist = idMapDict[famid]
        except KeyError:
            print >> sys.stderr, "famid %s not found in idMapDict"%(famid)
            continue

        if restrictIDSet != set([]):
            seqidlist = list(set(seqidlist) & restrictIDSet)
        if g_params['isOnlyPDB']:
            seqidlist = list(set(seqidlist) & g_params['uniprotidlist_with_pdb'])
# generate all-to-all pairs
        numseq = len(seqidlist)
        for i in xrange(numseq):
            for j in xrange(i+1, numseq):
                pair = (i,j)
                if fpout != None:
                    print >> fpout, "%s %s" % (seqidlist[pair[0]],
                            seqidlist[pair[1]])
                if fpout_withfamid != None:
                    print >> fpout_withfamid, "%s %s %s" % (seqidlist[pair[0]],
                            seqidlist[pair[1]], famid)
                if g_params['isOnlyPDB'] and fpout_withpdb != None:
                    print >> fpout_withpdb, "%s %s %s" % (seqidlist[pair[0]], seqidlist[pair[1]], famid),\
                            g_params['uniprot2pdbMap'][seqidlist[pair[0]]].keys(),\
                            g_params['uniprot2pdbMap'][seqidlist[pair[1]]].keys()


    return 0
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    outfile_with_famid = ""
    outfile_with_pdb = ""
    outfile_fam2seqmap = ""
    idListFile = ""
    mapfile = "%s%s%s"%(DATADIR3, os.sep, "wk/MPTopo/pfamAna_refpro/pfammap_from_uniprot/refpro20120604-celluar.selmaxlength-m1.nr100.filter.fragmented.clanid2seqid")
    restrictIDListFile = ""
    idList = []
    maxseq_for_fam = 200
    maxpair_for_fam = 300
    method = 0
    rand_seed = None
    pdbtospFile = ""
    isOnlyPDB = False

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            idList.append(argv[i])
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-o", "--o", "-outfile", "--outfile"]:
                outfile, i = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-outwithfamid","--outwithfamid"]:
                outfile_with_famid,i = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-outfam2seqmap","--outfam2seqmap"]:
                outfile_fam2seqmap,i = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-outwithpdb","--outwithpdb"]:
                outfile_with_pdb,i = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-tmprolist","--tmprolist", "-restrictlist",
                    "--restrictlist"]:
                restrictIDListFile,i = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-mapfile", "--mapfile"]:
                mapfile,i = myfunc.my_getopt_str(argv,i)
            elif (argv[i] in ["-pdbtosp", "--pdbtosp"]):
                pdbtospFile, i = myfunc.my_getopt_str(argv, i)
            elif sys.argv[i]  in [ "-seed" , "--seed"]:
                rand_seed, i = myfunc.my_getopt_int(argv,i)
            elif argv[i] in ["-l", "--l"] :
                idListFile, i = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-maxseq", "--maxseq"]:
                maxseq_for_fam, i = myfunc.my_getopt_int(argv,i)
            elif argv[i] in ["-maxpair", "--maxpair"]:
                maxpair_for_fam, i = myfunc.my_getopt_int(argv,i)
            elif argv[i] in ["-m", "--m", "-method", "--method"]:
                method, i = myfunc.my_getopt_int(argv,i)
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True; i += 1
            elif argv[i] in ["-onlypdb", "--onlypdb"]:
                g_params['isOnlyPDB'] = True; i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            idList.append(argv[i])
            i += 1

    if os.path.exists(idListFile):
        idList += myfunc.ReadIDList(idListFile)

    if len(idList) < 1:
        print >> sys.stderr, "no ID set. exit"
        return 1
    if myfunc.checkfile(mapfile, "idMapFile") != 0:
        return 1

    idMapDict = myfunc.ReadFam2SeqidMap(mapfile)

# Read in pdbtosp map
    if pdbtospFile != "":
        (pdb2uniprotMap, uniprot2pdbMap) =\
                myfunc.ReadPDBTOSP(pdbtospFile)
        g_params['uniprotidlist_with_pdb'] = set(uniprot2pdbMap.keys())
        g_params['uniprot2pdbMap'] = uniprot2pdbMap

    if g_params['isOnlyPDB'] == True:
        if pdbtospFile == "":
            print >> sys.stderr, "onlypdb is true but pdbtospFile is not set. exit."
            return 1
        elif g_params['uniprotidlist_with_pdb'] == set([]):
            print >> sys.stderr, "onlypdb is true but uniprotidlist_with_pdb is empty. exit."
            return 1

    restrictIDSet = set([])
    if restrictIDListFile != "":
        restrictIDSet = set(myfunc.ReadIDList(restrictIDListFile))

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    fpout_withfamid = myfunc.myopen(outfile_with_famid, None, "w", False)
    fpout_withpdb = myfunc.myopen(outfile_with_pdb, None, "w", False)
    fpout_fam2seqmap = myfunc.myopen(outfile_fam2seqmap, None, "w", False)


    if method == 0:
        GeneratePairWithinFam_m_0(idList, idMapDict, restrictIDSet,
                maxseq_for_fam, rand_seed, fpout, fpout_withfamid)
    elif method == 1:
        GeneratePairWithinFam_m_1(idList, idMapDict, restrictIDSet,
                maxpair_for_fam, rand_seed, fpout, fpout_withfamid,
                fpout_fam2seqmap)
    elif method == 2: #all to all
        GeneratePairWithinFam_m_2(idList, idMapDict, restrictIDSet,
                fpout, fpout_withfamid, fpout_withpdb)




    myfunc.myclose(fpout)
    myfunc.myclose(fpout_withfamid)
    myfunc.myclose(fpout_withpdb)
    myfunc.myclose(fpout_fam2seqmap)
    return 0
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['isOnlyPDB'] = False
    g_params['uniprotidlist_with_pdb'] = set([])
    g_params['uniprot2pdbMap'] = {}
    return g_params
#}}}
if __name__ == '__main__' :
    rundir = os.path.dirname(sys.argv[0])
    binpath = rundir
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))
