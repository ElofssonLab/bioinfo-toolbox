#!/usr/bin/env python
# Description:
import os
import sys
import myfunc
usage="""
usage:  drawseqmsa_topo.py  seqmsafile_in_fasta_format
                      
Description:
    Given sequence MSA file, draw topomsa image

Options:
  -outpath DIR    Set ouput path
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2012-03-21, updated 2012-03-21, Nanjiang Shu 
"""

def PrintHelp():
    print usage;

def DrawSeqMSA(seqmsafile, outpath):
    print "Remove gaps from sequence"
    (idList, annotationList, seqList) = myfunc.ReadFasta(seqmsafile)
    rootname=os.path.basename(os.path.splitext(seqmsafile)[0]);
    basename = os.path.basename(seqmsafile)
    seqfile = outpath + os.sep + rootname + '.fa'
    fpout = open(seqfile,"w")
    for i in xrange(len(idList)):
        fpout.write(">%s\n"%annotationList[i])
        fpout.write("%s\n"%seqList[i].replace("-","").replace(".",""))
    fpout.close()

    print "Predicting topologies..."
    scampi_exe = "%s/mySCAMPI_run.pl" %g_params['newscampiscriptpath']
    scampi_dir = g_params['scampi_dir']
    modhmm_bin = g_params['modhmm_bin']
    cmd = "%s %s --scampipath %s --modhmmpath %s --outpath %s" % (
            scampi_exe, seqfile, scampi_dir, modhmm_bin, outpath)
    os.system(cmd)
    os.system("rm -f %s/*.res"%outpath)

    print "Get topomsa"
    binpath = g_params['binpath']
    topofile = outpath + os.sep + rootname + '.fa.topo'
    topomsafile = outpath + os.sep + rootname + '.topomsa.fa'
    cmd = "%s/matchMSAtopo -msa %s -topo %s -o %s"%(
            binpath, seqmsafile, topofile, topomsafile)
    os.system(cmd)

    print "Draw topomsa"
    cmd = "python %s/drawMSATopo.py %s -text y -outpath %s -aaseq %s"%(
            binpath, topomsafile, outpath, seqfile)
    os.system(cmd)


def main(g_params):#{{{
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    isQuiet=False;
    outpath="./";
    seqmsafile_in_fasta_format = ""

    i = 1;
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            seqmsafile_in_fasta_format = sys.argv[i]
            isNonOptionArg=False;
            i += 1;
        elif sys.argv[i] == "--":
            isNonOptionArg=True;
            i += 1;
        elif sys.argv[i][0] == "-":
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp();
                return 1
            elif sys.argv[i] in ["-outpath", "--outpath"]:
                outpath=sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-q":
                isQuiet=True;
                i += 1;
            else:
                print >> sys.stderr, "Error! Wrong argument:", sys.argv[i];
                return 1
        else:
            seqmsafile_in_fasta_format = sys.argv[i]
            i += 1
    if not os.path.exists(outpath):
        os.system("mkdir %s"%outpath)
    g_params['outpath'] = outpath
    rundir = os.path.dirname(sys.argv[0])
    g_params['newscampiscriptpath'] = "%s/../../../program/newscampiscript"%rundir
    g_params['scampi_dir'] = "%s/../../../share/scampi"%rundir 
    g_params['modhmm_bin'] = "%s/../../../share/modhmm/bin"%rundir 
    g_params['binpath'] = rundir

    if (seqmsafile_in_fasta_format != "" and
            os.path.exists(seqmsafile_in_fasta_format)):
        DrawSeqMSA(seqmsafile_in_fasta_format, outpath)
#}}}

if __name__ == '__main__' :
    g_params = {}
    sys.exit(main(g_params));
