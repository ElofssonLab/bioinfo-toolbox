#!/usr/bin/env python

import sys,re,os
import comptopo as ct
import myfunc
import time
import libtopologycmp as lcmp
import subprocess

def PrintFuncName():
    print  ("this function name is: %s"% sys._getframe().f_code.co_name)
def ExtractHit(line):#{{{
    hit = {}
    strs = line.split()
    if len(strs) >= 11:
        evalue = float(strs[4])
        posQuery = strs[9].split("-")
        posQuery = [int(x) for x in posQuery]
        posTemplate = strs[10].split("-")
        posTemplate = [int(x) for x in posTemplate]
        hit['evalue'] = evalue
        hit['posQuery'] = posQuery
        hit['posTemplate'] = posTemplate
    return hit
#}}}
def IsDuplicatedByHHSearch(hhrfile):#{{{
    try:
        fpin = open(hhrfile,"r")
        lines = fpin.readlines()
        fpin.close()
        hitList = []
        numLine = len(lines)
        i = 0
        while i < numLine:
            line = lines[i]
            if line.find(" No Hit") == 0:
                j = 1
                while i+j < numLine and lines[i+j] != "":
                    hit = ExtractHit(lines[i+j])
                    if hit != {}:
                        hitList.append(hit)
                    else:
                        break
                    j += 1
                break
            else:
                i += 1
        if len(hitList) < 2:
            return False
        else:
            sortedHitList = sorted(hitList, key=lambda x:x['evalue'], reverse=False)
            hit1 = hitList[0]
            hit2 = hitList[1]
            if hit2['evalue'] > 1e-3:
                return False
            else:
                (b1, e1) = hit1['posTemplate']
                (b2, e2) = hit2['posTemplate']
                overlap = max(0, myfunc.coverage(b1, e1, b2, e2))
                if overlap / float(e1-b1) < 0.2 and overlap / float(e2-b2) < 0.2:
                    return True
                else:
                    return False
    except IOError:
        print >> sys.stderr, "Failed to read hhrfile %s"%hhrfile
        return False
#}}}
def GetAlignedRegion(pairaln):
    try:
        lengthAln = len(pairaln[0])
        b = [0]*2
        e = [lengthAln]*2
        for i in xrange(2):
            while pairaln[i][b[i]] == "-":
                b[i] += 1
            while pairaln[i][e[i]-1] == "-":
                e[i] -= 1
        return (max(b),min(e))
    except IndexError:
        return (0, 0)
def main():  #{{{
    if 0:#{{{
        strTop1 = "---MMMM-----i-i-i---MMM----MMMM-ooo"
        strTop2 = "----MMMM-----i-ii-----MMM---MMM--oo"
        strProtein1 = "id1"
        strProtein2 = "id2"
        fpLog=sys.stdout
        class_gapless, num1_gapless, num2_gapless = ct.CompareToposGaplesslyNew(strTop1, strTop2, strProtein1, strProtein2, fpLog)
        # Note: calling the int, float, string will not change their original value
        # calling the dict, list will change their original value
        print "strTop1:",strTop1
        print "strTop2:",strTop2
#}}}
    if 0:#{{{
        PrintFuncName()
        print  ("this file name is: %s"% __file__)
#}}}
    if 0:#{{{
        # filename="/nanjiang/data/blastdb/uniprot_KW181_idt50.fasta"
        filename = sys.argv[1]
        print filename
        fp = open(filename,"r")
        lines=fp.readlines()
        fp.close()
#}}}
    if 0:#{{{
        # filename="/nanjiang/data/blastdb/uniprot_KW181_idt50.fasta"
        filename = sys.argv[1]
        print filename
        BLOCK_SIZE = 100000
        fp = open(filename,"r")
        buff = fp.read(BLOCK_SIZE)
        while buff:
            buff = fp.read(BLOCK_SIZE)
        fp.close()
#}}}
    if 0:#{{{
        # filename="/nanjiang/data/blastdb/uniprot_KW181_idt50.fasta"
        filename = sys.argv[1]
        print filename
        fp = open(filename,"r")
        line = fp.readline()
        while line:
            line = fp.readline()
        fp.close();#}}}
    if 0:#{{{
        try:
            BLOCK_SIZE = 100000
            infile = sys.argv[1]
            fpin = open(infile, 'rb')
            unprocessedBuffer=""
            isEOFreached = False
            while 1:
                buff = fpin.read(BLOCK_SIZE)
                if len(buff) < BLOCK_SIZE:
                    isEOFreached=True
                buff = unprocessedBuffer + buff
                recordList = []
                unprocessedBuffer = myfunc.ReadFastaFromBuffer(buff,recordList, isEOFreached)
                if len(recordList) > 0: 
                    for record in recordList:
                        sys.stdout.write(">%s\n"%record[1])
                        sys.stdout.write("%s\n"%record[2])
                if isEOFreached == True:
                    break
            fpin.close()
        except IOError:
            raise;   #}}}
    if 0:#{{{
        try:
            infile = sys.argv[1]
            (annoList, seqList)  = myfunc.ReadFasta_without_id(infile)
            for i in xrange(len(seqList)):
                sys.stdout.write(">%s\n"%annoList[i])
                sys.stdout.write("%s\n"%seqList[i])
        except IOError:
            raise;   #}}}
    if 0:#{{{
        hhrfile = "hhsearch/A1RZ92-Q74DY9.hhr"
        if IsDuplicatedByHHSearch(hhrfile):
            print "yes"


#}}}
    if 0:#{{{
        import pairlistwithfamid2pairaln_by_msa
        seq1 = "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------MLSSTATTMLRAGVSRSSGALQPMLLRSAACPCSPFSMNTKLSQPTSV-----RPLSTSPSALVLRFRAQQQAQLAQQQLRRASSSSSSSSSSTRPRSDAELDANAAEAAAAAQSAAHAGEPVLDWNTFFKLRKTRRRVQLAFSVIMTLITSGAGGAVLSTGVADAMVAQVPLEPMFAVGLMTASFGALGWLMGPAMGGMVFNALKSKYRGQMEIKEGQFFARIKKHRVDPSASSMGNPVPDFYGEKISSVAGYRQWLKDQRAFNKKRTTFV"
        seq2 = "MDILLAVLEQGFIFSIVCFGVYITYKILDFPDLSVDGTFPLGAAVAAAFLVKGYSPVLSSLAALVAGAIAGGITGILHVKFKITNLLSGILVMVGLYSINLRIMGKSNIPLFNKIHLFSDTMNPIIIITVFLLICKITLDLFLKTKAGFILKATGDNEQLVLSLGVNKDLVKIMGLMLSNALVALGGALMAQYQGFSDVGMGTGIVVMGLASVIIGESLFGRIKALNATTRVLLGALVYKLSVSI---ALTVGLAP-------TDLKLVTAIIVVIALSLNKNPLKIITKQKTKEGGIL------NASNTKSAQSVQ-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
        seq1 = "---------------------------------------------------------------------------------------------------------------------------------------MALSSLFFTASALLLMFLAFLGGARNSNPLDRIYWLEAATGNIPGAPALSRWTYWNLCAVNSEGHNECGKSYPDYPFDPPSHRNFNTHVNIPAAFIGTRHYFLTSRFMFPFHIIALFFATCSLLTGFLAMCTRIGNWVSAFSAYFALTFQTITTCLMTAVYVQGRDKFNNNGQSSHLGVKAFAFMWTSVALLFLSCVIYCMGGAVGRKDGGYSGREQRRRGFFNSHRSGSLRSNKETAP"
        seq2 = "MRKIAAIGGIVFISFILTIVAMFTKLWISWSIGKFSYGIGIVPYHSNSAGWFTAASWMVFISFGLFIPLILVVLFTAYKVHHDGCCHSIRHCFNSICLICSIIAVLEIIAFVLMAVNASRYVKGASISEKKSLLQLGSSAYLDLVSAILIIVATVLSGHASHHDCH----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
        alignFactor = pairlistwithfamid2pairaln_by_msa.GetAlignmentFactorFromPairAlignment(seq1,seq2)
        print alignFactor
#}}}
    if 0:#{{{
        try:
            dbname = sys.argv[1]
            print dbname 
            from myfunc import MyDB
            cls = MyDB(dbname)
#            print cls.idList
            record =  cls.GetRecord("A0FGX9")
            if record:
                print record
        #             for rd in  cls.GetAllRecord():
        #                 print rd
                (seqid, anno, seq) = myfunc.ExtractFromSeqWithAnno(record)
                print (seqid, anno, seq)
        except IndexError:
            pass

#}}}
    if 0:#{{{
        import my_extractdb
        #miniking my_extractdb.py see which one is faster
        try:
            dbname = sys.argv[1]
            idlistfile = sys.argv[2]
            cls = myfunc.MyDB(dbname)
            if cls.failure:
                print >> sys.stderr, "MyDB init failed"
            else:
                idlist = open(idlistfile, "r").read().split("\n")
                fpout = sys.stdout
                for seqid in idlist:
                    if seqid:
                        record =  cls.GetRecord(seqid)
                        fpout.write(record)
            #             for rd in  cls.GetAllRecord():
            #                 print rd
#                (seqid, anno, seq) = myfunc.ExtractFromSeqWithAnno(record)
#                print (seqid, anno, seq)
        except IndexError:
            print "error"
            pass
#}}}
    if 0:#{{{ #test ReadLineByBlock
        try: 
            infile = sys.argv[1]
            from myfunc import ReadLineByBlock
            cls = ReadLineByBlock(infile)
            lines = cls.readlines()
            while lines != None:
                for line in lines:
                    print line
                lines = cls.readlines()

        except IndexError:
            pass
#}}}
    if 0:#{{{ #test speed of ReadLineByBlock
# ReadLineByBlock is about 3 times fater than file.readline()
        try:
            from myfunc import ReadLineByBlock
            infile = sys.argv[1]

            start = time.time()
            hdl = ReadLineByBlock(infile)
            lines = hdl.readlines()
            while lines != None:
                lines = hdl.readlines()
            hdl.close()
            end = time.time()
            msg =  "Reading %s by ReadLineByBlock costs %.3fs seconds"
            print  msg%(infile, (end-start))

            start = time.time()
            hdl  = open(infile, "r")
            line = hdl.readline()
            while line:
                line = hdl.readline()
            hdl.close()
            end = time.time()
            msg =  "Reading %s by readline() costs %.3fs seconds"
            print  msg%(infile, (end-start))

        except IndexError:
            pass
#}}}
    if 0:#{{{ #test readline
        try: 
            infile = sys.argv[1]
            fp = open(infile, "r")
            line = fp.readline()
            while line:
                print line
                line = fp.readline()
            fp.close()
        except IndexError:
            pass
#}}}
    if 0:#{{{ #test the speed of GetFirstWord
        try:
            nloop = int(sys.argv[1])
            string = "kjdafk jasdfj j"
            #string = "askdf askdf "
#            string = "kajsdfasdfsdfjakasjdfka"
#            string = "kajsdfasdf,sdfjakasjdfka"
            delimiter = " \t\r,.\n"
            delimiter = " "
            for i in xrange(nloop):
                #firstword = myfunc.GetFirstWord(string, delimiter)
                #firstword = string.split()[0]
                #firstword = string.partition(" ")[0]
                firstword = myfunc.GetFirstWord(string)
                #pass
                #print firstword
        except (IndexError, ValueError):
            pass
#}}}
    if 0:#{{{ # read seq by SeqIO
        from Bio import SeqIO
        try:
            seqfile = sys.argv[1]
            # 1. SeqIO ####################
            start = time.time()
            handle = open(seqfile, "rU")
            cnt = 0
            for record in SeqIO.parse(handle, "fasta") :
                cnt += 1
            handle.close()
            end = time.time()
            msg =  "Reading %d sequences by SeqIO costs %.3fs seconds"
            print  msg%(cnt, (end-start))

            # 2. ReadFasta ####################
            start = time.time()
            seqfile = sys.argv[1]
            (idList, annoList, seqList) = myfunc.ReadFasta(seqfile)
            end = time.time()
            msg =  "Reading %d sequences by ReadFasta costs %.3fs seconds"
            print  msg%(len(idList), (end-start))

            # 3. ReadFasta from buffer
            BLOCK_SIZE = 100000
            start = time.time()
            cnt = 0
            fpin = open(seqfile, 'rb')
            unprocessedBuffer=""
            isEOFreached = False
            while 1:
                buff = fpin.read(BLOCK_SIZE)
                if len(buff) < BLOCK_SIZE:
                    isEOFreached=True
                buff = unprocessedBuffer + buff
                recordList = []
                unprocessedBuffer = myfunc.ReadFastaFromBuffer(buff,recordList, isEOFreached)
                cnt += len(recordList)
                if isEOFreached == True:
                    break
            fpin.close()
            end = time.time()
            msg =  "Reading %d sequences by ReadFastaFromBuffer costs %.3fs seconds"
            print  msg%(cnt, (end-start))

            # 4. ReadFastaByBlock ####################
            start = time.time()
            seqfile = sys.argv[1]
            hdl = myfunc.ReadFastaByBlock(seqfile, 0, 0)
            if hdl.failure:
                print >> sys.stderr, "Failed to init ReadFastaByBlock"
                return 1
            recordList = hdl.readseq()
            cnt = 0
            while recordList != None:
                cnt += len(recordList)
#                 for rd in recordList:
#                     print ">%s"%rd.description
#                     print rd.seq
                recordList = hdl.readseq()
            hdl.close()
            end = time.time()
            msg =  "Reading %d sequences by ReadFastaByBlock costs %.3fs seconds"
            print  msg%(cnt, (end-start))
        except (IndexError, ValueError):
            pass
#}}}
    if 0:#{{{ #test RemoveUnnecessaryGap
        try:
            infile = sys.argv[1]
            start = time.time()
            (idList, seqList) = myfunc.ReadFasta_without_annotation(infile)
            seqList = lcmp.RemoveUnnecessaryGap_old(seqList)
            end = time.time()
            msg =  "Run RemoveUnnecessaryGap_old for %s costs %.3fs seconds"
            print >> sys.stderr,  msg%(infile, (end-start))
            for seq in seqList:
                print seq

            start = time.time()
            (idList, seqList) = myfunc.ReadFasta_without_annotation(infile)

            seqList = lcmp.RemoveUnnecessaryGap(seqList)
            end = time.time()
            msg =  "Run RemoveUnnecessaryGap for %s costs %.3fs seconds"
            print >> sys.stderr,  msg%(infile, (end-start))
            for seq in seqList:
                print seq

        except IndexError:
            pass
#}}}
    if 0:#{{{ #test ReadMPAByBlock
        try: 
            infile = sys.argv[1]
            hdl = myfunc.ReadMPAByBlock(infile)
            if hdl.failure:
                return 
            recordList = hdl.readseq()
            while recordList != None:
                for rd in recordList:
                    #print rd.seqid
                    print ">%s"%(rd.description)
                    print "%s"%(myfunc.mpa2seq(rd.mpa))
                recordList = hdl.readseq()
            hdl.close()
        except IndexError:
            pass
#}}}
    if 0:#{{{
        try:
            dbname = sys.argv[1]
            print dbname 
            from myfunc import MyDB
            cls = MyDB(dbname)
#            print cls.idList
            record =  cls.GetRecord("A0FGX9")
            if record:
                print record
        #             for rd in  cls.GetAllRecord():
        #                 print rd
                (seqid, anno, seq) = myfunc.ExtractFromSeqWithAnno(record)
                print (seqid, anno, seq)
        except IndexError:
            pass

#}}}
    if 0:#{{{ #test subprocess
        import glob
        #invoke shell explicitly, not very good, may have security problems
        subprocess.call("seq 10", shell=True)
        subprocess.call("echo wait for 2 seconds...; sleep 2", shell=True)
        subprocess.call("ls topo*.py", shell=True)
    if 1:#{{{ #test subprocess
        import glob
        #invoke shell implicitly, recommended way
        subprocess.call(["seq", "10"], shell=False)
        subprocess.call(["echo", "wait for 1 seconds..."])
        subprocess.call(["sleep", "1"])
        try:
            print subprocess.check_call(["ls", "topo*.py"]) #This will not work
        except subprocess.CalledProcessError, e:
            print "error message:", e
        subprocess.call(["ls"]+glob.glob("topo*.py"))
#}}}
#}}}

if __name__ == "__main__":
    sys.exit(main())
